"""Thin client for the NCI GDC API: paginated file/case queries and streamed downloads."""

from __future__ import annotations

import logging
import re
import shutil
import tarfile
import tempfile
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd
import requests
from tqdm.auto import tqdm

log = logging.getLogger(__name__)

API_URL = "https://api.gdc.cancer.gov"
BETA_DATA_TYPE = "Methylation Beta Value"

FILE_FIELDS = [
    "file_id",
    "file_name",
    "file_size",
    "md5sum",
    "platform",
    "cases.case_id",
    "cases.submitter_id",
    "cases.project.project_id",
    "cases.primary_site",
    "cases.samples.sample_id",
    "cases.samples.sample_type",
    "cases.samples.tissue_type",
]

CASE_FIELDS = [
    "case_id",
    "submitter_id",
    "project.project_id",
    "demographic.vital_status",
    "demographic.days_to_death",
    "demographic.gender",
    "diagnoses.age_at_diagnosis",
    "diagnoses.days_to_last_follow_up",
    "diagnoses.ajcc_pathologic_stage",
    "diagnoses.tumor_grade",
    "diagnoses.primary_diagnosis",
    "diagnoses.diagnosis_is_primary_disease",
    "follow_ups.days_to_follow_up",
]

_PLATFORM_ALIASES = {
    "450": "450K",
    "27": "27K",
    "epic v2": "EPICv2",
    "epic": "EPIC",
}


def normalize_platform(platform: str | None) -> str | None:
    """Map GDC platform strings ('Illumina Human Methylation 450', ...) to 450K/EPIC/EPICv2/27K."""
    if platform is None:
        return None
    p = platform.lower()
    for key, alias in _PLATFORM_ALIASES.items():
        if key in p:
            return alias
    return platform


def _in(field: str, values: Iterable[str]) -> dict:
    return {"op": "in", "content": {"field": field, "value": list(values)}}


def _and(*clauses: dict | None) -> dict:
    return {"op": "and", "content": [c for c in clauses if c is not None]}


def paginate(
    endpoint: str,
    filters: dict,
    fields: Sequence[str],
    page_size: int = 500,
    session: requests.Session | None = None,
) -> list[dict]:
    """Return every hit for a GDC search, following pagination (the API caps a page at `size`)."""
    http = session or requests.Session()
    hits: list[dict] = []
    offset = 0
    while True:
        body = {
            "filters": filters,
            "fields": ",".join(fields),
            "format": "JSON",
            "size": page_size,
            "from": offset,
        }
        resp = http.post(f"{API_URL}/{endpoint}", json=body, timeout=120)
        resp.raise_for_status()
        data = resp.json()["data"]
        hits.extend(data["hits"])
        total = data["pagination"]["total"]
        offset += len(data["hits"])
        if offset >= total or not data["hits"]:
            return hits


def query_methylation_files(
    primary_site: str | Sequence[str] | None = None,
    project_ids: Sequence[str] | None = None,
    session: requests.Session | None = None,
) -> pd.DataFrame:
    """All methylation beta-value files for a primary site and/or project, one row per file.

    No sample-type, platform or per-patient filtering happens here; that is
    `methylspan.cohort.select_cohort`'s job, so the raw composition stays visible.
    """
    if primary_site is None and project_ids is None:
        raise ValueError("Give a primary_site, project_ids, or both")
    sites = [primary_site] if isinstance(primary_site, str) else primary_site
    filters = _and(
        _in("files.data_type", [BETA_DATA_TYPE]),
        _in("cases.primary_site", sites) if sites else None,
        _in("cases.project.project_id", project_ids) if project_ids else None,
    )
    hits = paginate("files", filters, FILE_FIELDS, session=session)

    rows = []
    for h in hits:
        cases = h.get("cases") or [{}]
        if len(cases) > 1:
            log.warning("File %s is linked to %d cases; using the first", h["file_id"], len(cases))
        case = cases[0]
        samples = case.get("samples") or [{}]
        rows.append(
            {
                "file_id": h["file_id"],
                "file_name": h["file_name"],
                "file_size": h.get("file_size"),
                "md5sum": h.get("md5sum"),
                "platform": normalize_platform(h.get("platform")),
                "case_id": case.get("case_id"),
                "submitter_id": case.get("submitter_id"),
                "project_id": (case.get("project") or {}).get("project_id"),
                "sample_id": samples[0].get("sample_id"),
                "sample_type": samples[0].get("sample_type"),
                "tissue_type": samples[0].get("tissue_type"),
            }
        )
    cols = ["file_id", "file_name", "file_size", "md5sum", "platform", "case_id",
            "submitter_id", "project_id", "sample_id", "sample_type", "tissue_type"]
    return pd.DataFrame(rows, columns=cols)


def query_cases(
    case_ids: Sequence[str],
    fields: Sequence[str] = CASE_FIELDS,
    batch_size: int = 500,
    session: requests.Session | None = None,
) -> list[dict]:
    """Raw clinical records for the given case UUIDs."""
    ids = list(dict.fromkeys(case_ids))
    hits: list[dict] = []
    for i in range(0, len(ids), batch_size):
        hits.extend(paginate("cases", _in("case_id", ids[i : i + batch_size]), fields, session=session))
    return hits


def latest_file_id(file_id: str, session: requests.Session | None = None) -> str | None:
    """Current released UUID for a possibly superseded GDC file id, via the /history endpoint."""
    http = session or requests.Session()
    resp = http.get(f"{API_URL}/history/{file_id}", timeout=60)
    resp.raise_for_status()
    released = [h for h in resp.json() if h.get("file_change") == "released"]
    if not released:
        return None
    return max(released, key=lambda h: float(h["version"]))["uuid"]


def resolve_superseded(file_ids: Iterable[str], session: requests.Session | None = None) -> dict[str, str]:
    """Map each (old) file id to its current released id. Ids with no released version are omitted."""
    http = session or requests.Session()
    mapping = {}
    for fid in tqdm(list(file_ids), desc="Resolving file versions"):
        current = latest_file_id(fid, session=http)
        if current:
            mapping[fid] = current
    return mapping


def _safe_members(tar: tarfile.TarFile, dest: Path) -> list[tarfile.TarInfo]:
    dest = dest.resolve()
    members = []
    for m in tar.getmembers():
        target = (dest / m.name).resolve()
        if not target.is_relative_to(dest) or m.issym() or m.islnk() or m.isdev():
            raise ValueError(f"Refusing unsafe archive member: {m.name}")
        members.append(m)
    return members


def _extract(archive: Path, dest: Path) -> None:
    with tarfile.open(archive, "r:gz") as tar:
        members = _safe_members(tar, dest)
        if hasattr(tarfile, "data_filter"):
            tar.extractall(dest, members=members, filter="data")
        else:  # Python < 3.11.4
            tar.extractall(dest, members=members)
    (dest / "MANIFEST.txt").unlink(missing_ok=True)


def downloaded_path(raw_dir: Path, file_id: str, file_name: str) -> Path:
    """GDC layout used for both bulk tarballs and single downloads: <raw_dir>/<file_id>/<file_name>."""
    return Path(raw_dir) / file_id / file_name


def download_files(
    files: pd.DataFrame,
    raw_dir: str | Path,
    batch_size: int = 25,
    session: requests.Session | None = None,
) -> None:
    """Stream files to <raw_dir>/<file_id>/<file_name>, skipping ones already present.

    Downloads in batches so a failure partway through only loses one batch, and
    streams to disk rather than holding a multi-GB tarball in memory.
    """
    http = session or requests.Session()
    raw_dir = Path(raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)

    todo = files[[not downloaded_path(raw_dir, r.file_id, r.file_name).exists()
                  for r in files.itertuples()]]
    if todo.empty:
        log.info("All %d files already downloaded", len(files))
        return

    for i in tqdm(range(0, len(todo), batch_size), desc="GDC batches"):
        batch = todo.iloc[i : i + batch_size]
        ids = batch["file_id"].tolist()
        with http.post(f"{API_URL}/data", json={"ids": ids}, stream=True, timeout=600) as resp:
            resp.raise_for_status()
            disposition = resp.headers.get("Content-Disposition", "")
            match = re.search(r'filename="?([^";]+)"?', disposition)
            served_name = match.group(1) if match else ""
            with tempfile.NamedTemporaryFile(dir=raw_dir, delete=False) as tmp:
                for chunk in resp.iter_content(chunk_size=1 << 20):
                    tmp.write(chunk)
                tmp_path = Path(tmp.name)
        try:
            if served_name.endswith((".tar.gz", ".tgz")):
                _extract(tmp_path, raw_dir)
            else:  # a single-id request returns the bare file, not an archive
                row = batch.iloc[0]
                target = downloaded_path(raw_dir, row.file_id, row.file_name)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(tmp_path, target)
        finally:
            tmp_path.unlink(missing_ok=True)
