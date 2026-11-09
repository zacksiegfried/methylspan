"""Read GDC SeSAMe beta-value files into a single probes x samples matrix, cached as parquet."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from .gdc import downloaded_path

log = logging.getLogger(__name__)


def read_beta_file(path: str | Path, name: str | None = None) -> pd.Series:
    """One GDC `*.level3betas.txt` file (two tab-separated columns, no header) as a float32 Series."""
    s = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["probe", "beta"],
        index_col="probe",
        dtype={"beta": np.float32},
        na_values=["NA"],
    )["beta"]
    s.name = name
    return s


def _index_raw_dir(raw_dir: Path) -> dict[str, Path]:
    """file_name -> path for everything under raw_dir, to find files whose GDC file_id changed."""
    return {p.name: p for p in raw_dir.glob("*/*") if p.is_file()}


def unmatched_dirs(files: pd.DataFrame, raw_dir: str | Path) -> list[str]:
    """Folder names under raw_dir that aren't a file_id in `files` (e.g. superseded GDC versions)."""
    known = set(files["file_id"])
    return sorted(p.name for p in Path(raw_dir).iterdir() if p.is_dir() and p.name not in known)


def _single_file(folder: Path) -> Path | None:
    found = [p for p in folder.iterdir() if p.is_file() and not p.name.startswith(".")]
    return found[0] if len(found) == 1 else None


def locate_files(files: pd.DataFrame, raw_dir: str | Path, aliases: dict[str, str] | None = None) -> pd.Series:
    """Path on disk for each row of `files` (None if missing). Index is file_id.

    Looks in <raw_dir>/<file_id>/<file_name>, then for the same file_name
    anywhere under raw_dir, then in <raw_dir>/<old_id>/ where `aliases` maps
    current file_id -> old (superseded) file id from `gdc.resolve_superseded`.
    """
    raw_dir = Path(raw_dir)
    by_name: dict[str, Path] | None = None
    aliases = aliases or {}
    paths = {}
    for r in files.itertuples():
        p = downloaded_path(raw_dir, r.file_id, r.file_name)
        if not p.exists():
            by_name = by_name if by_name is not None else _index_raw_dir(raw_dir)
            p = by_name.get(r.file_name)
        if p is None and r.file_id in aliases and (raw_dir / aliases[r.file_id]).is_dir():
            p = _single_file(raw_dir / aliases[r.file_id])
        paths[r.file_id] = p
    return pd.Series(paths, dtype=object)


def build_beta_matrix(
    files: pd.DataFrame,
    raw_dir: str | Path,
    join: str = "inner",
    column: str = "case_id",
    aliases: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Probes x samples float32 matrix for the files in `files`.

    Columns are labeled by `column` (case_id by default, so they line up with
    clinical data; this requires one file per case). join="inner" keeps only
    probes present on every file, which is what you want when combining 450K
    and EPIC arrays.
    """
    if files[column].duplicated().any():
        raise ValueError(f"Duplicate {column} values; run select_cohort(one_per_case=True) first")
    paths = locate_files(files, raw_dir, aliases)
    missing = paths.isna()
    if missing.any():
        log.warning("%d of %d files not found under %s; skipping them", missing.sum(), len(paths), raw_dir)
    labels = files.set_index("file_id")[column]
    series = [
        read_beta_file(path, name=labels[fid])
        for fid, path in tqdm(paths[~missing].items(), total=(~missing).sum(), desc="Reading betas")
    ]
    if not series:
        raise FileNotFoundError(f"None of the requested files are under {raw_dir}")
    matrix = pd.concat(series, axis=1, join=join, copy=False)
    matrix.index.name = "probe"
    matrix.columns.name = column
    return matrix


def save_matrix(matrix: pd.DataFrame, path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    matrix.to_parquet(path)


def load_matrix(path: str | Path) -> pd.DataFrame:
    return pd.read_parquet(path)


def filter_probes(matrix: pd.DataFrame, max_missing: float = 0.1) -> pd.DataFrame:
    """Drop probes missing in more than `max_missing` of samples. Uses no outcome data."""
    keep = matrix.isna().mean(axis=1) <= max_missing
    log.info("Keeping %d of %d probes (<= %.0f%% missing)", keep.sum(), len(keep), max_missing * 100)
    return matrix.loc[keep]
