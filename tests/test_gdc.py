import io
import tarfile

import pandas as pd
import pytest

from methylspan import gdc


class FakeResponse:
    def __init__(self, payload=None, content=b"", headers=None):
        self._payload = payload
        self._content = content
        self.headers = headers or {}

    def raise_for_status(self):
        pass

    def json(self):
        return self._payload

    def iter_content(self, chunk_size):
        yield self._content

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class PagedSession:
    """Serves `hits` in pages according to the request's size/from."""

    def __init__(self, hits):
        self.hits = hits
        self.calls = []

    def post(self, url, json=None, **kw):
        self.calls.append(json)
        start, size = json["from"], json["size"]
        page = self.hits[start : start + size]
        return FakeResponse({"data": {"hits": page, "pagination": {"total": len(self.hits)}}})


def _file_hit(i, case="c1", sample_type="Primary Tumor", platform="Illumina Human Methylation 450"):
    return {
        "file_id": f"f{i}",
        "file_name": f"f{i}.txt",
        "file_size": 10,
        "platform": platform,
        "cases": [{"case_id": case, "project": {"project_id": "TCGA-X"},
                   "samples": [{"sample_type": sample_type}]}],
    }


def test_paginate_follows_all_pages():
    session = PagedSession([{"i": i} for i in range(7)])
    hits = gdc.paginate("files", {}, ["i"], page_size=3, session=session)
    assert [h["i"] for h in hits] == list(range(7))
    assert [c["from"] for c in session.calls] == [0, 3, 6]


def test_query_methylation_files_flattens_and_handles_missing_samples():
    hits = [_file_hit(1), {**_file_hit(2), "cases": [{"case_id": "c2", "project": {"project_id": "P"}}]}]
    df = gdc.query_methylation_files(primary_site="pancreas", session=PagedSession(hits))
    assert list(df["platform"]) == ["450K", "450K"]
    assert df.loc[1, "sample_type"] is None


@pytest.mark.parametrize("raw,expected", [
    ("Illumina Human Methylation 450", "450K"),
    ("Illumina Methylation Epic", "EPIC"),
    ("Illumina Methylation Epic v2", "EPICv2"),
    ("Illumina Human Methylation 27", "27K"),
    (None, None),
])
def test_normalize_platform(raw, expected):
    assert gdc.normalize_platform(raw) == expected


def _tarball(members: dict[str, bytes]) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return buf.getvalue()


class DownloadSession:
    def __init__(self, content, filename):
        self.content, self.filename = content, filename
        self.requested = []

    def post(self, url, json=None, **kw):
        self.requested.append(json["ids"])
        return FakeResponse(content=self.content,
                            headers={"Content-Disposition": f"attachment; filename={self.filename}"})


def test_download_extracts_tarball_and_skips_existing(tmp_path):
    files = pd.DataFrame({"file_id": ["a", "b"], "file_name": ["a.txt", "b.txt"]})
    content = _tarball({"a/a.txt": b"cg1\t0.5\n", "b/b.txt": b"cg1\t0.6\n", "MANIFEST.txt": b"x"})
    session = DownloadSession(content, "gdc_download.tar.gz")
    gdc.download_files(files, tmp_path, session=session)
    assert (tmp_path / "a" / "a.txt").read_bytes() == b"cg1\t0.5\n"
    assert not (tmp_path / "MANIFEST.txt").exists()

    gdc.download_files(files, tmp_path, session=session)
    assert len(session.requested) == 1  # second call found everything on disk


def test_download_single_file_is_not_a_tarball(tmp_path):
    files = pd.DataFrame({"file_id": ["a"], "file_name": ["a.txt"]})
    gdc.download_files(files, tmp_path, session=DownloadSession(b"cg1\t0.5\n", "a.txt"))
    assert (tmp_path / "a" / "a.txt").read_bytes() == b"cg1\t0.5\n"


def test_download_rejects_path_traversal(tmp_path):
    files = pd.DataFrame({"file_id": ["a", "b"], "file_name": ["a.txt", "b.txt"]})
    content = _tarball({"../evil.txt": b"x", "b/b.txt": b"y"})
    with pytest.raises(ValueError, match="unsafe"):
        gdc.download_files(files, tmp_path / "raw", session=DownloadSession(content, "x.tar.gz"))
    assert not (tmp_path / "evil.txt").exists()


class HistorySession:
    def __init__(self, histories):
        self.histories = histories

    def get(self, url, **kw):
        return FakeResponse(self.histories[url.rsplit("/", 1)[-1]])


def test_resolve_superseded_picks_latest_released_version():
    session = HistorySession({
        "old": [{"uuid": "v2", "version": "2", "file_change": "released"},
                {"uuid": "v3", "version": "3", "file_change": "released"},
                {"uuid": "old", "version": "1", "file_change": "superseded"}],
        "gone": [{"uuid": "gone", "version": "1", "file_change": "superseded"}],
    })
    assert gdc.resolve_superseded(["old", "gone"], session=session) == {"old": "v3"}
