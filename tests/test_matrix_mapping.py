import numpy as np
import pandas as pd
import pytest

from methylspan import mapping, matrix


def _write_beta(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{p}\t{v}\n" for p, v in rows))


def test_build_beta_matrix_inner_join_and_labels(tmp_path):
    _write_beta(tmp_path / "f1" / "a.txt", [("cg1", 0.1), ("cg2", "NA"), ("cg3", 0.3)])
    _write_beta(tmp_path / "f2" / "b.txt", [("cg1", 0.5), ("cg2", 0.6)])
    files = pd.DataFrame({"file_id": ["f1", "f2"], "file_name": ["a.txt", "b.txt"], "case_id": ["c1", "c2"]})
    m = matrix.build_beta_matrix(files, tmp_path)
    assert list(m.columns) == ["c1", "c2"]
    assert list(m.index) == ["cg1", "cg2"]
    assert m.dtypes.eq(np.float32).all()
    assert np.isnan(m.loc["cg2", "c1"])

    out = tmp_path / "m.parquet"
    matrix.save_matrix(m, out)
    pd.testing.assert_frame_equal(matrix.load_matrix(out), m)


def test_build_beta_matrix_finds_files_by_name_when_file_id_changed(tmp_path):
    _write_beta(tmp_path / "old-id" / "a.txt", [("cg1", 0.1)])
    files = pd.DataFrame({"file_id": ["new-id"], "file_name": ["a.txt"], "case_id": ["c1"]})
    assert matrix.build_beta_matrix(files, tmp_path).loc["cg1", "c1"] == pytest.approx(0.1)


def test_build_beta_matrix_rejects_duplicate_patients(tmp_path):
    files = pd.DataFrame({"file_id": ["f1", "f2"], "file_name": ["a", "b"], "case_id": ["c1", "c1"]})
    with pytest.raises(ValueError, match="Duplicate"):
        matrix.build_beta_matrix(files, tmp_path)


def test_filter_probes():
    m = pd.DataFrame({"a": [0.1, np.nan], "b": [0.2, np.nan]}, index=["cg1", "cg2"])
    assert list(matrix.filter_probes(m, 0.5).index) == ["cg1"]


MANIFEST = """[Heading]
Descriptor File Name,x
[Assay]
IlmnID,Name,CHR,MAPINFO,UCSC_RefGene_Name,UCSC_RefGene_Group,Relation_to_UCSC_CpG_Island
cg1,cg1,1,100,GENEA;GENEB,TSS200;Body,Island
cg2,cg2,1,200,GENEA;GENEA,TSS1500;5'UTR,Island
cg3,cg3,1,300,,,
cg4,cg4,1,400,GENEB,TSS200,
[Controls]
ctl,stuff
"""


@pytest.fixture
def manifest(tmp_path):
    p = tmp_path / "manifest.csv"
    p.write_text(MANIFEST)
    return mapping.load_manifest(p)


def test_load_manifest_handles_raw_illumina_sections(manifest):
    assert list(manifest.index) == ["cg1", "cg2", "cg3", "cg4"]


def test_probe_gene_table_splits_and_pairs_regions(manifest):
    pg = mapping.probe_gene_table(manifest)
    pairs = set(zip(pg["probe"], pg["gene"], pg["region"]))
    assert ("cg1", "GENEA", "TSS200") in pairs
    assert ("cg1", "GENEB", "Body") not in pairs  # gene-body hit excluded from promoter set
    assert ("cg4", "GENEB", "TSS200") in pairs
    assert len(mapping.probe_gene_table(manifest, regions=None)) == 5


def test_aggregate_by_gene_counts_each_probe_once(manifest):
    pg = mapping.probe_gene_table(manifest)
    betas = pd.DataFrame({"s1": [0.2, 0.4, 0.9, 0.8]}, index=["cg1", "cg2", "cg3", "cg4"])
    genes = mapping.aggregate_by_gene(betas, pg, min_probes=2)
    assert list(genes.index) == ["GENEA"]  # GENEB has only one promoter probe
    assert genes.loc["GENEA", "s1"] == pytest.approx(0.3)  # cg2 hits GENEA twice but counts once
    assert list(mapping.gene_probes(betas, pg, "GENEA").index) == ["cg1", "cg2"]


def test_build_beta_matrix_uses_superseded_aliases(tmp_path):
    _write_beta(tmp_path / "old-id" / "old-name.txt", [("cg1", 0.7)])
    files = pd.DataFrame({"file_id": ["new-id"], "file_name": ["new-name.txt"], "case_id": ["c1"]})
    assert matrix.unmatched_dirs(files, tmp_path) == ["old-id"]
    m = matrix.build_beta_matrix(files, tmp_path, aliases={"new-id": "old-id"})
    assert m.loc["cg1", "c1"] == pytest.approx(0.7)
