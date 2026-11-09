"""Map Illumina CpG probes to genes and summarize methylation per gene region."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd

# UCSC_RefGene_Group values that sit in or just upstream of a promoter.
PROMOTER_REGIONS = ("TSS1500", "TSS200", "5'UTR", "1stExon")

_MANIFEST_COLUMNS = ["IlmnID", "CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group",
                     "Relation_to_UCSC_CpG_Island"]


def load_manifest(path: str | Path) -> pd.DataFrame:
    """Read an Illumina 450K/EPIC manifest CSV, raw (with [Heading]/[Controls] sections) or cleaned.

    Returns a frame indexed by probe ID with the gene annotation columns.
    """
    path = Path(path)
    header_row, n_rows = None, None
    with path.open() as fh:
        for i, line in enumerate(fh):
            if header_row is None and line.startswith("IlmnID,"):
                header_row = i
            elif header_row is not None and line.startswith("[Controls]"):
                n_rows = i - header_row - 1
                break
    if header_row is None:
        raise ValueError(f"No 'IlmnID' header row found in {path}")
    df = pd.read_csv(path, skiprows=header_row, nrows=n_rows, low_memory=False,
                     usecols=lambda c: c in _MANIFEST_COLUMNS, dtype=str)
    return df.set_index("IlmnID")


def probe_gene_table(manifest: pd.DataFrame, regions: Iterable[str] | None = PROMOTER_REGIONS) -> pd.DataFrame:
    """Long table of (probe, gene, region), one row per distinct probe-gene-region.

    The manifest packs several transcripts into one cell ('TSPY4;FAM197Y2' with
    groups 'Body;TSS1500'); this splits them and pairs each gene with its own
    region, so a probe in one gene's promoter and another gene's body is only
    counted for the first when regions=PROMOTER_REGIONS. Pass regions=None to keep all.
    """
    ann = manifest[["UCSC_RefGene_Name", "UCSC_RefGene_Group"]].dropna(subset=["UCSC_RefGene_Name"])
    names = ann["UCSC_RefGene_Name"].str.split(";")
    groups = ann["UCSC_RefGene_Group"].fillna("").str.split(";")
    long = pd.DataFrame({"gene": names, "region": groups}).explode(["gene", "region"])
    long.index.name = "probe"
    long = long.reset_index().drop_duplicates()
    long = long[long["gene"].ne("")]
    if regions is not None:
        long = long[long["region"].isin(list(regions))]
    return long.reset_index(drop=True)


def aggregate_by_gene(
    betas: pd.DataFrame,
    probe_genes: pd.DataFrame,
    min_probes: int = 2,
    how: str = "mean",
) -> pd.DataFrame:
    """Genes x samples matrix from a probes x samples matrix.

    Each probe counts once per gene even if it hits several transcripts' regions.
    Genes with fewer than `min_probes` measured probes are dropped.
    """
    pairs = probe_genes[["probe", "gene"]].drop_duplicates()
    pairs = pairs[pairs["probe"].isin(betas.index)]
    counts = pairs.groupby("gene").size()
    pairs = pairs[pairs["gene"].isin(counts[counts >= min_probes].index)]
    stacked = betas.reindex(pairs["probe"]).set_axis(pd.Index(pairs["gene"], name="gene"))
    return stacked.groupby(level="gene").agg(how)


def gene_probes(betas: pd.DataFrame, probe_genes: pd.DataFrame, gene: str) -> pd.DataFrame:
    """The individual probes x samples rows for one gene."""
    probes = probe_genes.loc[probe_genes["gene"].eq(gene), "probe"].unique()
    return betas.loc[betas.index.intersection(probes)]
