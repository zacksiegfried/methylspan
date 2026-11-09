"""Turn a raw GDC file listing into an analysis cohort: one primary tumor per patient, one platform."""

from __future__ import annotations

import logging
from typing import Sequence

import pandas as pd

log = logging.getLogger(__name__)

# Lower rank wins when a patient has several files.
SAMPLE_TYPE_PRIORITY = {
    "Primary Tumor": 0,
    "Primary Blood Derived Cancer - Peripheral Blood": 1,
    "Primary Blood Derived Cancer - Bone Marrow": 1,
    "Recurrent Tumor": 2,
    "Metastatic": 3,
}


def describe(files: pd.DataFrame) -> pd.DataFrame:
    """Counts by project, platform and sample type. Look at this before modeling anything."""
    return (
        files.fillna({"sample_type": "(none)", "platform": "(none)"})
        .groupby(["project_id", "platform", "sample_type"])
        .agg(files=("file_id", "size"), patients=("case_id", "nunique"))
    )


def select_cohort(
    files: pd.DataFrame,
    sample_types: Sequence[str] | None = ("Primary Tumor",),
    platforms: Sequence[str] | None = None,
    projects: Sequence[str] | None = None,
    one_per_case: bool = True,
) -> pd.DataFrame:
    """Filter a `query_methylation_files` table down to a modeling cohort.

    Defaults to primary tumors only and one file per patient. Survival models
    assume each row is an independent patient; a patient with a tumor and a
    normal-tissue sample would otherwise appear twice and could land in both
    the training and test folds.

    Mixing platforms (450K and EPIC) is allowed but logged, because platform
    differences show up as large principal components.
    """
    out = files
    steps = [("start", len(out))]
    if sample_types is not None:
        out = out[out["sample_type"].isin(sample_types)]
        steps.append((f"sample_type in {list(sample_types)}", len(out)))
    if platforms is not None:
        out = out[out["platform"].isin(platforms)]
        steps.append((f"platform in {list(platforms)}", len(out)))
    if projects is not None:
        out = out[out["project_id"].isin(projects)]
        steps.append((f"project in {list(projects)}", len(out)))
    if one_per_case:
        rank = out["sample_type"].map(SAMPLE_TYPE_PRIORITY).fillna(99)
        out = (
            out.assign(_rank=rank)
            .sort_values(["case_id", "_rank", "file_id"])  # file_id tiebreak keeps it deterministic
            .drop_duplicates("case_id")
            .drop(columns="_rank")
        )
        steps.append(("one file per patient", len(out)))

    for label, n in steps:
        log.info("%-40s %5d files", label, n)
    mixed = out["platform"].nunique()
    if mixed > 1:
        log.warning("Cohort mixes %d platforms %s; expect batch effects", mixed,
                    sorted(out["platform"].dropna().unique()))
    return out.reset_index(drop=True)
