"""Turn GDC case records into a one-row-per-patient survival table."""

from __future__ import annotations

import re

import numpy as np
import pandas as pd


def _primary_diagnosis(case: dict) -> dict:
    """The diagnosis flagged as the primary disease; GDC does not guarantee it is listed first."""
    diagnoses = case.get("diagnoses") or []
    for d in diagnoses:
        if d.get("diagnosis_is_primary_disease"):
            return d
    return diagnoses[0] if diagnoses else {}


def _max_or_nan(values) -> float:
    vals = [v for v in values if v is not None]
    return float(max(vals)) if vals else np.nan


def simplify_stage(stage: str | None) -> str | None:
    """'Stage IIB' -> 'II'. Returns None for missing or unparseable values."""
    if not stage:
        return None
    m = re.match(r"stage\s+(IV|I{1,3})", stage, flags=re.IGNORECASE)
    return m.group(1).upper() if m else None


def survival_table(cases: list[dict]) -> pd.DataFrame:
    """Index case_id; columns time (days), event (bool), and clinical covariates.

    time is days_to_death for patients who died, otherwise the latest follow-up
    found anywhere in the record (the diagnosis's days_to_last_follow_up or any
    follow_ups entry). Censored patients are kept; dropping them, or treating
    their follow-up time as a death time, biases survival estimates.
    """
    rows = []
    for c in cases:
        demo = c.get("demographic") or {}
        dx = _primary_diagnosis(c)
        vital = demo.get("vital_status")
        follow_up = _max_or_nan(
            [dx.get("days_to_last_follow_up")]
            + [f.get("days_to_follow_up") for f in c.get("follow_ups") or []]
        )
        if vital == "Dead":
            event, time = True, demo.get("days_to_death")
            time = np.nan if time is None else float(time)
        elif vital == "Alive":
            event, time = False, follow_up
        else:
            event, time = None, np.nan
        age = dx.get("age_at_diagnosis")
        rows.append(
            {
                "case_id": c["case_id"],
                "submitter_id": c.get("submitter_id"),
                "project_id": (c.get("project") or {}).get("project_id"),
                "time": time,
                "event": event,
                "age": np.nan if age is None else age / 365.25,
                "gender": demo.get("gender"),
                "stage": simplify_stage(dx.get("ajcc_pathologic_stage")),
                "grade": dx.get("tumor_grade"),
                "primary_diagnosis": dx.get("primary_diagnosis"),
            }
        )
    return pd.DataFrame(rows).set_index("case_id")


def usable_for_survival(table: pd.DataFrame) -> pd.DataFrame:
    """Rows with a known event status and a positive survival time."""
    ok = table["event"].notna() & table["time"].gt(0)
    return table.loc[ok].astype({"event": bool})
