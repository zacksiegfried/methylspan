import numpy as np
import pandas as pd

from methylspan import clinical, cohort


def _files():
    return pd.DataFrame(
        {
            "file_id": ["f1", "f2", "f3", "f4", "f5"],
            "case_id": ["c1", "c1", "c2", "c3", "c3"],
            "project_id": ["P"] * 5,
            "platform": ["450K", "450K", "EPIC", "450K", "450K"],
            "sample_type": ["Solid Tissue Normal", "Primary Tumor", "Primary Tumor",
                            "Metastatic", "Primary Tumor"],
        }
    )


def test_select_cohort_defaults_to_one_primary_tumor_per_patient():
    out = cohort.select_cohort(_files())
    assert sorted(out["file_id"]) == ["f2", "f3", "f5"]
    assert out["case_id"].is_unique


def test_select_cohort_platform_filter():
    out = cohort.select_cohort(_files(), platforms=["450K"])
    assert sorted(out["case_id"]) == ["c1", "c3"]


def test_select_cohort_prefers_primary_tumor_when_keeping_all_types():
    out = cohort.select_cohort(_files(), sample_types=None)
    assert dict(zip(out["case_id"], out["file_id"])) == {"c1": "f2", "c2": "f3", "c3": "f5"}


def test_survival_table_uses_primary_diagnosis_and_latest_follow_up():
    cases = [
        {   # alive: follow-up comes from the follow_ups list, not the first diagnosis
            "case_id": "a",
            "demographic": {"vital_status": "Alive"},
            "diagnoses": [
                {"diagnosis_is_primary_disease": False, "age_at_diagnosis": 1},
                {"diagnosis_is_primary_disease": True, "age_at_diagnosis": 365.25 * 60,
                 "days_to_last_follow_up": 100.0, "ajcc_pathologic_stage": "Stage IIB"},
            ],
            "follow_ups": [{"days_to_follow_up": 400}, {"days_to_follow_up": None}],
        },
        {
            "case_id": "d",
            "demographic": {"vital_status": "Dead", "days_to_death": 250},
            "diagnoses": [{"ajcc_pathologic_stage": "Stage IV"}],
        },
        {"case_id": "u", "demographic": {}, "diagnoses": []},
    ]
    t = clinical.survival_table(cases)
    assert t.loc["a", "time"] == 400 and t.loc["a", "event"] is False
    assert t.loc["a", "age"] == 60 and t.loc["a", "stage"] == "II"
    assert t.loc["d", "time"] == 250 and t.loc["d", "event"] is True
    assert np.isnan(t.loc["u", "time"])
    assert list(clinical.usable_for_survival(t).index) == ["a", "d"]


def test_simplify_stage():
    assert clinical.simplify_stage("Stage IA") == "I"
    assert clinical.simplify_stage("Stage III") == "III"
    assert clinical.simplify_stage("Stage IV") == "IV"
    assert clinical.simplify_stage("Not Reported") is None
    assert clinical.simplify_stage(None) is None
