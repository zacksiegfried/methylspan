import numpy as np
import pandas as pd
import pytest

from methylspan import survival


@pytest.fixture
def synthetic():
    """120 patients; probe cg0 drives hazard, stage adds a little, 300 noise probes."""
    rng = np.random.default_rng(0)
    n, p = 120, 300
    betas = pd.DataFrame(rng.uniform(0, 1, (n, p)), columns=[f"cg{i}" for i in range(p)],
                         index=[f"case{i}" for i in range(n)]).astype(np.float32)
    betas.iloc[rng.choice(n, 10), 5] = np.nan
    stage = rng.choice(["I", "II", "III", None], n)
    risk = 3 * betas["cg0"].to_numpy() + (stage == "III")
    event_time = rng.exponential(500 / np.exp(risk))
    censor_time = rng.uniform(100, 1500, n)
    clinical = pd.DataFrame(
        {"time": np.minimum(event_time, censor_time) + 1, "event": event_time < censor_time,
         "age": rng.normal(65, 8, n), "stage": stage},
        index=betas.index,
    )
    return clinical, betas


def test_assemble_aligns_and_rejects_duplicates(synthetic):
    clinical, betas = synthetic
    X, y = survival.assemble(clinical, betas.iloc[:100])
    assert len(X) == len(y) == 100
    assert {"age", "stage", "cg0"} <= set(X.columns)
    with pytest.raises(ValueError, match="Duplicate"):
        survival.assemble(clinical, pd.concat([betas.iloc[:3], betas.iloc[:1]]))


def test_top_variance_selector_respects_missingness():
    X = np.array([[0.0, 0.0, 1.0], [1.0, np.nan, 1.0], [0.0, 9.0, 1.0], [1.0, np.nan, 1.0]])
    sel = survival.TopVarianceSelector(k=1, max_missing=0.25).fit(X)
    assert list(sel.support_) == [0]  # column 1 is most variable but half missing


def test_methylation_model_finds_signal_and_reports_all_metrics(synthetic):
    clinical, betas = synthetic
    X, y = survival.assemble(clinical, betas)
    models = {
        "clinical": survival.make_model(use_methylation=False),
        "methylation_only": survival.make_model(use_clinical=False, n_top=50),
        "clinical+pca": survival.make_model(n_top=300, n_components=5),
        "clinical+coxnet": survival.tune(
            survival.make_model(n_top=100, estimator=survival.coxnet()), alphas=[0.05, 0.2], inner_splits=2),
    }
    results = survival.compare_models(models, X, y, n_splits=3)
    summary = survival.summarize(results)
    assert list(summary.index) == list(models)
    assert results[["harrell_c", "uno_c", "mean_auc", "ibs"]].notna().all().all()
    # cg0 is a strong, high-variance-agnostic signal; coxnet over 100 probes should pick it up when present
    assert summary.loc["methylation_only", ("harrell_c", "mean")] > 0.5
    assert "best_params" in results.columns


def test_methylation_only_ignores_clinical_columns(synthetic):
    clinical, betas = synthetic
    X, y = survival.assemble(clinical, betas)
    model = survival.make_model(use_clinical=False, n_top=10).fit(X, y)
    name, _, cols = model.named_steps["features"].transformers_[-1]
    assert name == "remainder"
    cols = [X.columns[c] if isinstance(c, (int, np.integer)) else c for c in cols]
    assert len(cols) == betas.shape[1]
    assert "stage" not in cols and "age" not in cols


class _Exploding(survival.CoxPHSurvivalAnalysis):
    def fit(self, X, y):
        raise ArithmeticError("weights are too large")


def test_failed_fold_is_recorded_not_raised(synthetic):
    clinical, betas = synthetic
    X, y = survival.assemble(clinical)
    models = {"ok": survival.make_model(use_methylation=False),
              "bad": survival.make_model(use_methylation=False, estimator=_Exploding())}
    results = survival.compare_models(models, X, y, n_splits=3, verbose=False)
    summary = survival.summarize(results)
    assert summary.loc["bad", ("failed_folds", "")] == 3
    assert summary.loc["ok", ("failed_folds", "")] == 0
    assert results.loc[results.model == "bad", "error"].str.contains("weights").all()
