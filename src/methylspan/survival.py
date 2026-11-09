"""Cox models on clinical and methylation features, evaluated with all preprocessing inside each CV fold."""

from __future__ import annotations

import time
import warnings
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin, clone
from sklearn.compose import ColumnTransformer
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.model_selection import GridSearchCV, KFold, RepeatedStratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sksurv.linear_model import CoxnetSurvivalAnalysis, CoxPHSurvivalAnalysis
from sksurv.metrics import (
    concordance_index_censored,
    concordance_index_ipcw,
    cumulative_dynamic_auc,
    integrated_brier_score,
)
from sksurv.util import Surv

from .clinical import usable_for_survival

CLINICAL_NUMERIC = ("age",)
CLINICAL_CATEGORICAL = ("stage",)


def make_y(table: pd.DataFrame) -> np.ndarray:
    """Structured survival array (event, time) for scikit-survival."""
    return Surv.from_arrays(event=table["event"].to_numpy(bool), time=table["time"].to_numpy(float))


def assemble(
    clinical: pd.DataFrame,
    features: pd.DataFrame | None = None,
    clinical_columns: Sequence[str] = CLINICAL_NUMERIC + CLINICAL_CATEGORICAL,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Join clinical covariates with a samples x features matrix (e.g. `betas.T`) on case_id.

    Keeps only patients present in both with usable survival data. Returns (X, y).
    """
    table = usable_for_survival(clinical)
    X = table[list(clinical_columns)]
    if features is not None:
        X = X.join(features, how="inner")
    if not X.index.is_unique:
        raise ValueError("Duplicate patients in X; each row must be a distinct case")
    return X, make_y(table.loc[X.index])


class TopVarianceSelector(BaseEstimator, TransformerMixin):
    """Keep the `k` most variable columns among those with at most `max_missing` NaNs.

    Fit on the training fold only, so test samples never influence which probes are used.
    """

    def __init__(self, k: int = 5000, max_missing: float = 0.1):
        self.k = k
        self.max_missing = max_missing

    def fit(self, X, y=None):
        X = np.asarray(X, dtype=np.float32)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN columns
            var = np.nanvar(X, axis=0)
        eligible = (np.isnan(X).mean(axis=0) <= self.max_missing) & np.isfinite(var)
        score = np.where(eligible, var, -np.inf)
        k = min(self.k, int(eligible.sum()))
        if k == 0:
            raise ValueError("No columns pass the missingness filter")
        self.support_ = np.sort(np.argpartition(score, -k)[-k:])
        self.n_features_in_ = X.shape[1]
        return self

    def transform(self, X):
        return np.asarray(X, dtype=np.float32)[:, self.support_]


def methylation_pipeline(n_top: int = 5000, n_components: int | None = None, max_missing: float = 0.1) -> Pipeline:
    steps = [
        ("select", TopVarianceSelector(n_top, max_missing)),
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
    ]
    if n_components:
        steps.append(("pca", PCA(n_components=n_components)))
    return Pipeline(steps)


def make_model(
    use_clinical: bool = True,
    use_methylation: bool = True,
    clinical_numeric: Sequence[str] = CLINICAL_NUMERIC,
    clinical_categorical: Sequence[str] = CLINICAL_CATEGORICAL,
    n_top: int = 5000,
    n_components: int | None = None,
    estimator: BaseEstimator | None = None,
) -> Pipeline:
    """Feature pipeline + Cox estimator. Every column not named as clinical is treated as methylation.

    The default estimator is a lightly ridge-penalized Cox model; pass
    `coxnet()` for an elastic net over many probes.
    """
    if not (use_clinical or use_methylation):
        raise ValueError("Model needs at least one feature group")
    clinical_cols = list(clinical_numeric) + list(clinical_categorical)
    if use_clinical:
        transformers = [
            ("num", Pipeline([("impute", SimpleImputer(strategy="median")),
                              ("scale", StandardScaler())]), list(clinical_numeric)),
            ("cat", Pipeline([("impute", SimpleImputer(strategy="constant", fill_value="missing")),
                              ("onehot", OneHotEncoder(handle_unknown="ignore", sparse_output=False))]),
             list(clinical_categorical)),
        ]
    else:
        transformers = [("clinical", "drop", clinical_cols)]
    remainder = methylation_pipeline(n_top, n_components) if use_methylation else "drop"
    features = ColumnTransformer(transformers, remainder=remainder)
    return Pipeline([("features", features), ("model", estimator or CoxPHSurvivalAnalysis(alpha=0.1))])


def coxnet(l1_ratio: float = 0.5, alpha: float = 0.1) -> CoxnetSurvivalAnalysis:
    """Elastic-net Cox with a single alpha, so it can be tuned with `tune`."""
    return CoxnetSurvivalAnalysis(l1_ratio=l1_ratio, alphas=[alpha], fit_baseline_model=True)


def tune(model: Pipeline, alphas: Sequence[float] = tuple(np.logspace(-1, 0.5, 6)), inner_splits: int = 3,
         random_state: int = 0) -> GridSearchCV:
    """Choose a coxnet alpha by inner CV (Harrell's C). Wrap this, not the bare model, in `cross_validate`.

    With thousands of standardized probes and ~150 patients, alphas much below
    0.1 leave the model nearly unpenalized: coxnet then iterates for minutes
    before failing to converge, so the default grid starts at 0.1.
    """
    return GridSearchCV(
        model,
        {"model__alphas": [[a] for a in alphas]},
        cv=KFold(inner_splits, shuffle=True, random_state=random_state),
        error_score=0.5,  # coxnet can fail to converge at tiny alphas; score those as random
    )


def _eval_times(y_train: np.ndarray, y_test: np.ndarray, n: int = 15) -> np.ndarray:
    """Time grid inside the follow-up range of both folds, from the test fold's event times.

    Capped at the last training-fold death, beyond which a Cox model's baseline
    survival (and so its predicted survival curves) is undefined.
    """
    event_times = y_test["time"][y_test["event"]]
    train_event_times = y_train["time"][y_train["event"]]
    if len(event_times) < 2 or len(train_event_times) == 0:
        return np.array([])
    grid = np.unique(np.quantile(event_times, np.linspace(0.1, 0.8, n)))
    lo = max(y_test["time"].min(), y_train["time"].min())
    hi = min(y_test["time"].max(), y_train["time"].max(), train_event_times.max())
    return grid[(grid > lo) & (grid < hi)]


def _survival_matrix(model, X: pd.DataFrame, times: np.ndarray) -> np.ndarray | None:
    est = getattr(model, "best_estimator_", model)
    if not hasattr(est, "predict_survival_function"):
        return None
    try:
        return np.vstack([fn(times) for fn in est.predict_survival_function(X)])
    except (ValueError, AttributeError):
        return None


def cross_validate(
    model,
    X: pd.DataFrame,
    y: np.ndarray,
    n_splits: int = 5,
    n_repeats: int = 1,
    random_state: int = 0,
) -> pd.DataFrame:
    """Per-fold Harrell's C, Uno's C, mean time-dependent AUC and integrated Brier score.

    The whole model (probe selection, imputation, scaling, PCA, any tuning) is
    refit on each training fold. Folds are stratified on the event indicator.
    """
    if not X.index.is_unique:
        raise ValueError("Duplicate patients in X; each row must be a distinct case")
    cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=random_state)
    rows = []
    for fold, (tr, te) in enumerate(cv.split(X, y["event"])):
        y_tr, y_te = y[tr], y[te]
        row = {
            "fold": fold,
            "n_test": len(te),
            "events_test": int(y_te["event"].sum()),
            "harrell_c": np.nan,
            "uno_c": np.nan,
            "mean_auc": np.nan,
            "ibs": np.nan,
            "error": None,
        }
        try:
            fitted = clone(model).fit(X.iloc[tr], y_tr)
        except (ArithmeticError, ValueError, np.linalg.LinAlgError) as exc:
            # One non-converging fold shouldn't discard the others; it shows up as NaN + error.
            row["error"] = f"{type(exc).__name__}: {exc}"
            rows.append(row)
            continue
        risk = fitted.predict(X.iloc[te])
        times = _eval_times(y_tr, y_te)
        row["harrell_c"] = concordance_index_censored(y_te["event"], y_te["time"], risk)[0]
        if len(times):
            row["uno_c"] = concordance_index_ipcw(y_tr, y_te, risk, tau=times[-1])[0]
            # IPCW weights come from the training fold's censoring distribution, which is
            # undefined past its last follow-up, so the few test patients followed longer
            # are left out of the time-dependent AUC and Brier score only.
            inside = y_te["time"] < y_tr["time"].max()
            row["mean_auc"] = cumulative_dynamic_auc(y_tr, y_te[inside], risk[inside], times)[1]
            surv = _survival_matrix(fitted, X.iloc[te[inside]], times)
            if surv is not None and len(times) >= 2:
                row["ibs"] = integrated_brier_score(y_tr, y_te[inside], surv, times)
        if hasattr(fitted, "best_params_"):
            row["best_params"] = str(fitted.best_params_)
        rows.append(row)
    return pd.DataFrame(rows)


def compare_models(models: Mapping[str, object], X: pd.DataFrame, y: np.ndarray, verbose: bool = True,
                   **cv_kwargs) -> pd.DataFrame:
    """Cross-validate several models on identical folds. Returns per-fold results with a `model` column."""
    results = []
    for name, m in models.items():
        start = time.perf_counter()
        res = cross_validate(m, X, y, **cv_kwargs).assign(model=name)
        if verbose:
            failed = res["error"].notna().sum()
            note = f", {failed} folds failed" if failed else ""
            print(f"{name}: C={res['harrell_c'].mean():.3f} ({time.perf_counter() - start:.0f}s{note})", flush=True)
        results.append(res)
    return pd.concat(results, ignore_index=True)


def summarize(results: pd.DataFrame) -> pd.DataFrame:
    """Mean and SD of each metric per model, plus folds that failed to fit.

    Higher C/AUC is better (0.5 = random); lower IBS is better.
    """
    metrics = ["harrell_c", "uno_c", "mean_auc", "ibs"]
    grouped = results.groupby("model", sort=False)
    summary = grouped[metrics].agg(["mean", "std"]).round(3)
    summary[("failed_folds", "")] = grouped["error"].apply(lambda e: int(e.notna().sum()))
    return summary
