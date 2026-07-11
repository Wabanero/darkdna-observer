"""Residual anomaly scoring after classical covariate control."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from darkdna.utils.optional_deps import optional_import
from darkdna.utils.progress import ProgressReporter
from darkdna.utils.stats import empirical_p_value
from darkdna.views.primitive_scores import PRIMITIVE_SCORE_COLUMNS


def numeric_covariates(covariates: pd.DataFrame, target: str) -> list[str]:
    target_tokens = set(target.replace("_score", "").split("_"))
    cols = []
    for col in covariates.columns:
        if col in {"region_id", "chrom"}:
            continue
        if col == target or target.replace("_score", "") in col:
            continue
        if target_tokens and len(target_tokens.intersection(set(col.split("_")))) >= 2:
            continue
        if pd.api.types.is_bool_dtype(covariates[col]) or pd.api.types.is_numeric_dtype(covariates[col]):
            cols.append(col)
    return cols


def fit_predict(X: pd.DataFrame, y: pd.Series, method: str = "linear") -> tuple[np.ndarray, str]:
    Xn = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    yn = pd.to_numeric(y, errors="coerce").fillna(float(pd.to_numeric(y, errors="coerce").mean()))
    if len(yn) < 3 or Xn.shape[1] == 0:
        return np.repeat(float(yn.mean()), len(yn)), "mean"
    if method == "robust_linear":
        sm, warning = optional_import("statsmodels.api")
        if sm is not None:
            try:  # pragma: no cover - statsmodels may be absent.
                Xc = sm.add_constant(Xn, has_constant="add")
                model = sm.RLM(yn, Xc).fit()
                return np.asarray(model.predict(Xc)), "robust_linear"
            except Exception:
                pass
    sklearn_ensemble, _ = optional_import("sklearn.ensemble")
    if method == "random_forest" and sklearn_ensemble is not None:
        model = sklearn_ensemble.RandomForestRegressor(n_estimators=100, random_state=13, min_samples_leaf=2)
        model.fit(Xn, yn)
        return model.predict(Xn), "random_forest"
    if method == "gradient_boosting" and sklearn_ensemble is not None:
        model = sklearn_ensemble.GradientBoostingRegressor(random_state=13)
        model.fit(Xn, yn)
        return model.predict(Xn), "gradient_boosting"
    if method == "lightgbm":
        lightgbm, _ = optional_import("lightgbm")
        if lightgbm is not None:
            try:  # pragma: no cover - optional dependency.
                model = lightgbm.LGBMRegressor(random_state=13, verbosity=-1)
                model.fit(Xn, yn)
                return model.predict(Xn), "lightgbm"
            except Exception:
                pass
    if method == "xgboost":
        xgboost, _ = optional_import("xgboost")
        if xgboost is not None:
            try:  # pragma: no cover - optional dependency.
                model = xgboost.XGBRegressor(random_state=13)
                model.fit(Xn, yn)
                return model.predict(Xn), "xgboost"
            except Exception:
                pass
    return pure_numpy_linear_predict(Xn, yn), "linear_gradient_descent"


def pure_numpy_linear_predict(X: pd.DataFrame, y: pd.Series, steps: int = 800, learning_rate: float = 0.05) -> np.ndarray:
    """Small linear-regression fallback that avoids LAPACK-backed lstsq.

    Some Windows scientific stacks can crash in native BLAS/LAPACK calls. This
    optimizer is modest but stable for MVP residualization and smoke tests.
    """

    x = X.to_numpy(dtype=float)
    target = y.to_numpy(dtype=float)
    target_mean = float(np.nanmean(target))
    target_std = float(np.nanstd(target))
    if not np.isfinite(target_std) or target_std == 0.0:
        return np.repeat(target_mean if np.isfinite(target_mean) else 0.0, len(target))
    target_z = (target - target_mean) / target_std
    means = x.mean(axis=0)
    stds = x.std(axis=0)
    stds[stds == 0] = 1.0
    z = (x - means) / stds
    z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)
    design = np.c_[np.ones(len(z)), z]
    weights = np.zeros(design.shape[1], dtype=float)
    for _ in range(steps):
        pred = design @ weights
        grad = (design.T @ (pred - target_z)) / max(1, len(target_z))
        grad = np.clip(np.nan_to_num(grad, nan=0.0, posinf=1e3, neginf=-1e3), -1e3, 1e3)
        weights -= learning_rate * grad
    pred_z = np.nan_to_num(design @ weights, nan=0.0, posinf=0.0, neginf=0.0)
    return pred_z * target_std + target_mean


def residualize_scores(
    scores: pd.DataFrame,
    covariates: pd.DataFrame,
    nulls: pd.DataFrame | None = None,
    method: str = "linear",
    *,
    progress: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    data = scores.merge(covariates, on="region_id", how="left", suffixes=("", "_cov"))
    rows = []
    summaries = []
    primitives = [col for col in PRIMITIVE_SCORE_COLUMNS if col in scores.columns]
    reporter = ProgressReporter("residualize", total=len(primitives)) if progress else None
    if reporter:
        reporter.start(f"residualizing primitives method={method}")
    for idx, primitive in enumerate(primitives, start=1):
        cov_cols = numeric_covariates(covariates, primitive)
        pred, used_method = fit_predict(data[cov_cols], data[primitive], method=method)
        obs = pd.to_numeric(data[primitive], errors="coerce").fillna(0.0).to_numpy()
        pred = np.nan_to_num(pred, nan=float(np.nanmean(obs) if len(obs) else 0.0), posinf=0.0, neginf=0.0)
        residual = obs - pred
        residual_std = np.std(residual, ddof=1) if len(residual) > 1 else 0.0
        residual_z = residual / residual_std if residual_std and np.isfinite(residual_std) else np.zeros_like(residual)
        ss_res = float(np.sum((obs - pred) ** 2))
        ss_tot = float(np.sum((obs - obs.mean()) ** 2))
        r2 = 0.0 if ss_tot == 0 else max(0.0, 1.0 - ss_res / ss_tot)
        null_lookup = {}
        if nulls is not None and not nulls.empty:
            subset = nulls[nulls["primitive"] == primitive]
            null_lookup = subset.set_index("region_id").to_dict(orient="index")
        for row_idx, region_id in enumerate(data["region_id"].astype(str)):
            null_row = null_lookup.get(region_id, {})
            rows.append(
                {
                    "region_id": region_id,
                    "primitive": primitive,
                    "observed_score": float(obs[row_idx]),
                    "predicted_classical_score": float(pred[row_idx]),
                    "residual_score": float(residual[row_idx]),
                    "residual_zscore": float(residual_z[row_idx]),
                    "matched_null_zscore": float(null_row.get("null_zscore", np.nan)),
                    "empirical_p_value": float(null_row.get("empirical_p_value", empirical_p_value(obs[row_idx], obs, higher=True))),
                    "classical_explanation_fraction": r2,
                    "covariates_used": ",".join(cov_cols),
                    "model_method": used_method,
                }
            )
        summaries.append(
            {
                "primitive": primitive,
                "method": used_method,
                "n_regions": len(data),
                "n_covariates": len(cov_cols),
                "covariates_used": ",".join(cov_cols),
                "classical_explanation_fraction": r2,
                "residual_std": float(residual_std),
            }
        )
        if reporter:
            reporter.update(idx, message=f"{primitive} covariates={len(cov_cols)}")
    if reporter:
        reporter.finish()
    return pd.DataFrame(rows), pd.DataFrame(summaries)


def write_residual_outputs(residuals: pd.DataFrame, summary: pd.DataFrame, outdir: str | Path) -> dict[str, Path]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    residual_path = out / "residual_scores.parquet"
    summary_path = out / "residualization_summary.json"
    residuals.to_parquet(residual_path, index=False)
    summary_path.write_text(json.dumps(summary.to_dict(orient="records"), indent=2), encoding="utf-8")
    return {"residuals": residual_path, "summary": summary_path}
