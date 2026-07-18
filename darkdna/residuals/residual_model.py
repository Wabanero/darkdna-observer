"""Residual anomaly scoring after classical covariate control."""

from __future__ import annotations

import json
import math
from pathlib import Path
from statistics import NormalDist

import numpy as np
import pandas as pd

from darkdna.utils.optional_deps import optional_import
from darkdna.utils.progress import ProgressReporter
from darkdna.views.primitive_scores import PRIMITIVE_SCORE_COLUMNS


def numeric_covariates(covariates: pd.DataFrame, target: str) -> list[str]:
    target_tokens = set(target.replace("_score", "").split("_"))
    cols = []
    for col in covariates.columns:
        if col in {"region_id", "chrom", "chrom_cov", "start", "end", "block_id", "chromosome_cv_fold"}:
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


def fit_predict_holdout(X_train: pd.DataFrame, y_train: pd.Series, X_test: pd.DataFrame, method: str = "linear") -> tuple[np.ndarray, str]:
    Xn = X_train.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    Xh = X_test.reindex(columns=Xn.columns).apply(pd.to_numeric, errors="coerce").fillna(0.0)
    yn = pd.to_numeric(y_train, errors="coerce").fillna(float(pd.to_numeric(y_train, errors="coerce").mean()))
    if len(yn) < 3 or Xn.shape[1] == 0:
        return np.repeat(float(yn.mean() if len(yn) else 0.0), len(Xh)), "blocked_mean"
    if method == "robust_linear":
        sm, warning = optional_import("statsmodels.api")
        if sm is not None:
            try:  # pragma: no cover - statsmodels may be absent.
                train = sm.add_constant(Xn, has_constant="add")
                test = sm.add_constant(Xh, has_constant="add")
                model = sm.RLM(yn, train).fit()
                return np.asarray(model.predict(test)), "robust_linear_block_cv"
            except Exception:
                pass
    sklearn_ensemble, _ = optional_import("sklearn.ensemble")
    if method == "random_forest" and sklearn_ensemble is not None:
        model = sklearn_ensemble.RandomForestRegressor(n_estimators=100, random_state=13, min_samples_leaf=2)
        model.fit(Xn, yn)
        return model.predict(Xh), "random_forest_block_cv"
    if method == "gradient_boosting" and sklearn_ensemble is not None:
        model = sklearn_ensemble.GradientBoostingRegressor(random_state=13)
        model.fit(Xn, yn)
        return model.predict(Xh), "gradient_boosting_block_cv"
    if method == "lightgbm":
        lightgbm, _ = optional_import("lightgbm")
        if lightgbm is not None:
            try:  # pragma: no cover - optional dependency.
                model = lightgbm.LGBMRegressor(random_state=13, verbosity=-1)
                model.fit(Xn, yn)
                return model.predict(Xh), "lightgbm_block_cv"
            except Exception:
                pass
    if method == "xgboost":
        xgboost, _ = optional_import("xgboost")
        if xgboost is not None:
            try:  # pragma: no cover - optional dependency.
                model = xgboost.XGBRegressor(random_state=13)
                model.fit(Xn, yn)
                return model.predict(Xh), "xgboost_block_cv"
            except Exception:
                pass
    return pure_numpy_linear_predict(Xn, yn, X_predict=Xh), "linear_gradient_descent_block_cv"


def select_cv_group(data: pd.DataFrame) -> str | None:
    for col in ["chrom", "chrom_cov"]:
        if col in data.columns and data[col].fillna("").astype(str).nunique() >= 2:
            return col
    if "block_id" in data.columns and data["block_id"].fillna("").astype(str).nunique() >= 2:
        return "block_id"
    return None


def fit_predict_blocked(X: pd.DataFrame, y: pd.Series, groups: pd.Series, method: str = "linear") -> tuple[np.ndarray, str]:
    group_values = groups.fillna("missing").astype(str)
    if group_values.nunique() < 2:
        return fit_predict(X, y, method=method)
    predictions = np.zeros(len(group_values), dtype=float)
    used_methods: list[str] = []
    y_numeric = pd.to_numeric(y, errors="coerce").fillna(float(pd.to_numeric(y, errors="coerce").mean()))
    for group in group_values.drop_duplicates():
        test_mask = group_values == group
        train_mask = ~test_mask
        if int(train_mask.sum()) < 3:
            pred = np.repeat(float(y_numeric.loc[train_mask].mean() if int(train_mask.sum()) else y_numeric.mean()), int(test_mask.sum()))
            used = "blocked_mean"
        else:
            pred, used = fit_predict_holdout(X.loc[train_mask], y_numeric.loc[train_mask], X.loc[test_mask], method=method)
        predictions[np.flatnonzero(test_mask.to_numpy())] = np.nan_to_num(pred, nan=float(y_numeric.mean()), posinf=0.0, neginf=0.0)
        used_methods.append(used)
    unique_methods = sorted(set(used_methods))
    used_method = unique_methods[0] if len(unique_methods) == 1 else "mixed_block_cv"
    return predictions, used_method


def pure_numpy_linear_predict(
    X: pd.DataFrame,
    y: pd.Series,
    X_predict: pd.DataFrame | None = None,
    steps: int = 800,
    learning_rate: float = 0.05,
) -> np.ndarray:
    """Small linear-regression fallback that avoids LAPACK-backed lstsq.

    Some Windows scientific stacks can crash in native BLAS/LAPACK calls. This
    optimizer is modest but stable for MVP residualization and smoke tests.
    """

    x = X.to_numpy(dtype=float)
    xp = X_predict.reindex(columns=X.columns).to_numpy(dtype=float) if X_predict is not None else x
    target = y.to_numpy(dtype=float)
    target_mean = float(np.nanmean(target))
    target_std = float(np.nanstd(target))
    if not np.isfinite(target_std) or target_std == 0.0:
        return np.repeat(target_mean if np.isfinite(target_mean) else 0.0, len(xp))
    target_z = (target - target_mean) / target_std
    means = x.mean(axis=0)
    stds = x.std(axis=0)
    stds[stds == 0] = 1.0
    z = (x - means) / stds
    zp = (xp - means) / stds
    z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)
    zp = np.nan_to_num(zp, nan=0.0, posinf=0.0, neginf=0.0)
    design = np.c_[np.ones(len(z)), z]
    predict_design = np.c_[np.ones(len(zp)), zp]
    weights = np.zeros(design.shape[1], dtype=float)
    for _ in range(steps):
        pred = design @ weights
        grad = (design.T @ (pred - target_z)) / max(1, len(target_z))
        grad = np.clip(np.nan_to_num(grad, nan=0.0, posinf=1e3, neginf=-1e3), -1e3, 1e3)
        weights -= learning_rate * grad
    pred_z = np.nan_to_num(predict_design @ weights, nan=0.0, posinf=0.0, neginf=0.0)
    return pred_z * target_std + target_mean


def _robust_residual_zscores(residuals: np.ndarray) -> np.ndarray:
    values = np.asarray(residuals, dtype=float)
    median = float(np.nanmedian(values))
    mad = float(np.nanmedian(np.abs(values - median)) / 0.6745)
    if not np.isfinite(mad) or mad <= 1e-12:
        mad = float(np.nanstd(values, ddof=1)) if len(values) > 1 else math.nan
    if not np.isfinite(mad) or mad <= 1e-12:
        return np.full(len(values), np.nan, dtype=float)
    return (values - median) / mad


def _conditional_residual_calibration(
    X: pd.DataFrame,
    residuals: np.ndarray,
    groups: pd.Series | None,
    method: str,
) -> tuple[np.ndarray, np.ndarray, str, str]:
    if groups is None or groups.fillna("missing").astype(str).nunique() < 2:
        unavailable = np.full(len(residuals), np.nan, dtype=float)
        return unavailable, unavailable, "unavailable", "Conditional variance calibration requires at least two chromosome/block groups."
    absolute = pd.Series(np.abs(residuals), index=X.index, dtype=float)
    predicted_absolute, used = fit_predict_blocked(X, absolute, groups, method=method)
    # E|N(0,sigma)| = sigma * sqrt(2/pi).  This is a descriptive variance
    # adapter; its method and assumptions remain explicit in the output.
    sigma = np.asarray(predicted_absolute, dtype=float) * math.sqrt(math.pi / 2.0)
    floor = float(np.nanmedian(np.abs(residuals)) * 0.1)
    floor = max(floor, 1e-9)
    sigma = np.where(np.isfinite(sigma) & (sigma > floor), sigma, floor)
    conditional_z = np.asarray(residuals, dtype=float) / sigma
    return conditional_z, sigma, f"available:{used}", "Cross-fitted absolute-residual model converted to a conditional scale estimate."


def _heldout_quantile_and_conformal(
    observed: np.ndarray,
    predicted: np.ndarray,
    residuals: np.ndarray,
    groups: pd.Series | None,
    *,
    alpha: float = 0.05,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, str, str]:
    count = len(residuals)
    quantile_residual = np.full(count, np.nan, dtype=float)
    lower = np.full(count, np.nan, dtype=float)
    upper = np.full(count, np.nan, dtype=float)
    if groups is None or groups.fillna("missing").astype(str).nunique() < 2:
        return quantile_residual, lower, upper, "unavailable", "Held-out quantile and conformal residuals require at least two chromosome/block groups."
    group_values = groups.fillna("missing").astype(str)
    for group in group_values.drop_duplicates():
        test_mask = (group_values == group).to_numpy()
        train_mask = ~test_mask
        training = np.asarray(residuals[train_mask], dtype=float)
        training = training[np.isfinite(training)]
        if len(training) < 5:
            continue
        training_abs = np.abs(training)
        conformal_quantile = float(np.quantile(training_abs, min(1.0, math.ceil((len(training_abs) + 1) * (1 - alpha)) / len(training_abs))))
        lower[test_mask] = predicted[test_mask] - conformal_quantile
        upper[test_mask] = predicted[test_mask] + conformal_quantile
        for index in np.flatnonzero(test_mask):
            rank_probability = (float(np.sum(training <= residuals[index])) + 0.5) / (len(training) + 1.0)
            rank_probability = float(np.clip(rank_probability, 1e-6, 1.0 - 1e-6))
            quantile_residual[index] = float(NormalDist().inv_cdf(rank_probability))
    status = "available" if np.isfinite(quantile_residual).all() else "partial"
    if not np.isfinite(quantile_residual).any():
        status = "unavailable"
    return quantile_residual, lower, upper, status, "Ranks and conformal intervals were calibrated only against other held-out groups' cross-fitted residuals."


def _heteroscedasticity_diagnostic(predicted: np.ndarray, residuals: np.ndarray) -> tuple[float, float, str]:
    predicted_values = np.asarray(predicted, dtype=float)
    absolute_residuals = np.abs(np.asarray(residuals, dtype=float))
    finite = np.isfinite(predicted_values) & np.isfinite(absolute_residuals)
    if int(finite.sum()) < 8 or float(np.nanstd(predicted_values[finite])) <= 1e-12:
        return math.nan, math.nan, "unavailable_requires_prediction_variation_and_at_least_8_rows"
    left = pd.Series(predicted_values[finite]).rank(method="average").to_numpy(dtype=float)
    right = pd.Series(absolute_residuals[finite]).rank(method="average").to_numpy(dtype=float)
    left_centered = left - float(left.mean())
    right_centered = right - float(right.mean())
    denominator = math.sqrt(float(np.sum(left_centered**2)) * float(np.sum(right_centered**2)))
    if denominator <= 1e-12:
        return math.nan, math.nan, "unavailable_degenerate_rank_distribution"
    correlation = float(np.sum(left_centered * right_centered) / denominator)
    return correlation, math.nan, "available_descriptive_spearman_abs_residual_vs_prediction_p_value_unavailable"


def _block_bootstrap_abs_residual_ci(
    residuals: np.ndarray,
    groups: pd.Series | None,
    *,
    n_bootstrap: int = 500,
    seed: int = 13,
) -> tuple[float, float, str]:
    if groups is None or groups.fillna("missing").astype(str).nunique() < 2:
        return math.nan, math.nan, "unavailable_requires_at_least_two_blocks"
    frame = pd.DataFrame({"group": groups.fillna("missing").astype(str), "absolute_residual": np.abs(residuals)})
    block_means = frame.groupby("group")["absolute_residual"].mean().to_numpy(dtype=float)
    if len(block_means) < 2:
        return math.nan, math.nan, "unavailable_requires_at_least_two_blocks"
    rng = np.random.default_rng(seed)
    estimates = [float(rng.choice(block_means, size=len(block_means), replace=True).mean()) for _ in range(max(1, n_bootstrap))]
    low, high = np.percentile(estimates, [2.5, 97.5])
    return float(low), float(high), "available_block_resampled_mean_absolute_residual"


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
    cv_group_col = select_cv_group(data)
    cv_strategy = f"leave_one_{cv_group_col}_out" if cv_group_col else "in_sample_classical_fit_no_independent_calibration"
    reporter = ProgressReporter("residualize", total=len(primitives)) if progress else None
    if reporter:
        reporter.start(f"residualizing primitives method={method}")
    for idx, primitive in enumerate(primitives, start=1):
        cov_cols = numeric_covariates(covariates, primitive)
        observed_series = pd.to_numeric(data[primitive], errors="coerce")
        if not observed_series.notna().any():
            for row_idx, region_id in enumerate(data["region_id"].astype(str)):
                rows.append(
                    {
                        "region_id": region_id,
                        "primitive": primitive,
                        "observed_score": math.nan,
                        "predicted_classical_score": math.nan,
                        "residual_score": math.nan,
                        "residual_zscore": math.nan,
                        "residual_zscore_method": "unavailable",
                        "matched_null_zscore": math.nan,
                        "empirical_p_value": math.nan,
                        "empirical_p_value_status": "unavailable",
                        "empirical_p_value_reason": "The primitive screen itself was unavailable; no residual or p-value was computed.",
                        "null_panel_status": "unavailable_input_score",
                        "classical_model_global_r2": math.nan,
                        "classical_explanation_fraction": math.nan,
                        "classical_explanation_fraction_deprecation_warning": (
                            "Deprecated alias for classical_model_global_r2; no value is available for this primitive."
                        ),
                        "covariates_used": ",".join(cov_cols),
                        "model_method": "unavailable_input_score",
                        "cv_strategy": cv_strategy,
                        "cv_group_col": cv_group_col or "",
                        "cv_group_id": str(data[cv_group_col].iloc[row_idx]) if cv_group_col and cv_group_col in data.columns else "",
                    }
                )
            summaries.append(
                {
                    "primitive": primitive,
                    "method": "unavailable_input_score",
                    "n_regions": len(data),
                    "n_covariates": len(cov_cols),
                    "covariates_used": ",".join(cov_cols),
                    "classical_model_global_r2": math.nan,
                    "classical_explanation_fraction": math.nan,
                    "residual_calibration_status": "unavailable",
                    "residual_calibration_reason": "All primitive scores were NA; unavailable values were not converted to zero.",
                    "cv_strategy": cv_strategy,
                    "cv_group_col": cv_group_col or "",
                }
            )
            if reporter:
                reporter.update(idx, message=f"{primitive} unavailable")
            continue
        if cv_group_col and cv_group_col in data.columns:
            pred, used_method = fit_predict_blocked(data[cov_cols], data[primitive], data[cv_group_col], method=method)
        else:
            pred, used_method = fit_predict(data[cov_cols], data[primitive], method=method)
        observed_mean = float(observed_series.mean()) if observed_series.notna().any() else 0.0
        obs = observed_series.fillna(observed_mean).to_numpy(dtype=float)
        pred = np.nan_to_num(pred, nan=float(np.nanmean(obs) if len(obs) else 0.0), posinf=0.0, neginf=0.0)
        residual = obs - pred
        robust_residual_z = _robust_residual_zscores(residual)
        group_series = data[cv_group_col] if cv_group_col and cv_group_col in data.columns else None
        conditional_z, conditional_sigma, conditional_status, conditional_reason = _conditional_residual_calibration(
            data[cov_cols], residual, group_series, method
        )
        quantile_residual, conformal_lower, conformal_upper, conformal_status, conformal_reason = _heldout_quantile_and_conformal(
            obs, pred, residual, group_series
        )
        residual_z = np.where(np.isfinite(conditional_z), conditional_z, robust_residual_z)
        residual_z_method = np.where(
            np.isfinite(conditional_z),
            "cross_fitted_conditional_variance_zscore",
            "descriptive_robust_unconditional_zscore_not_significance",
        )
        residual_std = np.std(residual, ddof=1) if len(residual) > 1 else math.nan
        ss_res = float(np.sum((obs - pred) ** 2))
        ss_tot = float(np.sum((obs - obs.mean()) ** 2))
        r2 = math.nan if ss_tot == 0 else float(1.0 - ss_res / ss_tot)
        hetero_correlation, hetero_p_value, hetero_status = _heteroscedasticity_diagnostic(pred, residual)
        bootstrap_low, bootstrap_high, bootstrap_status = _block_bootstrap_abs_residual_ci(residual, group_series)
        null_lookup = {}
        if nulls is not None and not nulls.empty:
            subset = nulls[nulls["primitive"] == primitive]
            null_lookup = subset.set_index("region_id").to_dict(orient="index")
        for row_idx, region_id in enumerate(data["region_id"].astype(str)):
            null_row = null_lookup.get(region_id, {})
            null_p_value = pd.to_numeric(pd.Series([null_row.get("empirical_p_value", np.nan)]), errors="coerce").iloc[0]
            null_p_status = str(null_row.get("empirical_p_value_status", ""))
            if not np.isfinite(null_p_value):
                null_p_value = math.nan
                null_p_status = null_p_status or "unavailable"
                null_p_reason = str(
                    null_row.get(
                        "empirical_p_value_reason",
                        "No explicit valid null distribution was supplied; observed genomic ranks are not substituted as p-values.",
                    )
                )
            else:
                null_p_status = null_p_status or "available_explicit_matched_control_null"
                null_p_reason = str(null_row.get("empirical_p_value_reason", "One-sided empirical tail under the named matched-control null."))
            centered_observed = float(obs[row_idx] - obs.mean())
            local_ratio = abs(float(pred[row_idx] - obs.mean())) / (abs(centered_observed) + 1e-12)
            rows.append(
                {
                    "region_id": region_id,
                    "primitive": primitive,
                    "observed_score": float(obs[row_idx]),
                    "predicted_classical_score": float(pred[row_idx]),
                    "residual_score": float(residual[row_idx]),
                    "residual_zscore": float(residual_z[row_idx]),
                    "residual_zscore_method": str(residual_z_method[row_idx]),
                    "robust_unconditional_residual_zscore": float(robust_residual_z[row_idx]),
                    "conditional_residual_scale": float(conditional_sigma[row_idx]),
                    "conditional_residual_zscore": float(conditional_z[row_idx]),
                    "conditional_variance_status": conditional_status,
                    "conditional_variance_reason": conditional_reason,
                    "quantile_residual": float(quantile_residual[row_idx]),
                    "conformal_interval_lower": float(conformal_lower[row_idx]),
                    "conformal_interval_upper": float(conformal_upper[row_idx]),
                    "conformal_outlier": bool(
                        np.isfinite(conformal_lower[row_idx])
                        and (obs[row_idx] < conformal_lower[row_idx] or obs[row_idx] > conformal_upper[row_idx])
                    ),
                    "conformal_residual_status": conformal_status,
                    "conformal_residual_reason": conformal_reason,
                    "matched_null_zscore": float(null_row.get("null_zscore", np.nan)),
                    "empirical_p_value": float(null_p_value),
                    "empirical_p_value_status": null_p_status,
                    "empirical_p_value_reason": null_p_reason,
                    "null_model_id": str(null_row.get("null_model_id", "")),
                    "null_panel_status": str(null_row.get("null_panel_status", "")),
                    "available_null_models": str(null_row.get("available_null_models", "")),
                    "missing_null_models": str(null_row.get("missing_null_models", null_row.get("missing_or_partial_null_models", ""))),
                    "missing_or_partial_null_models": str(null_row.get("missing_or_partial_null_models", "")),
                    "null_model_count": int(null_row.get("null_model_count", 0) or 0),
                    "null_model_agreement": float(null_row.get("null_model_agreement", np.nan)),
                    "null_model_conflict": bool(null_row.get("null_model_conflict", False)),
                    "classical_model_global_r2": r2,
                    "classical_explanation_fraction": r2,
                    "classical_explanation_fraction_deprecation_warning": (
                        "Deprecated alias for classical_model_global_r2; this model-level cross-validated R2 is repeated per row and is not a per-region explanation fraction."
                    ),
                    "local_prediction_contribution_ratio": float(local_ratio),
                    "local_prediction_contribution_ratio_interpretation": "descriptive_not_causal_variance_decomposition",
                    "covariates_used": ",".join(cov_cols),
                    "model_method": used_method,
                    "cv_strategy": cv_strategy,
                    "cv_group_col": cv_group_col or "",
                    "cv_group_id": str(data[cv_group_col].iloc[row_idx]) if cv_group_col and cv_group_col in data.columns else "",
                }
            )
        summaries.append(
            {
                "primitive": primitive,
                "method": used_method,
                "n_regions": len(data),
                "n_covariates": len(cov_cols),
                "covariates_used": ",".join(cov_cols),
                "classical_model_global_r2": r2,
                "classical_explanation_fraction": r2,
                "classical_explanation_fraction_deprecation_warning": (
                    "Deprecated alias for classical_model_global_r2; it is not a per-region explanation fraction."
                ),
                "unconditional_residual_std": float(residual_std),
                "residual_std": float(residual_std),
                "conditional_variance_status": conditional_status,
                "conditional_variance_reason": conditional_reason,
                "conformal_residual_status": conformal_status,
                "heteroscedasticity_spearman_abs_residual_vs_prediction": hetero_correlation,
                "heteroscedasticity_diagnostic_p_value": hetero_p_value,
                "heteroscedasticity_diagnostic_status": hetero_status,
                "block_bootstrap_mean_absolute_residual_ci_low": bootstrap_low,
                "block_bootstrap_mean_absolute_residual_ci_high": bootstrap_high,
                "block_bootstrap_status": bootstrap_status,
                "cv_strategy": cv_strategy,
                "cv_group_col": cv_group_col or "",
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
