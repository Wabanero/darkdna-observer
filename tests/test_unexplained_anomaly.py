import numpy as np
import pandas as pd

from darkdna.views.unexplained_anomaly import cross_fitted_unexplained_outlierness


def test_unexplained_anomaly_is_cross_fitted_multivariate_outlierness_not_mean():
    rng = np.random.default_rng(13)
    matrix = pd.DataFrame(rng.normal(size=(30, 3)), columns=["known_a", "known_b", "known_c"])
    matrix.loc[29] = [9.0, -8.0, 7.0]
    groups = pd.Series(["b1"] * 10 + ["b2"] * 10 + ["b3"] * 10)

    scores, audit = cross_fitted_unexplained_outlierness(matrix, groups)

    assert audit["method"] == "cross_fitted_shrinkage_robust_mahalanobis"
    assert audit["status"] == "available"
    assert scores.notna().all()
    assert scores.iloc[29] > scores.iloc[:-1].quantile(0.95)
    assert not np.isclose(scores.iloc[29], matrix.iloc[29].mean())


def test_unexplained_anomaly_refuses_random_row_split_without_groups():
    matrix = pd.DataFrame({"a": range(10), "b": range(10, 20)})
    scores, audit = cross_fitted_unexplained_outlierness(matrix, None)

    assert audit["status"] == "unavailable"
    assert scores.isna().all()
    assert "random row splitting is forbidden" in audit["reason"]

