import math

from darkdna.utils.stats import finite_mean, optional_float, optional_row_float, row_score
from darkdna.views.boundary_conditions import compute_boundary_condition_view
from darkdna.views.entropy_noise import compute_entropy_noise_view
from darkdna.views.physical_susceptibility import compute_physical_susceptibility_view


def test_optional_float_preserves_measured_zero_and_keeps_missing_as_na():
    assert optional_float(0.0) == 0.0
    assert optional_float("0") == 0.0
    assert math.isnan(optional_float(None))
    assert math.isnan(optional_float(""))
    assert math.isnan(optional_float(float("nan")))
    assert math.isnan(optional_row_float({}, "entropy_boundary_score"))
    assert optional_row_float({"entropy_boundary_score": 0.0}, "entropy_boundary_score") == 0.0
    assert math.isnan(row_score({}, ["a", "b"]))
    assert row_score({"a": 2.0}, ["a", "b"]) == 2.0
    assert math.isnan(finite_mean([float("nan"), float("nan")]))


def test_entropy_and_boundary_views_do_not_convert_missing_inputs_to_zero():
    empty = compute_entropy_noise_view({})
    assert math.isnan(empty["decoherence_boundary_candidate_score"])
    assert math.isnan(empty["entropy_cliff_score"])

    measured_zero = compute_entropy_noise_view({"entropy_boundary_score": 0.0})
    assert measured_zero["entropy_cliff_score"] == 0.0
    assert measured_zero["decoherence_boundary_candidate_score"] == 0.0
    assert math.isnan(measured_zero["feature_void_score"])

    boundary = compute_boundary_condition_view({"entropy_boundary_score": 1.5})
    assert boundary["entropy_transition_score"] == 1.5
    assert boundary["boundary_condition_candidate_score"] == 1.5
    assert math.isnan(boundary["sequence_regime_boundary_score"])
    assert math.isnan(compute_boundary_condition_view({})["boundary_condition_candidate_score"])


def test_physical_susceptibility_view_omits_missing_predictors_from_aggregates():
    empty = compute_physical_susceptibility_view({})
    assert math.isnan(empty["fork_texture_score"])
    assert math.isnan(empty["charge_oxidation_susceptibility_score"])
    assert math.isnan(empty["nonB_physical_susceptibility_score"])
    assert math.isnan(empty["physical_view_G4_susceptibility"])

    g4_only = compute_physical_susceptibility_view({"G4_sequence_potential": 0.8})
    assert g4_only["physical_view_G4_susceptibility"] == 0.8
    assert g4_only["charge_oxidation_susceptibility_score"] == 0.8
    assert g4_only["nonB_physical_susceptibility_score"] == 0.8
    assert math.isnan(g4_only["physical_view_Z_DNA_sequence_potential"])
    assert math.isnan(g4_only["fork_texture_score"])
