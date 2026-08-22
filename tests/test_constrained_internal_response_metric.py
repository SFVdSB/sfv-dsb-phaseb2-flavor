import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/constrained_internal_response_metric"


def summary():
    return json.loads((OUT / "constrained_internal_response_metric_summary.json").read_text())


def test_unique_traceless_generator():
    s = summary()
    assert abs(s["generator"]["trace_T"]) < 1e-14
    assert abs(s["generator"]["trace_T2"] - 4 / 35) < 1e-14


def test_radial_closed_form():
    r = summary()["bounce_radial_driver"]
    assert abs(r["identity_error"]) < 1e-12


def test_exact_action_minimum_and_determinant():
    m = summary()["action_minimum"]["minimum"]
    assert abs(m["A_ln_det_C21"] - m["expected_A"]) < 1e-12
    assert abs(m["B_ln_det_C15"] - m["expected_B"]) < 1e-12
    assert abs(m["det_product"] - 1) < 1e-12


def test_exact_compensator_factors():
    s = summary()
    sigma = s["bounce_radial_driver"]["Sigma_closed_form"]
    m = s["action_minimum"]["minimum"]
    assert abs(m["Z0"] - np.exp(sigma / 21)) < 1e-14
    assert abs(m["ZF_two_legs"] - np.exp(-2 * sigma / 15)) < 1e-14


def test_frozen_flavor_remains_subpercent():
    assert summary()["frozen_flavor_prediction"]["max_error_pct"] < 1.0


def test_finite_heavy_modulus_viable():
    f = summary()["finite_mass_locking"]
    assert f["action_FWHM"]["minimum_m_over_p_below_1pct"] < 2.0
    assert f["gradient_FWHM"]["minimum_m_over_p_below_1pct"] < 2.0


def test_no_extra_dimensions_or_required_particle():
    v = summary()["scientific_verdict"]
    assert not v["extra_spacetime_dimensions_added"]
    assert not v["new_propagating_field_required"]
