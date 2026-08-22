import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/geometric_modulus_embedding"


def load(name):
    return json.loads((OUT / name).read_text())


def test_fluctuation_spectrum_has_one_negative_and_translation_zero():
    s = load("geometric_modulus_embedding_summary.json")
    m = s["o4_fluctuation_spectrum"]
    assert m["ell0_lowest_eigenvalue"] < -0.05
    assert m["ell0_second_eigenvalue"] > 0.0
    assert abs(m["ell1_lowest_eigenvalue"]) < 2e-3
    assert m["ell1_zero_corr_translation"] > 0.99


def test_residual_is_not_identical_to_breathing_mode():
    m = load("geometric_modulus_embedding_summary.json")["o4_fluctuation_spectrum"]
    assert 0.6 < m["residual_breathing_corr_min"]
    assert m["residual_breathing_corr_max"] < 0.95


def test_ordinary_embedding_does_not_force_weights():
    s = load("geometric_modulus_embedding_summary.json")
    assert s["ordinary_embedding_geometry"]["weyl_pair_exact_match_found"] is False
    assert s["scientific_verdict"]["embedding_only"].startswith("fails")


def test_internal_volume_modulus_exactly_reproduces_factors():
    s = load("geometric_modulus_embedding_summary.json")
    req = s["required_compensator"]
    mod = s["minimal_added_symmetry_candidate"]
    assert abs(mod["Z0"] - req["Z0_required"]) < 1e-14
    assert abs(mod["ZF_two_legs"] - req["ZF_required"]) < 1e-14
    assert abs(mod["det_product"] - 1.0) < 1e-12


def test_shell_near_match_is_not_crosswall_law():
    n = load("geometric_modulus_embedding_summary.json")["numerical_shell_coincidence"]
    assert abs(n["relative_to_required_pct"]) < 1.0
    assert n["all51_d3_Rmix_max_abs_relative_error_pct"] > 10.0
