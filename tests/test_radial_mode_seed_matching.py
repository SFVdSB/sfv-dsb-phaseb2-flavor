import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/radial_mode_seed_matching"


def test_radial_mode_matching_outputs_exist():
    assert (OUT / "radial_mode_seed_matching_summary.json").exists()
    assert (OUT / "radial_identity_all51.csv").exists()


def test_radial_mode_matching_verification():
    tests = json.loads((OUT / "verification_tests.json").read_text())
    assert all(tests.values())


def test_zero_fit_below_one_percent():
    d = json.loads((OUT / "radial_mode_seed_matching_summary.json").read_text())
    assert d["zero_fit_flavor"]["max_error_pct"] < 1.0
    assert abs(d["canonical_anharmonic_matching"]["lambda_radial_sqrt_Vharm_over_Vexact"] - 1.07185503) < 1e-6
