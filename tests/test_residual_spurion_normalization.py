import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "results/residual_spurion_normalization/residual_spurion_normalization_summary.json"


def load():
    return json.loads(SUMMARY.read_text())


def test_unit_common_is_subpercent():
    s = load()
    assert s["benchmark"]["unit_common_coefficient"]["max_error_pct"] < 1.0


def test_best_common_near_unity_and_subpercent():
    s = load()
    b = s["benchmark"]["best_one_common_coefficient"]
    assert abs(b["lambda_res"] - 1.0) < 0.1
    assert b["max_error_pct"] < 1.0


def test_claim_boundary_not_parameter_free():
    s = load()
    assert "not yet" in s["main_conclusion"].lower()
