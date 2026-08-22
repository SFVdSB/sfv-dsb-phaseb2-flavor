import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUMMARY = json.loads((ROOT / "results/explicit_radial_uv_matching/explicit_radial_uv_matching_summary.json").read_text())
TESTS = json.loads((ROOT / "results/explicit_radial_uv_matching/verification_tests.json").read_text())


def test_internal_verification():
    assert all(TESTS.values())


def test_uv_response_not_unique():
    assert "independent seed/source functions" in SUMMARY["exact_uv_elimination"]["result"]
    assert SUMMARY["claim_boundary"]["scientific_classification"].startswith("No uniqueness theorem")


def test_radial_matching_remains_subpercent():
    r = SUMMARY["candidate_results"]["canonical_radial_anharmonic_R"]
    assert r["max_error_pct"] < 1.0


def test_unit_lowest_dimension_portal_fails():
    r = SUMMARY["candidate_results"]["lowest_dimension_U1_invariant_chi"]
    assert r["max_error_pct"] > 1.0


def test_exact_radial_identity():
    b = SUMMARY["bounce_invariants"]
    assert abs(b["radial_chi_weight"] - b["two_u"]) < 1e-12
