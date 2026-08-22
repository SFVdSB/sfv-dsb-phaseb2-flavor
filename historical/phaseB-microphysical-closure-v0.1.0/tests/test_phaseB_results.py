from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def test_summary_claim_and_formula_accuracy():
    s = json.loads((RESULTS / "phaseB_summary.json").read_text())
    assert "post-hoc" in s["status"]
    assert s["benchmark_formula"]["continuous_fit_count"] == 0
    assert s["benchmark_formula"]["max_abs_percent_error"] < 0.61
    assert "not a prediction" in s["benchmark_formula"]["classification"]


def test_locking_relation_across_51_walls():
    s = json.loads((RESULTS / "phaseB_summary.json").read_text())
    r = s["accepted_empirical_relation"]
    assert r["all_51_wall_correlation"] > 0.96
    assert r["all_51_wall_mean_abs_relative_error_pct"] < 0.7
    assert r["all_51_wall_max_abs_relative_error_pct"] < 3.0


def test_constrained_models():
    df = pd.read_csv(RESULTS / "constrained_model_results.csv").set_index("model")
    assert df.loc["locking_interval_fixed", "max_abs_percent_error"] < 0.07
    assert df.loc["four_formula_controls_fixed", "max_abs_percent_error"] < 0.28
    assert df.loc["all_seven_formula_controls_fixed", "max_abs_percent_error"] < 0.61


def test_spectral_chirality():
    df = pd.read_csv(RESULTS / "posthoc_formula_partner_spectrum.csv")
    desired = df[df.operator == "desired"]
    opposite = df[df.operator == "opposite"]
    assert len(desired) == 9
    assert len(opposite) == 9
    assert (desired.near_zero == 1).all()
    assert (opposite.near_zero == 0).all()
    assert opposite.eigenvalue_0.min() > 0.65


def test_balanced_internal_holdouts():
    s = json.loads((RESULTS / "phaseB_summary.json").read_text())
    h = s["internal_balanced_holdout_test"]
    assert h["number_of_calibration_sets"] == 12
    assert h["sets_with_heldout_max_below_1pct"] >= 10
    assert h["median_heldout_max_percent_error"] < 0.5


def test_simple_universal_charge_models_fail():
    df = pd.read_csv(RESULTS / "compact_charge_model_results.csv")
    assert df.max_abs_percent_error.min() > 10.0


def test_all_wall_tables_have_51_cases():
    assert len(pd.read_csv(RESULTS / "refit_controls_all_walls.csv")) == 51
    assert len(pd.read_csv(RESULTS / "formula_closure_all_walls.csv")) == 51
