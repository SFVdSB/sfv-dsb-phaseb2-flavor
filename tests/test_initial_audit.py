import json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def test_audit_outputs_exist_and_are_finite():
    summary = json.loads((ROOT / "results/initial_audit_summary.json").read_text())
    assert summary["main_positive_results"]["local_kinetic_connection_correlation_with_E"] > 0.97
    assert summary["main_positive_results"]["Hessian_Berry_connection_correlation_with_E"] > 0.90


def test_locking_formula_reproduced():
    d = json.loads((ROOT / "results/berry_locking_diagnostics.json").read_text())
    assert abs(d["alpha_delta_R"] - 0.48790861203019914) < 1e-10
    assert 0.30 < d["integrated_abs_connection_between_landmarks_rad"] < 0.40


def test_formula_basis_dependence_is_exposed():
    df = pd.read_csv(ROOT / "results/frozen_formula_basis_dependence.csv")
    canonical = float(df.loc[df.basis == "canonical_symmetrized_O_E", "max_abs_percent_error"].iloc[0])
    local = float(df.loc[df.basis == "strictly_local_T_raw_G", "max_abs_percent_error"].iloc[0])
    assert canonical < 1.0
    assert local > 10.0
