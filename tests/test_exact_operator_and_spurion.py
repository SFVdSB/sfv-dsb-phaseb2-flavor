import json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def test_exact_route1_embedding():
    summary = json.loads((ROOT / "results/exact_operator_and_two_spurion_summary.json").read_text())
    assert summary["canonical_route1_embedding"]["max_B_reconstruction_abs_error"] < 1e-12


def test_isotropic_kernel_correction_and_berry_covariance():
    summary = json.loads((ROOT / "results/exact_operator_and_two_spurion_summary.json").read_text())
    assert "isotropic" in summary["key_correction"]
    berry = summary["berry_connection"]
    assert berry["connection_diagonal_is_zero"] is True
    assert berry["rotating_basis_with_connection_relative_max_error"] < 1e-6
    assert berry["rotating_basis_without_connection_relative_max_error"] > 0.1


def test_minimal_two_spurion_does_not_fake_success():
    df = pd.read_csv(ROOT / "results/two_spurion_hypercharge_results.csv").set_index("model")
    assert df.loc["hypercharge_two_spurion_natural_bound_3", "max_abs_percent_error"] > 10.0
    assert df.loc["hypercharge_two_spurion_relaxed_bound_5", "max_abs_percent_error"] > 5.0
    assert df.loc["hypercharge_two_spurion_unbounded_seeded", "h_max"] > 100.0


def test_anisotropy_warning_is_quantified():
    df = pd.read_csv(ROOT / "results/exact_anisotropic_embedding.csv")
    assert df["B_reconstruction_max_abs_error"].max() < 1e-12
    assert df["A_max_over_min"].max() > 1e4
