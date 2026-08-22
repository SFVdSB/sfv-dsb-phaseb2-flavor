#!/usr/bin/env python3
"""Canonical radial-mode matching for the residual seed coefficient (v1.4.0).

This module extracts the incomplete-settling mode directly from the corrected
zero-bias two-field bounce.  It tests whether the one remaining residual seed
coefficient can be identified with the finite-amplitude normalization that
converts the true-vacuum harmonic radial coordinate into the exact quartic
bounce displacement.

No flavor observable is used to determine lambda_radial.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from absolute_amplitude_closure import context  # noqa: E402
from residual_spurion_normalization import amplitudes, eval_amplitudes  # noqa: E402

OBS = ["ct", "ut", "sb", "db", "Vus", "Vcb", "Vub"]


def potential_terms(P: float, q: float, rho: float, lam: float,
                    lam_b: float, mu2: float, g: float) -> dict[str, float]:
    return {
        "bulk_quartic": lam / 4.0 * (P * P - rho * rho) ** 2,
        "brane_quartic": lam_b / 4.0 * (q * q - 1.0) ** 2,
        "brane_mass": -mu2 * q * q,
        "portal": g * P * P * q * q,
    }


def main() -> None:
    out = ROOT / "results/radial_mode_seed_matching"
    out.mkdir(parents=True, exist_ok=True)

    profile = pd.read_csv(ROOT / "data/background_profile_O4_regular_robin_full.csv")
    hessian = pd.read_csv(ROOT / "data/baseline_local_hessian_and_mixing.csv")
    all51 = pd.read_csv(
        "/mnt/data/work_phaseA/sfv-dsb-microphysical-dictionary-phaseA-v0.1.0/"
        "results/microphysics_dictionary_all_points.csv"
    )
    old = json.loads(
        (ROOT / "results/residual_spurion_normalization/"
         "residual_spurion_normalization_summary.json").read_text()
    )
    amp_old = json.loads(
        (ROOT / "results/absolute_amplitude_closure/"
         "absolute_amplitude_closure_summary.json").read_text()
    )

    base = all51.loc[all51.name == "baseline_zero_bias"].iloc[0]
    rho = float(base.rho)
    lam = float(base.lambda_Phi)
    lam_b = float(base.lambda_brane)
    mu2 = float(base.mu2_tilde)
    g = float(base.g)
    P0 = float(profile.Phi.iloc[0])
    q0 = float(profile.phi.iloc[0])

    true = np.array([rho, 0.0])
    center = np.array([P0, q0])
    displacement = center - true
    sigma = float(np.linalg.norm(displacement))
    u = float((rho - P0) / rho)

    H_true = np.diag([
        2.0 * lam * rho * rho,
        -lam_b - 2.0 * mu2 + 2.0 * g * rho * rho,
    ])
    eval_true, evec_true = np.linalg.eigh(H_true)
    proj_true = evec_true.T @ displacement
    purity_true = float(proj_true[0] ** 2 / np.dot(displacement, displacement))

    H_center = np.array([
        [float(hessian.H_PhiPhi.iloc[0]), float(hessian.H_cross.iloc[0])],
        [float(hessian.H_cross.iloc[0]), float(hessian.H_phiphi.iloc[0])],
    ])
    eval_center, evec_center = np.linalg.eigh(H_center)
    proj_center = evec_center.T @ displacement
    purity_center = float(proj_center[0] ** 2 / np.dot(displacement, displacement))

    terms = potential_terms(P0, q0, rho, lam, lam_b, mu2, g)
    V_true = lam_b / 4.0
    V_false = lam / 4.0 * rho ** 4 + lam_b - 3.0 * mu2
    deltaV = V_false - V_true
    residual_terms = dict(terms)
    residual_terms["brane_quartic"] -= V_true
    residual_energy = float(sum(residual_terms.values()))
    epsilon_exact = residual_energy / deltaV
    epsilon_bulk_identity = (1.0 - (P0 / rho) ** 2) ** 2

    m_true2 = float(eval_true[0])
    V_harmonic = 0.5 * m_true2 * (rho - P0) ** 2
    anharmonic_energy_ratio = residual_terms["bulk_quartic"] / V_harmonic
    lambda_radial = float(np.sqrt(V_harmonic / residual_terms["bulk_quartic"]))
    lambda_closed = float(1.0 / (1.0 - u / 2.0))

    # Exact identity over all 51 walls.
    all51 = all51.copy()
    all51["epsilon_radial_identity"] = (
        1.0 - all51.center_bulk_completion_fraction ** 2
    ) ** 2
    all51["epsilon_identity_abs_error"] = np.abs(
        all51.epsilon_radial_identity - all51.center_energy_excess_fraction
    )
    all51["u_center"] = 1.0 - all51.center_bulk_completion_fraction
    all51["lambda_radial"] = 1.0 / (1.0 - all51.u_center / 2.0)
    all51[[
        "name", "design_phase", "x", "y", "z",
        "center_bulk_completion_fraction", "center_energy_excess_fraction",
        "epsilon_radial_identity", "epsilon_identity_abs_error",
        "u_center", "lambda_radial",
    ]].to_csv(out / "radial_identity_all51.csv", index=False)

    # Flavor evaluation without determining lambda from flavor data.
    b, H, G, C, gmax, env, target = context()
    inputs = amp_old["inputs"]
    amp = amplitudes(
        float(inputs["center_energy_excess_fraction"]),
        float(inputs["m_true_Phi_dimless"]),
        float(inputs["L_lock"]),
        float(inputs["g_L"]),
        float(inputs["g_R"]),
        lambda_radial,
        lambda_radial,
    )
    flavor = eval_amplitudes(amp, gmax, b, H, C, env, target)

    # Compare only predeclared, bounce-derived normalization choices.
    candidates = {
        "unit_energy_operator": 1.0,
        "canonical_field_to_exact_energy_amplitude": lambda_radial,
        "energy_density_conversion_squared": lambda_radial ** 2,
        "true_to_center_soft_mass_ratio": float(np.sqrt(eval_true[0] / eval_center[0])),
        "sqrt_field_path_tortuosity": float(np.sqrt(base.field_path_tortuosity)),
        "field_path_tortuosity": float(base.field_path_tortuosity),
    }
    rows = []
    for name, value in candidates.items():
        aa = amplitudes(
            float(inputs["center_energy_excess_fraction"]),
            float(inputs["m_true_Phi_dimless"]),
            float(inputs["L_lock"]),
            float(inputs["g_L"]),
            float(inputs["g_R"]),
            value, value,
        )
        rr = eval_amplitudes(aa, gmax, b, H, C, env, target)
        rows.append({
            "normalization": name,
            "lambda_value": value,
            "max_error_pct": rr["max_error_pct"],
            "rms_error_pct": rr["rms_error_pct"],
        })
    pd.DataFrame(rows).to_csv(out / "bounce_normalization_candidates.csv", index=False)

    center_ledger = pd.DataFrame([
        {"component": key, "residual_energy": value,
         "fraction_of_center_residual": value / residual_energy}
        for key, value in residual_terms.items()
    ])
    center_ledger.to_csv(out / "center_residual_energy_decomposition.csv", index=False)

    result = {
        "version": "1.4.0",
        "radial_mode": {
            "true_vacuum": true.tolist(),
            "bounce_center": center.tolist(),
            "canonical_displacement_sigma": sigma,
            "fractional_displacement_u": u,
            "true_hessian_eigenvalues": eval_true.tolist(),
            "center_hessian_eigenvalues": eval_center.tolist(),
            "soft_mode_purity_true_basis": purity_true,
            "soft_mode_purity_center_basis": purity_center,
            "center_residual_energy": residual_energy,
            "false_true_gap": deltaV,
            "epsilon_exact": epsilon_exact,
            "epsilon_bulk_identity": epsilon_bulk_identity,
            "bulk_fraction_of_center_residual": residual_terms["bulk_quartic"] / residual_energy,
        },
        "canonical_anharmonic_matching": {
            "V_harmonic": V_harmonic,
            "V_exact_bulk": residual_terms["bulk_quartic"],
            "V_exact_over_V_harmonic": anharmonic_energy_ratio,
            "lambda_radial_sqrt_Vharm_over_Vexact": lambda_radial,
            "lambda_radial_closed_form": lambda_closed,
            "formula": "lambda_radial = sqrt(V_harmonic/V_exact) = 1/(1-u/2)",
            "interpretation": (
                "Conversion from the true-vacuum harmonic radial amplitude to the exact "
                "finite-displacement quartic radial amplitude. This is bounce-derived and "
                "uses no flavor observable. Its use as the seed Wilson normalization is a "
                "minimal canonical-radial matching principle, not a symmetry theorem."
            ),
        },
        "all51_identity": {
            "max_abs_epsilon_identity_error": float(all51.epsilon_identity_abs_error.max()),
            "mean_abs_epsilon_identity_error": float(all51.epsilon_identity_abs_error.mean()),
            "lambda_radial_min": float(all51.lambda_radial.min()),
            "lambda_radial_max": float(all51.lambda_radial.max()),
            "lambda_radial_mean": float(all51.lambda_radial.mean()),
        },
        "zero_fit_flavor": {
            "lambda_res": lambda_radial,
            "amplitudes": amp,
            "controls": flavor["controls"].tolist(),
            "errors_pct": dict(zip(OBS, flavor["errors"].tolist())),
            "max_error_pct": flavor["max_error_pct"],
            "rms_error_pct": flavor["rms_error_pct"],
        },
        "comparison": {
            "unit_lambda_max_error_pct": old["benchmark"]["unit_common_coefficient"]["max_error_pct"],
            "best_rms_fit_lambda": old["benchmark"]["best_one_common_coefficient"]["lambda_res"],
            "best_rms_fit_max_error_pct": old["benchmark"]["best_one_common_coefficient"]["max_error_pct"],
            "radial_lambda_relative_to_best_fit_pct": 100.0 * (
                lambda_radial / old["benchmark"]["best_one_common_coefficient"]["lambda_res"] - 1.0
            ),
        },
        "claim_boundary": {
            "exactly_derived": [
                "The center displacement is a pure true-bulk soft/radial mode to better than 1e-12 in norm fraction.",
                "epsilon_c=(1-(Phi_c/rho)^2)^2 across all 51 walls to numerical precision.",
                "The finite-amplitude anharmonic conversion lambda_radial=1/(1-u/2).",
            ],
            "conditional": [
                "Identifying lambda_res with lambda_radial requires the seed threshold to be normalized to the canonical radial displacement rather than an independently normalized composite energy operator.",
                "The current symmetries allow an independent O(1) Wilson coefficient, so the matching principle is not uniquely forced by symmetry alone.",
            ],
        },
    }
    (out / "radial_mode_seed_matching_summary.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )

    tests = {
        "center_mode_pure_soft_true_basis": purity_true > 1.0 - 1e-10,
        "center_mode_pure_soft_center_basis": purity_center > 1.0 - 1e-10,
        "center_residual_bulk_dominated": residual_terms["bulk_quartic"] / residual_energy > 1.0 - 1e-9,
        "baseline_radial_energy_identity": abs(epsilon_exact - epsilon_bulk_identity) < 1e-8,
        "all51_radial_energy_identity": float(all51.epsilon_identity_abs_error.max()) < 1e-7,
        "closed_form_matches_energy_ratio": abs(lambda_radial - lambda_closed) < 1e-12,
        "radial_lambda_no_fit_below_1pct": flavor["max_error_pct"] < 1.0,
        "radial_lambda_near_best_common_fit": abs(lambda_radial / old["benchmark"]["best_one_common_coefficient"]["lambda_res"] - 1.0) < 0.01,
    }
    (out / "verification_tests.json").write_text(json.dumps(tests, indent=2) + "\n")
    if not all(tests.values()):
        raise SystemExit("radial matching verification failed")


if __name__ == "__main__":
    main()
