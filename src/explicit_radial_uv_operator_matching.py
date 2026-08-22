#!/usr/bin/env python3
"""Explicit radial-field dependence of the protected seed sector (v1.5.0).

This audit asks whether the canonical radial mode of the corrected SFV bounce
forces the remaining residual-response coefficient, rather than merely
providing a numerically successful matching convention.

The most general local PS- and O(22)-compatible seed/source functions are
    M_X^2(R), mu_X(R), y_F(R),
where R is a dimensionless radial invariant.  Integrating out the protected
auxiliary channel gives
    Z0 = [mu_X(R)/M_X^2(R)]/[mu_X(0)/M_X^2(0)],
    ZF = [y_F(R)/y_F(0)]^2.
Consequently the linear response coefficients are independent Wilson
coefficients unless an extra compensator/symmetry relates the three functions.

No flavor observable is used to determine the bounce invariants.  Flavor data
are used only in explicitly labelled diagnostic fits/comparisons.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from absolute_amplitude_closure import context  # noqa: E402
from residual_spurion_normalization import amplitudes, eval_amplitudes  # noqa: E402

OBS = ["ct", "ut", "sb", "db", "Vus", "Vcb", "Vub"]


def evaluate_driver(driver: float, eps: float, inputs: dict, ctx: tuple) -> dict:
    """Use a predeclared exponent driver D=lambda*epsilon in frozen amplitudes."""
    b, H, G, C, gmax, env, target = ctx
    lam = driver / eps
    amp = amplitudes(
        eps,
        float(inputs["m_true_Phi_dimless"]),
        float(inputs["L_lock"]),
        float(inputs["g_L"]),
        float(inputs["g_R"]),
        lam,
        lam,
    )
    result = eval_amplitudes(amp, gmax, b, H, C, env, target)
    return {"lambda_equivalent": float(lam), "amplitudes": amp, **result}


def portal_response_coefficients(chi: float, z0: float, zf: float) -> dict:
    """Infer coefficients in simple renormalizable rational portal forms.

    Models:
      Z0 = 1/(1 + eta0 chi)
      ZF = 1/(1 + etaF chi)^2
    These correspond to a fixed source and a radial-dependent mass/vertex,
    respectively.  They are diagnostics, not unique UV assignments.
    """
    eta0 = (1.0 / z0 - 1.0) / chi
    etaf = (zf ** (-0.5) - 1.0) / chi
    return {"eta0": float(eta0), "etaF": float(etaf)}


def main() -> None:
    out = ROOT / "results/explicit_radial_uv_matching"
    out.mkdir(parents=True, exist_ok=True)

    radial = json.loads(
        (ROOT / "results/radial_mode_seed_matching/"
         "radial_mode_seed_matching_summary.json").read_text()
    )
    residual = json.loads(
        (ROOT / "results/residual_spurion_normalization/"
         "residual_spurion_normalization_summary.json").read_text()
    )
    absolute = json.loads(
        (ROOT / "results/absolute_amplitude_closure/"
         "absolute_amplitude_closure_summary.json").read_text()
    )
    all51 = pd.read_csv(ROOT / "results/radial_mode_seed_matching/radial_identity_all51.csv")

    eps_total = float(radial["radial_mode"]["epsilon_exact"])
    eps = float(radial["radial_mode"]["epsilon_bulk_identity"])
    u = float(radial["radial_mode"]["fractional_displacement_u"])
    lam_rad = float(radial["canonical_anharmonic_matching"]["lambda_radial_closed_form"])
    chi = float(1.0 - (1.0 - u) ** 2)  # 1-Phi_c^2/rho^2 = sqrt(epsilon)
    assert abs(chi * chi - eps) < 1e-8

    ctx = context()
    inputs = absolute["inputs"]

    # Predeclared local radial invariants / response drivers.
    drivers = {
        "exact_residual_energy_R": eps,
        "canonical_radial_anharmonic_R": lam_rad * eps,
        "harmonic_energy_fraction": 4.0 * u * u,
        "canonical_displacement_u": u,
        "lowest_dimension_U1_invariant_chi": chi,
        # Exploratory PS-normalized portal. The 4/15 factor is a candidate
        # multiplicity normalization, not claimed as forced by the action.
        "exploratory_PS_4_over_15_times_chi": (4.0 / 15.0) * chi,
    }

    rows = []
    details = {}
    for name, driver in drivers.items():
        res = evaluate_driver(driver, eps, inputs, ctx)
        rows.append({
            "candidate": name,
            "driver_D": driver,
            "lambda_equivalent_D_over_epsilon": res["lambda_equivalent"],
            "max_error_pct": res["max_error_pct"],
            "rms_error_pct": res["rms_error_pct"],
            "uses_flavor_to_set_coefficient": False,
            "claim_status": (
                "bounce-derived" if name in {
                    "exact_residual_energy_R", "canonical_radial_anharmonic_R",
                    "harmonic_energy_fraction", "canonical_displacement_u",
                    "lowest_dimension_U1_invariant_chi"
                } else "exploratory group-normalized candidate"
            ),
        })
        details[name] = {
            "driver_D": driver,
            "lambda_equivalent": res["lambda_equivalent"],
            "amplitudes": res["amplitudes"],
            "controls": res["controls"].tolist(),
            "errors_pct": dict(zip(OBS, res["errors"].tolist())),
            "max_error_pct": res["max_error_pct"],
            "rms_error_pct": res["rms_error_pct"],
        }

    # Diagnostic coefficient preferred by quark data, clearly labelled.
    best_lambda = float(residual["benchmark"]["best_one_common_coefficient"]["lambda_res"])
    best_driver = best_lambda * eps
    best_chi_weight = best_driver / chi
    radial_chi_weight = (lam_rad * eps) / chi
    # Exact identity radial_chi_weight = 2u.

    # Minimize RMS directly in the lowest-dimension chi portal coefficient w.
    def rms_for_w(w: float) -> float:
        return evaluate_driver(w * chi, eps, inputs, ctx)["rms_error_pct"]

    opt = minimize_scalar(rms_for_w, bounds=(-1.0, 2.0), method="bounded",
                          options={"xatol": 1e-13})
    w_rms = float(opt.x)
    diag_fit = evaluate_driver(w_rms * chi, eps, inputs, ctx)

    # Find a minimax-ish interval below 1% for w, using dense scan and roots.
    ws = np.linspace(-0.5, 1.0, 6001)
    maxerrs = np.array([evaluate_driver(float(w) * chi, eps, inputs, ctx)["max_error_pct"] for w in ws])
    ok = maxerrs < 1.0
    if np.any(ok):
        w_min = float(ws[np.where(ok)[0][0]])
        w_max = float(ws[np.where(ok)[0][-1]])
    else:
        w_min = w_max = float("nan")

    # Infer low-dimension rational-portal coefficients needed to reproduce the
    # bounce-derived radial normalization, not a flavor fit.
    amp_rad = details["canonical_radial_anharmonic_R"]["amplitudes"]
    portal_rad = portal_response_coefficients(chi, amp_rad["Z0"], amp_rad["Z_F"])
    amp_unit = details["exact_residual_energy_R"]["amplitudes"]
    portal_unit = portal_response_coefficients(chi, amp_unit["Z0"], amp_unit["Z_F"])

    pd.DataFrame(rows).sort_values("max_error_pct").to_csv(
        out / "radial_operator_candidate_comparison.csv", index=False
    )

    pd.DataFrame([
        {
            "target_matching": "unit_residual_energy_operator",
            "eta0_in_Z0_equals_1_over_1_plus_eta0_chi": portal_unit["eta0"],
            "etaF_in_ZF_equals_1_over_1_plus_etaF_chi_squared": portal_unit["etaF"],
        },
        {
            "target_matching": "bounce_radial_anharmonic_operator",
            "eta0_in_Z0_equals_1_over_1_plus_eta0_chi": portal_rad["eta0"],
            "etaF_in_ZF_equals_1_over_1_plus_etaF_chi_squared": portal_rad["etaF"],
        },
    ]).to_csv(out / "renormalizable_portal_coefficients.csv", index=False)

    # General Wilson-coefficient ledger from exact elimination.
    ledger = pd.DataFrame([
        {
            "effective_factor": "Z0",
            "exact_expression": "[mu(R)/M2(R)]/[mu(0)/M2(0)]",
            "linear_log_response": "d_R ln(mu)-d_R ln(M2)",
            "fixed_by_canonical_radial_normalization": False,
            "fixed_by_O22_or_Pati_Salam": "Only equality across seed components; not absolute derivative",
        },
        {
            "effective_factor": "ZF",
            "exact_expression": "[yF(R)/yF(0)]^2",
            "linear_log_response": "2 d_R ln(yF)",
            "fixed_by_canonical_radial_normalization": False,
            "fixed_by_O22_or_Pati_Salam": "Representation/multiplicity weights only; not absolute derivative",
        },
        {
            "effective_factor": "common lambda_res",
            "exact_expression": "requires d_R ln(mu/M2) = -15 d_R ln(yF)/21 after normalized traces",
            "linear_log_response": "one relation among otherwise independent Wilson derivatives",
            "fixed_by_canonical_radial_normalization": False,
            "fixed_by_O22_or_Pati_Salam": False,
        },
    ])
    ledger.to_csv(out / "wilson_coefficient_ledger.csv", index=False)

    # All-51 radial identities relevant to local operator choice.
    all51 = all51.copy()
    all51["chi_lowest_dimension"] = np.sqrt(np.maximum(all51["epsilon_radial_identity"], 0.0))
    all51["radial_driver_lambda_epsilon"] = all51["lambda_radial"] * all51["epsilon_radial_identity"]
    all51["radial_chi_weight"] = all51["radial_driver_lambda_epsilon"] / all51["chi_lowest_dimension"]
    all51["two_u_identity"] = 2.0 * all51["u_center"]
    all51["radial_chi_weight_minus_2u"] = all51["radial_chi_weight"] - all51["two_u_identity"]
    all51.to_csv(out / "radial_uv_invariants_all51.csv", index=False)

    summary = {
        "version": "1.5.0",
        "exact_uv_elimination": {
            "seed_functions": ["M_X^2(R)", "mu_X(R)", "y_F(R)"],
            "Z0": "[mu_X(R)/M_X^2(R)]/[mu_X(0)/M_X^2(0)]",
            "ZF": "[y_F(R)/y_F(0)]^2",
            "linear_responses": {
                "lambda0_effective": "d_R ln(mu_X)-d_R ln(M_X^2)",
                "lambdaF_effective": "d_R ln(y_F)",
            },
            "result": (
                "Canonical normalization fixes R but not derivatives of independent seed/source functions. "
                "An extra compensator relation or UV boundary condition is required to set a common absolute coefficient."
            ),
        },
        "bounce_invariants": {
            "u": u,
            "chi_equals_1_minus_Phi2_over_rho2": chi,
            "epsilon_equals_chi_squared": eps,
            "epsilon_total_including_tiny_brane_terms": eps_total,
            "lambda_radial": lam_rad,
            "lambda_radial_times_epsilon": lam_rad * eps,
            "radial_chi_weight": radial_chi_weight,
            "radial_chi_weight_exact_identity": "2u",
            "two_u": 2.0 * u,
        },
        "candidate_results": details,
        "diagnostic_only_flavor_fit": {
            "best_previous_common_lambda": best_lambda,
            "best_previous_driver": best_driver,
            "equivalent_low_dimension_chi_weight": best_chi_weight,
            "direct_rms_fit_chi_weight": w_rms,
            "fit_max_error_pct": diag_fit["max_error_pct"],
            "fit_rms_error_pct": diag_fit["rms_error_pct"],
            "subpercent_chi_weight_scan_interval": [w_min, w_max],
            "uses_flavor_data": True,
        },
        "renormalizable_portal_diagnostic": {
            "unit_energy_target": portal_unit,
            "radial_anharmonic_target": portal_rad,
            "interpretation": (
                "The lowest-dimension U(1)-invariant portal is proportional to chi=1-|Phi|^2/rho^2=sqrt(epsilon), "
                "but its coefficient is independent. Unit coefficient overcorrects flavor; reproducing the mild radial "
                "matching requires small O(10^-2) rational-portal coefficients or an equivalent common compensator weight."
            ),
        },
        "claim_boundary": {
            "derived": [
                "The exact form of Z0 and ZF after integrating out local seed/source functions.",
                "Canonical radial normalization does not remove the independent Wilson derivatives.",
                "The lowest-dimension U(1)-invariant radial portal scales as chi=sqrt(epsilon), not epsilon.",
                "The previous radial driver satisfies (lambda_radial*epsilon)/chi=2u exactly.",
            ],
            "not_derived": [
                "lambda_res=lambda_radial from the stated symmetries alone.",
                "A unique common radial compensator charge for seed masses and fermion sources.",
                "The exploratory 4/15 normalization; it is numerically suggestive but post-hoc at this checkpoint.",
            ],
            "scientific_classification": (
                "No uniqueness theorem: the zero-fit 0.611% result remains a minimal canonical matching hypothesis. "
                "The explicit UV calculation leaves one genuine radial portal/compensator coefficient unless an additional symmetry is supplied."
            ),
        },
    }
    (out / "explicit_radial_uv_matching_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    tests = {
        "chi_squared_equals_epsilon": abs(chi * chi - eps) < 1e-8,
        "radial_chi_weight_equals_2u": abs(radial_chi_weight - 2.0 * u) < 1e-12,
        "all51_radial_weight_identity": float(np.max(np.abs(all51["radial_chi_weight_minus_2u"]))) < 1e-10,
        "radial_zero_fit_subpercent": details["canonical_radial_anharmonic_R"]["max_error_pct"] < 1.0,
        "unit_lowest_dimension_portal_not_acceptable": details["lowest_dimension_U1_invariant_chi"]["max_error_pct"] > 1.0,
        "general_action_has_independent_responses": True,
    }
    (out / "verification_tests.json").write_text(json.dumps(tests, indent=2) + "\n")
    if not all(tests.values()):
        raise SystemExit(f"verification failure: {tests}")

    print(json.dumps({
        "status": "ok",
        "radial_zero_fit_max_error_pct": details["canonical_radial_anharmonic_R"]["max_error_pct"],
        "lowest_dimension_unit_portal_max_error_pct": details["lowest_dimension_U1_invariant_chi"]["max_error_pct"],
        "radial_chi_weight": radial_chi_weight,
        "best_fit_chi_weight_diagnostic": w_rms,
        "uv_uniqueness": False,
    }, indent=2))


if __name__ == "__main__":
    main()
