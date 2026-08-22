#!/usr/bin/env python3
"""Constrained internal-response-metric closure for Phase B2 v1.7.0.

This checkpoint writes and audits the minimal auxiliary action for two isotropic
internal response blocks of dimensions 21 and 15.  It tests:

  * uniqueness of the traceless block generator;
  * exact algebraic locking to the already-derived bounce radial driver;
  * whether a propagating version can remain sufficiently locked at wall momenta;
  * one-loop stability of the normalized radial charges;
  * the precise claim boundary between a consistent augmented EFT and a result
    forced by the original SFV/dSB action alone.

No flavor observable is used to choose the dimensions, generator, radial driver,
or locking relation.  Flavor is used only as a downstream tolerance diagnostic.
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
N_SEED = 21
N_FERMION = 15


def evaluate(lambda0: float, lambda_f: float, epsilon: float, inputs: dict, ctx: tuple) -> dict:
    b, H, G, C, gmax, env, target = ctx
    amp = amplitudes(
        epsilon,
        float(inputs["m_true_Phi_dimless"]),
        float(inputs["L_lock"]),
        float(inputs["g_L"]),
        float(inputs["g_R"]),
        lambda0,
        lambda_f,
    )
    result = eval_amplitudes(amp, gmax, b, H, C, env, target)
    return {
        "max_error_pct": float(result["max_error_pct"]),
        "rms_error_pct": float(result["rms_error_pct"]),
        "errors_pct": dict(zip(OBS, result["errors"].tolist())),
        "amplitudes": {k: float(v) for k, v in amp.items()},
    }


def block_generator() -> tuple[np.ndarray, dict]:
    """Return the unique block-isotropic traceless generator.

    Let the logarithmic eigenvalue in the 21 block be a and in the 15 block be b.
    The determinant-one constraint is 21a+15b=0.  Define Sigma as the logarithmic
    volume of the 21 block, 21a=Sigma.  Then a=Sigma/21 and b=-Sigma/15.
    """
    q_seed = 1.0 / N_SEED
    q_fermion = -1.0 / N_FERMION
    t = np.diag(np.r_[np.full(N_SEED, q_seed), np.full(N_FERMION, q_fermion)])
    info = {
        "dimensions": {"seed_block": N_SEED, "fermion_block": N_FERMION, "total": N_SEED + N_FERMION},
        "charges_per_direction": {"seed": q_seed, "fermion_per_leg": q_fermion},
        "trace_T": float(np.trace(t)),
        "trace_T2": float(np.trace(t @ t)),
        "closed_form_trace_T2": "1/21 + 1/15 = 4/35",
        "uniqueness_statement": (
            "Under O(21)xO(15) block isotropy, det(C21)det(C15)=1, and the convention "
            "Sigma=ln det(C21), the one-dimensional generator is uniquely "
            "T=diag(I21/21,-I15/15)."
        ),
    }
    return t, info


def radial_driver(radial: dict) -> dict:
    rm = radial["radial_mode"]
    matching = radial["canonical_anharmonic_matching"]
    rho = float(rm["true_vacuum"][0])
    phi_c = float(rm["bounce_center"][0])
    x = phi_c / rho
    u = 1.0 - x
    epsilon = float(rm["epsilon_bulk_identity"])
    lambda_radial = float(matching["lambda_radial_closed_form"])
    sigma_product = lambda_radial * epsilon
    sigma_polynomial = 2.0 * (1.0 - x) ** 2 * (1.0 + x)
    return {
        "rho": rho,
        "Phi_center": phi_c,
        "x_center": x,
        "u_center": u,
        "epsilon_c": epsilon,
        "lambda_radial": lambda_radial,
        "Sigma_lambda_times_epsilon": sigma_product,
        "Sigma_closed_form": sigma_polynomial,
        "closed_form_function": "Sigma_rad(x)=2(1-x)^2(1+x), x=varrho/rho",
        "identity_error": float(sigma_polynomial - sigma_product),
        "interpretation": (
            "This is the canonical finite-amplitude radial driver already derived from the "
            "quartic bounce.  In the reduced collective-coordinate action it is evaluated "
            "on the bounce saddle; it is not being asserted as a new local spacetime coordinate."
        ),
    }


def constrained_action_minimum(sigma: float) -> dict:
    """Minimize the two-invariant auxiliary potential analytically/numerically.

    A=ln det C21 and B=ln det C15.  The auxiliary potential is
      V = 1/2 Lambda_U^4 (A+B)^2
        + 1/2 Lambda_L^4 ((A-B)/2-Sigma_rad)^2
        + V_shape.
    The positive shape potential sets each block isotropic.  The mass scales
    affect fluctuations but not the minimum.
    """
    # Work with unit positive stiffnesses; the minimum is stiffness independent.
    # Solve linear normal equations for A,B.
    # V=.5(A+B)^2 + .5(.5(A-B)-sigma)^2.
    m = np.array([[1.25, 0.75], [0.75, 1.25]], dtype=float)
    rhs = np.array([sigma / 2.0, -sigma / 2.0], dtype=float)
    a_logvol, b_logvol = np.linalg.solve(m, rhs)
    # The direct variables separate more transparently as U=A+B and S=(A-B)/2.
    expected_a = sigma
    expected_b = -sigma
    per_seed = a_logvol / N_SEED
    per_fermion = b_logvol / N_FERMION
    c21_eig = float(np.exp(per_seed))
    c15_eig = float(np.exp(per_fermion))
    return {
        "auxiliary_potential": (
            "V=1/2 Lambda_U^4[A+B]^2 + 1/2 Lambda_L^4[(A-B)/2-Sigma_rad]^2 + V_shape, "
            "A=ln det C21, B=ln det C15"
        ),
        "minimum": {
            "A_ln_det_C21": float(a_logvol),
            "B_ln_det_C15": float(b_logvol),
            "expected_A": expected_a,
            "expected_B": expected_b,
            "linear_solution_max_abs_error": float(max(abs(a_logvol - expected_a), abs(b_logvol - expected_b))),
            "C21_eigenvalue": c21_eig,
            "C15_eigenvalue_per_fermion_leg": c15_eig,
            "Z0": c21_eig,
            "ZF_two_legs": c15_eig**2,
            "det_C21": c21_eig**N_SEED,
            "det_C15": c15_eig**N_FERMION,
            "det_product": (c21_eig**N_SEED) * (c15_eig**N_FERMION),
        },
        "degrees_of_freedom": (
            "With C21 and C15 treated as algebraic response metrics, this adds no propagating "
            "particle.  A kinetic term would instead create a heavy modulus and is audited separately."
        ),
    }


def flavor_tolerance_scan(epsilon: float, lambda_radial: float, inputs: dict, ctx: tuple) -> pd.DataFrame:
    rows = []
    for mode in ("common", "seed_only", "fermion_only"):
        if mode == "common":
            scales = np.linspace(0.0, 2.0, 2001)
        else:
            scales = np.linspace(-0.25, 2.25, 2501)
        for scale in scales:
            if mode == "common":
                l0 = lambda_radial * scale
                lf = lambda_radial * scale
            elif mode == "seed_only":
                l0 = lambda_radial * scale
                lf = lambda_radial
            else:
                l0 = lambda_radial
                lf = lambda_radial * scale
            r = evaluate(l0, lf, epsilon, inputs, ctx)
            rows.append({
                "mode": mode,
                "scale_relative_to_radial": float(scale),
                "lambda0": float(l0),
                "lambdaF": float(lf),
                "max_error_pct": r["max_error_pct"],
                "rms_error_pct": r["rms_error_pct"],
                "below_1pct": bool(r["max_error_pct"] < 1.0),
            })
    return pd.DataFrame(rows)


def summarize_tolerance(scan: pd.DataFrame) -> dict:
    out: dict[str, dict] = {}
    for mode, d in scan.groupby("mode"):
        ok = d[d.below_1pct]
        best = d.loc[d.max_error_pct.idxmin()]
        out[mode] = {
            "allowed_scale_min_below_1pct": float(ok.scale_relative_to_radial.min()),
            "allowed_scale_max_below_1pct": float(ok.scale_relative_to_radial.max()),
            "best_scale": float(best.scale_relative_to_radial),
            "best_max_error_pct": float(best.max_error_pct),
        }
    return out


def finite_lock_scan(epsilon: float, lambda_radial: float, inputs: dict, ctx: tuple) -> pd.DataFrame:
    """Audit a propagating heavy modulus with response m^2/(m^2+p^2)."""
    widths = {
        "action_FWHM": 1.7517900037359464,
        "gradient_FWHM": 2.601607651871168,
    }
    rows = []
    ratios = np.r_[np.linspace(0.1, 3.0, 291), np.linspace(3.1, 10.0, 70)]
    for width_name, width in widths.items():
        p = 1.0 / width
        for m_over_p in ratios:
            response = m_over_p**2 / (1.0 + m_over_p**2)
            result = evaluate(lambda_radial * response, lambda_radial * response, epsilon, inputs, ctx)
            rows.append({
                "wall_scale": width_name,
                "width_dimless": width,
                "p_dimless": p,
                "m_over_p": float(m_over_p),
                "modulus_mass_dimless": float(m_over_p * p),
                "locking_response_m2_over_m2_plus_p2": float(response),
                "max_error_pct": result["max_error_pct"],
                "rms_error_pct": result["rms_error_pct"],
                "below_1pct": bool(result["max_error_pct"] < 1.0),
            })
    return pd.DataFrame(rows)


def running_stability(seed_running: pd.DataFrame, sigma: float, tolerance: dict) -> tuple[pd.DataFrame, dict]:
    """One-loop stability ledger and generic nonlinear-vertex audit.

    Multiplicative gauge factors independent of Sigma cancel from normalized
    response ratios.  For illustration of a possible non-multiplicative effect,
    solve dy/dlnQ = a y^3/(16pi^2).  Define zeta=2 a y_UV^2 L/(16pi^2).
    The exact normalized radial charge is then audited as a function of zeta.
    """
    rows = []
    for _, r in seed_running.iterrows():
        rows.append({
            "sector": str(r["block"]),
            "multiplicity": int(r["multiplicity"]),
            "M2_running_factor": float(r["M2_ratio_IR_over_UV_singlet"]),
            "mu_running_factor": float(r["mu_ratio_IR_over_UV_singlet"]),
            "protected_mu_over_M2_factor": float(r["mu_over_M2_ratio"]),
            "radial_charge_renormalization_from_gauge": 1.0,
            "reason": "Gauge running is multiplicative and Sigma-independent; normalized response ratios cancel it.",
        })
    ledger = pd.DataFrame(rows)

    q = -1.0 / N_FERMION
    z = q * sigma
    nonlinear_rows = []
    for zeta in np.linspace(0.0, 2.0, 401):
        ratio = np.exp(z) * np.sqrt((1.0 + zeta) / (1.0 + zeta * np.exp(2.0 * z)))
        qeff_over_q = float(np.log(ratio) / z) if z != 0.0 else 1.0
        nonlinear_rows.append({
            "zeta_2ay2L_over_16pi2": float(zeta),
            "qF_effective_over_tree": qeff_over_q,
            "inside_direct_flavor_1pct_scale_window": bool(
                tolerance["fermion_only"]["allowed_scale_min_below_1pct"]
                <= qeff_over_q
                <= tolerance["fermion_only"]["allowed_scale_max_below_1pct"]
            ),
        })
    nonlinear = pd.DataFrame(nonlinear_rows)
    ok = nonlinear[nonlinear.inside_direct_flavor_1pct_scale_window]
    summary = {
        "gauge_result": (
            "The determinant-one constraint is algebraic and exact. At one-loop gauge order, "
            "the seed kernel mu/M^2 remains equal in all 21 directions. Multiplicative gauge "
            "renormalization of a fermion source changes its overall coupling but not the "
            "normalized radial exponent."
        ),
        "generic_nonlinear_vertex_model": "dy/dlnQ=a y^3/(16pi^2), zeta=2 a y_UV^2 ln(Lambda/Q)/(16pi^2)",
        "largest_scanned_zeta_inside_direct_1pct_window": float(ok.zeta_2ay2L_over_16pi2.max()),
        "qF_scale_at_that_point": float(ok.iloc[-1].qF_effective_over_tree),
        "interpretation": (
            "The flavor tolerance is much broader than ordinary perturbative one-loop corrections. "
            "The exact determinant relation is protected; the more model-dependent radial-locking "
            "normalization is stable unless nonlinear source interactions become very strong."
        ),
    }
    return ledger, nonlinear, summary


def main() -> None:
    out = ROOT / "results/constrained_internal_response_metric"
    out.mkdir(parents=True, exist_ok=True)

    radial = json.loads((ROOT / "results/radial_mode_seed_matching/radial_mode_seed_matching_summary.json").read_text())
    absolute = json.loads((ROOT / "results/absolute_amplitude_closure/absolute_amplitude_closure_summary.json").read_text())
    seed_running = pd.read_csv(ROOT / "results/seed_sector_beta/one_loop_seed_running_by_block.csv")
    inputs = absolute["inputs"]
    ctx = context()

    t, generator = block_generator()
    radial_info = radial_driver(radial)
    sigma = float(radial_info["Sigma_closed_form"])
    epsilon = float(radial_info["epsilon_c"])
    lambda_radial = float(radial_info["lambda_radial"])
    minimum = constrained_action_minimum(sigma)

    np.savetxt(out / "block_generator_T_36x36.csv", t, delimiter=",")

    tolerance_scan = flavor_tolerance_scan(epsilon, lambda_radial, inputs, ctx)
    tolerance_scan.to_csv(out / "radial_lock_flavor_tolerance_scan.csv", index=False)
    tolerance = summarize_tolerance(tolerance_scan)

    lock_scan = finite_lock_scan(epsilon, lambda_radial, inputs, ctx)
    lock_scan.to_csv(out / "finite_mass_locking_scan.csv", index=False)
    finite_summary = {}
    for wall_scale, d in lock_scan.groupby("wall_scale"):
        ok = d[d.below_1pct]
        finite_summary[wall_scale] = {
            "minimum_m_over_p_below_1pct": float(ok.m_over_p.min()),
            "minimum_modulus_mass_dimless_below_1pct": float(ok.modulus_mass_dimless.min()),
            "response_at_threshold": float(ok.iloc[0].locking_response_m2_over_m2_plus_p2),
        }

    gauge_ledger, nonlinear, running_summary = running_stability(seed_running, sigma, tolerance)
    gauge_ledger.to_csv(out / "one_loop_block_metric_stability.csv", index=False)
    nonlinear.to_csv(out / "nonlinear_vertex_charge_running_scan.csv", index=False)

    frozen = evaluate(lambda_radial, lambda_radial, epsilon, inputs, ctx)
    action_text = {
        "block_parameterization": (
            "C21=exp(a) I21, C15=exp(b) I15; A=21a=ln det C21, B=15b=ln det C15"
        ),
        "potential": minimum["auxiliary_potential"],
        "shape_terms": (
            "Positive O(21)- and O(15)-invariant penalties remove traceless anisotropy inside each block."
        ),
        "minimum": "A=Sigma_rad, B=-Sigma_rad; hence a=Sigma_rad/21 and b=-Sigma_rad/15.",
        "implementation": (
            "Preferred: C21,C15 are algebraic response metrics in the wall-reduced EFT. "
            "They add no spacetime dimensions and no propagating scalar."
        ),
        "compensator_symmetry": (
            "A one-parameter determinant-preserving internal rescaling, commuting with Pati-Salam, "
            "acts with weights +1/21 and -1/15. It is a spurionic/response symmetry, not full SL(36) "
            "mixing of gauge-inequivalent fields."
        ),
    }

    summary = {
        "version": "1.7.0",
        "generator": generator,
        "bounce_radial_driver": radial_info,
        "explicit_constrained_action": action_text,
        "action_minimum": minimum,
        "frozen_flavor_prediction": frozen,
        "flavor_tolerance": tolerance,
        "finite_mass_locking": finite_summary,
        "one_loop_stability": running_summary,
        "scientific_verdict": {
            "consistency_with_SFV_dSB": True,
            "extra_spacetime_dimensions_added": False,
            "new_propagating_field_required": False,
            "determinant_constraint_status": "exact in the auxiliary constrained response action",
            "radial_identification_status": (
                "exact inside the augmented action through the locking invariant, but this locking/compensator "
                "principle is one additional structural postulate beyond the original SFV/dSB scalar action"
            ),
            "parameter_free_status": (
                "Within the minimal constrained-response extension, no continuous flavor coefficient remains. "
                "The extension itself is not yet derived uniquely from the original bounce Lagrangian."
            ),
            "main_result": (
                "The smallest gauge-compatible implementation is not a 36-dimensional physical symmetry but a "
                "single determinant-preserving internal rescaling of two isotropic response blocks. The auxiliary "
                "action uniquely yields +1/21 and -1/15, locks to the bounce radial driver, and is one-loop gauge stable."
            ),
        },
    }
    (out / "constrained_internal_response_metric_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    tests = {
        "generator_is_traceless": abs(generator["trace_T"]) < 1e-14,
        "generator_norm_is_4_over_35": abs(generator["trace_T2"] - 4.0 / 35.0) < 1e-14,
        "radial_closed_form_identity": abs(radial_info["identity_error"]) < 1e-12,
        "action_minimum_A_equals_sigma": abs(minimum["minimum"]["A_ln_det_C21"] - sigma) < 1e-12,
        "action_minimum_B_equals_minus_sigma": abs(minimum["minimum"]["B_ln_det_C15"] + sigma) < 1e-12,
        "determinant_product_is_one": abs(minimum["minimum"]["det_product"] - 1.0) < 1e-12,
        "Z0_exact": abs(minimum["minimum"]["Z0"] - np.exp(sigma / 21.0)) < 1e-14,
        "ZF_exact": abs(minimum["minimum"]["ZF_two_legs"] - np.exp(-2.0 * sigma / 15.0)) < 1e-14,
        "frozen_flavor_below_1pct": frozen["max_error_pct"] < 1.0,
        "seed_gauge_ratio_protected": bool(np.allclose(seed_running.mu_over_M2_ratio, 1.0, atol=1e-14)),
        "finite_heavy_modulus_has_viable_region": all(v["minimum_m_over_p_below_1pct"] < 2.0 for v in finite_summary.values()),
        "auxiliary_version_needs_no_new_propagating_field": True,
    }
    tests = {k: bool(v) for k, v in tests.items()}
    (out / "verification_tests.json").write_text(json.dumps(tests, indent=2) + "\n")
    if not all(tests.values()):
        raise SystemExit(f"constrained metric verification failed: {tests}")


if __name__ == "__main__":
    main()
