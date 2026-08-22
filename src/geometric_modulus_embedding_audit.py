#!/usr/bin/env python3
"""Geometric-modulus and brane-breathing audit for Phase B2 v1.6.0.

This checkpoint tests whether ordinary codimension-one brane embedding geometry
can force the residual compensator factors

    Z0 = exp(+Sigma/21),   ZF = exp(-2 Sigma/15)

with Sigma=lambda_radial*epsilon_c, without an additional internal symmetry.

The calculation has four parts:
  1. Solve the two-field O(4) fluctuation operator and verify the genuine
     negative breathing mode and l=1 translation zero mode.
  2. Derive universal Weyl/canonical-normalization weights for a uniform
     normal displacement of a d-dimensional wall.
  3. Compare natural geometric drivers with the bounce-derived Sigma at the
     benchmark and across all 51 walls.
  4. Exhibit the minimal internal block-volume (unimodular) modulus that does
     generate 1/21 and -2/15 exactly, while keeping its claim boundary clear.

No flavor observable is used to set a geometric coefficient. Flavor is used
only to propagate predeclared candidate drivers through the frozen model.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.sparse import bmat, diags
from scipy.sparse.linalg import eigsh

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from absolute_amplitude_closure import context  # noqa: E402
from residual_spurion_normalization import amplitudes, eval_amplitudes  # noqa: E402

OBS = ["ct", "ut", "sb", "db", "Vus", "Vcb", "Vub"]


def normalized(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n == 0.0:
        raise ValueError("zero vector")
    return v / n


def fluctuation_modes(profile: pd.DataFrame, n_grid: int = 2200,
                      r_max: float = 30.0) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    """Solve the coupled O(4) radial fluctuation operator for ell=0,1.

    With q=r^{-3/2}u, the radial operator is
      -d2/dr2 + H(r) + [ell(ell+2)+3/4]/r2.
    The ell=0 negative state is the bubble breathing mode.  The ell=1
    translation state should be numerically near zero and align with Phi'.
    """
    r = np.linspace(0.002, r_max, n_grid)
    P = np.interp(r, profile.r, profile.Phi)
    q = np.interp(r, profile.r, profile.phi)
    Pp = np.interp(r, profile.r, profile.Phi_prime)
    qp = np.interp(r, profile.r, profile.phi_prime)

    lam_phi = 0.1
    rho = 2.357142857142857
    bias = 0.0
    portal = 2.313019
    lam_b = 1.0e-8
    mu2 = 1.0e-8

    h11 = lam_phi * (3.0 * P * P - rho * rho) + 2.0 * bias + 2.0 * portal * q * q
    h22 = lam_b * (3.0 * q * q - 1.0) - 2.0 * mu2 + 2.0 * portal * P * P
    h12 = 4.0 * portal * P * q

    ri = r[1:-1]
    n = len(ri)
    dr = float(r[1] - r[0])
    T = diags(
        [np.full(n - 1, -1.0 / dr**2), np.full(n, 2.0 / dr**2),
         np.full(n - 1, -1.0 / dr**2)],
        [-1, 0, 1], format="csr",
    )

    mode_rows: list[dict] = []
    vectors: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for ell in (0, 1):
        centrifugal = (ell * (ell + 2.0) + 0.75) / (ri * ri)
        A = T + diags(h11[1:-1] + centrifugal)
        D = T + diags(h22[1:-1] + centrifugal)
        B = diags(h12[1:-1])
        L = bmat([[A, B], [B, D]], format="csr")
        vals, vecs = eigsh(L, k=8, which="SA", tol=1e-9, maxiter=250000)
        order = np.argsort(vals)
        vals = vals[order]
        vecs = vecs[:, order]
        vectors[ell] = (vals, vecs)

        translation = normalized(np.concatenate([
            ri**1.5 * Pp[1:-1], ri**1.5 * qp[1:-1]
        ]))
        dilation = normalized(np.concatenate([
            ri**2.5 * Pp[1:-1], ri**2.5 * qp[1:-1]
        ]))
        for j in range(8):
            v = normalized(vecs[:, j])
            mode_rows.append({
                "ell": ell,
                "mode_index": j,
                "eigenvalue": float(vals[j]),
                "abs_corr_translation_tangent": float(abs(v @ translation)),
                "abs_corr_dilation_tangent": float(abs(v @ dilation)),
            })

    # Audit how much a natural interior-settling deformation resembles the
    # ell=0 breathing mode.  The window choice is varied and reported.
    vals0, vecs0 = vectors[0]
    neg = normalized(vecs0[:, 0])
    residual_rows = []
    radii = {
        "gradient_peak": 5.207425309387891,
        "action_peak": 5.860202508437851,
        "mixing_peak": 5.913579732488854,
    }
    for radius_name, radius in radii.items():
        for width in (0.15, 0.30, 0.60, 1.00):
            window = 1.0 / (1.0 + np.exp(np.clip((ri - radius) / width, -60.0, 60.0)))
            residual = normalized(np.concatenate([
                ri**1.5 * (rho - P[1:-1]) * window,
                ri**1.5 * (-q[1:-1]) * window,
            ]))
            residual_rows.append({
                "window_center": radius_name,
                "window_width": width,
                "abs_corr_with_ell0_negative_breathing_mode": float(abs(neg @ residual)),
            })

    summary = {
        "ell0_lowest_eigenvalue": float(vectors[0][0][0]),
        "ell0_second_eigenvalue": float(vectors[0][0][1]),
        "ell1_lowest_eigenvalue": float(vectors[1][0][0]),
        "ell0_negative_corr_translation": float(
            pd.DataFrame(mode_rows).query("ell == 0 and mode_index == 0")[
                "abs_corr_translation_tangent"
            ].iloc[0]
        ),
        "ell0_negative_corr_dilation": float(
            pd.DataFrame(mode_rows).query("ell == 0 and mode_index == 0")[
                "abs_corr_dilation_tangent"
            ].iloc[0]
        ),
        "ell1_zero_corr_translation": float(
            pd.DataFrame(mode_rows).query("ell == 1 and mode_index == 0")[
                "abs_corr_translation_tangent"
            ].iloc[0]
        ),
        "residual_breathing_corr_min": float(
            min(row["abs_corr_with_ell0_negative_breathing_mode"] for row in residual_rows)
        ),
        "residual_breathing_corr_max": float(
            max(row["abs_corr_with_ell0_negative_breathing_mode"] for row in residual_rows)
        ),
    }
    return pd.DataFrame(mode_rows), summary, pd.DataFrame(residual_rows)


def weyl_weight_table(sigma: float) -> pd.DataFrame:
    """Canonical-normalization weights under gamma->e^(2 omega) gamma.

    X is a canonical scalar seed.  S is either nondynamical (auxiliary) or a
    canonical propagating scalar.  The protected seed coefficient is mu/M2;
    the fermion factor is the square of the source vertex y.
    """
    desired = np.array([sigma / 21.0, -2.0 * sigma / 15.0])
    rows = []
    for d in (3, 4):
        cases = {
            "auxiliary_S": np.array([0.0, 2.0]),
            "canonical_scalar_S": np.array([(2.0 - d) / 2.0, 4.0 - d]),
        }
        for name, weights in cases.items():
            denom = float(weights @ weights)
            omega_best = float((weights @ desired) / denom) if denom > 0 else 0.0
            pred = weights * omega_best
            residual = pred - desired
            rows.append({
                "worldvolume_dimension_d": d,
                "mediator_treatment": name,
                "dlnZ0_domega": float(weights[0]),
                "dlnZF_domega": float(weights[1]),
                "best_omega_for_required_pair": omega_best,
                "predicted_lnZ0": float(pred[0]),
                "predicted_lnZF": float(pred[1]),
                "required_lnZ0": float(desired[0]),
                "required_lnZF": float(desired[1]),
                "relative_L2_mismatch": float(np.linalg.norm(residual) / np.linalg.norm(desired)),
                "exact_pair_match": bool(np.linalg.norm(residual) < 1e-12),
            })
    return pd.DataFrame(rows)


def evaluate_driver(driver: float, eps: float, inputs: dict, ctx: tuple) -> dict:
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
    return {
        "lambda_equivalent": float(lam),
        "max_error_pct": float(result["max_error_pct"]),
        "rms_error_pct": float(result["rms_error_pct"]),
        "errors_pct": dict(zip(OBS, result["errors"].tolist())),
    }


def geometric_driver_table(radial: dict, absolute: dict, phase_a: pd.DataFrame,
                           profile: pd.DataFrame) -> pd.DataFrame:
    eps = float(radial["radial_mode"]["epsilon_bulk_identity"])
    u = float(radial["radial_mode"]["fractional_displacement_u"])
    sigma_field = float(radial["radial_mode"]["canonical_displacement_sigma"])
    m_true = float(absolute["inputs"]["m_true_Phi_dimless"])
    required = float(radial["canonical_anharmonic_matching"]["lambda_radial_closed_form"]) * eps
    base = phase_a[phase_a.name == "baseline_zero_bias"].iloc[0]
    radii = {
        "action_peak": float(base.R_peak_dimless),
        "gradient_peak": float(base.R_gradient_peak_dimless),
        "mixing_peak": float(base.hessian_mixing_max_radius_dimless),
    }
    speed = np.sqrt(profile.Phi_prime**2 + profile.phi_prime**2)
    speed_at_mix = float(np.interp(radii["mixing_peak"], profile.r, speed))
    displacements = {
        "field_fraction_u_as_coordinate_shift_unproven": u,
        "soft_correlation_length_u_over_m": u / m_true,
        "level_set_shift_sigma_over_speed_at_mix": sigma_field / speed_at_mix,
    }
    ctx = context()
    rows = []
    for d in (3, 4):
        for rname, radius in radii.items():
            for shift_name, delta_r in displacements.items():
                driver = d * np.log1p(delta_r / radius)
                flavor = evaluate_driver(float(driver), eps, absolute["inputs"], ctx)
                rows.append({
                    "worldvolume_or_shell_dimension": d,
                    "radius_definition": rname,
                    "normal_shift_definition": shift_name,
                    "radius": radius,
                    "delta_r": float(delta_r),
                    "geometric_driver_d_log_1_plus_deltaR_over_R": float(driver),
                    "required_radial_driver": required,
                    "relative_driver_error_pct": float(100.0 * (driver / required - 1.0)),
                    "lambda_equivalent": flavor["lambda_equivalent"],
                    "max_flavor_error_pct": flavor["max_error_pct"],
                    "rms_flavor_error_pct": flavor["rms_error_pct"],
                    "claim_status": (
                        "diagnostic only: field-space displacement is not a derived coordinate shift"
                        if shift_name == "field_fraction_u_as_coordinate_shift_unproven"
                        else "bounce-derived scale mapping but not a symmetry derivation"
                    ),
                })
    return pd.DataFrame(rows)


def crosswall_shell_audit(radial51: pd.DataFrame, phase_a: pd.DataFrame) -> pd.DataFrame:
    merged = radial51.merge(
        phase_a[["name", "R_peak_dimless", "hessian_mixing_max_radius_dimless"]],
        on="name", how="inner",
    )
    merged["required_driver"] = merged.lambda_radial * merged.epsilon_radial_identity
    for d in (3, 4):
        for tag, col in (
            ("action_peak", "R_peak_dimless"),
            ("mixing_peak", "hessian_mixing_max_radius_dimless"),
        ):
            cand = d * np.log1p(merged.u_center / merged[col])
            merged[f"driver_d{d}_{tag}"] = cand
            merged[f"relative_error_pct_d{d}_{tag}"] = 100.0 * (cand / merged.required_driver - 1.0)
    return merged


def internal_volume_modulus(sigma: float) -> dict:
    """Exact algebra of the minimal O(21)xO(15) block-volume compensator."""
    seed_per_direction = float(np.exp(sigma / 21.0))
    fermion_per_leg = float(np.exp(-sigma / 15.0))
    zf = fermion_per_leg**2
    det_seed_response = seed_per_direction**21
    det_fermion_leg_response = fermion_per_leg**15
    return {
        "ansatz": {
            "C21": "exp(Sigma/21) I_21",
            "C15_per_fermion_leg": "exp(-Sigma/15) I_15",
            "unimodular_constraint": "det(C21) det(C15)=1",
        },
        "Sigma": sigma,
        "Z0": seed_per_direction,
        "fermion_leg_factor": fermion_per_leg,
        "ZF_two_legs": zf,
        "det_C21": det_seed_response,
        "det_C15": det_fermion_leg_response,
        "det_product": det_seed_response * det_fermion_leg_response,
        "interpretation": (
            "O(21) and O(15) isotropy distribute one logarithmic block-volume modulus "
            "equally over 21 seed directions and 15 SU(4) directions. Two fermion legs "
            "give the factor -2/15. This is an added internal response-space symmetry, "
            "not a consequence of ordinary brane embedding geometry."
        ),
    }


def main() -> None:
    out = ROOT / "results/geometric_modulus_embedding"
    out.mkdir(parents=True, exist_ok=True)

    profile = pd.read_csv(ROOT / "data/background_profile_O4_regular_robin_full.csv")
    phase_a = pd.read_csv(ROOT / "data/phaseA_amplitude_invariants_all51.csv")
    radial51 = pd.read_csv(ROOT / "results/radial_mode_seed_matching/radial_identity_all51.csv")
    radial = json.loads((ROOT / "results/radial_mode_seed_matching/radial_mode_seed_matching_summary.json").read_text())
    absolute = json.loads((ROOT / "results/absolute_amplitude_closure/absolute_amplitude_closure_summary.json").read_text())

    eps = float(radial["radial_mode"]["epsilon_bulk_identity"])
    lambda_rad = float(radial["canonical_anharmonic_matching"]["lambda_radial_closed_form"])
    sigma_required = lambda_rad * eps

    modes, mode_summary, residual_scan = fluctuation_modes(profile)
    modes.to_csv(out / "o4_fluctuation_modes.csv", index=False)
    residual_scan.to_csv(out / "breathing_residual_overlap_window_scan.csv", index=False)

    weyl = weyl_weight_table(sigma_required)
    weyl.to_csv(out / "brane_breathing_weyl_weight_audit.csv", index=False)

    geom = geometric_driver_table(radial, absolute, phase_a, profile)
    geom.to_csv(out / "geometric_driver_candidates.csv", index=False)

    cross = crosswall_shell_audit(radial51, phase_a)
    cross.to_csv(out / "geometric_shell_driver_all51.csv", index=False)

    internal = internal_volume_modulus(sigma_required)
    (out / "internal_block_volume_modulus.json").write_text(json.dumps(internal, indent=2) + "\n")

    # Select benchmark diagnostics highlighted in the report.
    s3_mix = geom[
        (geom.worldvolume_or_shell_dimension == 3)
        & (geom.radius_definition == "mixing_peak")
        & (geom.normal_shift_definition == "field_fraction_u_as_coordinate_shift_unproven")
    ].iloc[0]
    d4_mix = geom[
        (geom.worldvolume_or_shell_dimension == 4)
        & (geom.radius_definition == "mixing_peak")
        & (geom.normal_shift_definition == "field_fraction_u_as_coordinate_shift_unproven")
    ].iloc[0]

    rel3 = cross["relative_error_pct_d3_mixing_peak"]
    rel3_local = cross[cross.design_phase != "corridor"]["relative_error_pct_d3_mixing_peak"]

    summary = {
        "version": "1.6.0",
        "required_compensator": {
            "epsilon_c": eps,
            "lambda_radial": lambda_rad,
            "Sigma_required": sigma_required,
            "Z0_required": float(np.exp(sigma_required / 21.0)),
            "ZF_required": float(np.exp(-2.0 * sigma_required / 15.0)),
        },
        "o4_fluctuation_spectrum": mode_summary,
        "ordinary_embedding_geometry": {
            "uniform_breathing_metric": "gamma_mn -> exp(2 omega) gamma_mn",
            "protected_auxiliary_seed_ratio": "mu/M^2 is invariant under the universal measure/canonical rescaling",
            "trace_average_result": (
                "For O_N=(1/N) sum_A X_A^2, a common geometric rescaling of every X_A "
                "rescales O_N by the common field weight, not by 1/N. Normalized traces do "
                "not generate 1/21 or 1/15 Weyl charges."
            ),
            "weyl_pair_exact_match_found": bool(weyl.exact_pair_match.any()),
            "conclusion": (
                "Ordinary normal-displacement/extrinsic-curvature geometry supplies a real "
                "breathing mode but does not force the opposite representation-normalized "
                "weights +1/21 and -2/15. Curvature portals can reproduce them only with "
                "additional sector-specific coefficients."
            ),
        },
        "numerical_shell_coincidence": {
            "O4_S3_shell_driver_using_deltaR_equals_u_at_Rmix": float(
                s3_mix.geometric_driver_d_log_1_plus_deltaR_over_R
            ),
            "relative_to_required_pct": float(s3_mix.relative_driver_error_pct),
            "frozen_flavor_max_error_pct": float(s3_mix.max_flavor_error_pct),
            "physical_3plus1_worldvolume_d4_analogue_driver": float(
                d4_mix.geometric_driver_d_log_1_plus_deltaR_over_R
            ),
            "d4_relative_to_required_pct": float(d4_mix.relative_driver_error_pct),
            "d4_frozen_flavor_max_error_pct": float(d4_mix.max_flavor_error_pct),
            "all51_d3_Rmix_mean_abs_relative_error_pct": float(np.mean(np.abs(rel3))),
            "all51_d3_Rmix_max_abs_relative_error_pct": float(np.max(np.abs(rel3))),
            "local33_d3_Rmix_mean_abs_relative_error_pct": float(np.mean(np.abs(rel3_local))),
            "local33_d3_Rmix_max_abs_relative_error_pct": float(np.max(np.abs(rel3_local))),
            "claim_boundary": (
                "The S3 benchmark near-match is numerically interesting but assumes without "
                "derivation that the field-space fraction u is a coordinate displacement. "
                "It is not stable enough across the 51 walls to establish a geometric law."
            ),
        },
        "minimal_added_symmetry_candidate": internal,
        "scientific_verdict": {
            "embedding_only": "fails to derive the compensator charges",
            "breathing_mode_exists": True,
            "residual_has_substantial_breathing_overlap_but_is_not_identical": True,
            "internal_unimodular_block_volume_symmetry": (
                "exactly generates +1/21 and -2/15 without compact spacetime dimensions"
            ),
            "remaining_condition": (
                "The internal block-volume modulus must be induced/locked to the existing "
                "bounce radial driver Sigma=lambda_radial epsilon_c. That identification is "
                "the added geometric-compensator principle; ordinary brane embedding does not force it."
            ),
        },
    }
    (out / "geometric_modulus_embedding_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    tests = {
        "one_ell0_negative_mode": mode_summary["ell0_lowest_eigenvalue"] < -0.05
        and mode_summary["ell0_second_eigenvalue"] > 0.0,
        "ell1_translation_near_zero": abs(mode_summary["ell1_lowest_eigenvalue"]) < 2.0e-3,
        "ell1_translation_profile_verified": mode_summary["ell1_zero_corr_translation"] > 0.99,
        "ell0_is_breathing_like": mode_summary["ell0_negative_corr_translation"] > 0.90,
        "residual_not_identical_to_breathing": mode_summary["residual_breathing_corr_max"] < 0.95,
        "ordinary_weyl_cases_do_not_exactly_match_pair": not bool(weyl.exact_pair_match.any()),
        "internal_volume_modulus_exact_Z0": abs(internal["Z0"] - np.exp(sigma_required / 21.0)) < 1e-14,
        "internal_volume_modulus_exact_ZF": abs(internal["ZF_two_legs"] - np.exp(-2.0 * sigma_required / 15.0)) < 1e-14,
        "internal_volume_constraint": abs(internal["det_product"] - 1.0) < 1e-12,
        "s3_benchmark_near_match_but_not_declared_derivation": abs(float(s3_mix.relative_driver_error_pct)) < 1.0,
        "s3_crosswall_not_universal": float(np.max(np.abs(rel3))) > 10.0,
    }
    tests = {k: bool(v) for k, v in tests.items()}
    (out / "verification_tests.json").write_text(json.dumps(tests, indent=2) + "\n")
    if not all(tests.values()):
        raise SystemExit(f"geometric modulus verification failed: {tests}")


if __name__ == "__main__":
    main()
