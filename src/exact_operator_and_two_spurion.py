#!/usr/bin/env python3
"""Exact operator reduction and minimal two-spurion test for SFV/dSB Phase B2.

The script establishes the canonically normalized zero-mode equations for a
4D-Lorentz-symmetric but transversely anisotropic fermion kinetic action,
checks the scalar-Hessian connection as a matrix connection, and evaluates a
predeclared generation-plus-hypercharge two-spurion model.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from chiral_localization import (  # noqa: E402
    CoordinateMap,
    build_wall_basis,
    derive_baseline_normalization,
    flavor_observables,
    higgs_profile,
    load_profile,
    overlap_matrix_trapezoid,
    route1_parameters,
)

OBS = ["ct", "ut", "sb", "db", "Vus", "Vcb", "Vub"]
SECTORS = ["QL", "uR", "dR"]
NGEN = np.array([-1.0, 0.0, 1.0])
HYPERCHARGE = {"QL": 1.0 / 6.0, "uR": 2.0 / 3.0, "dR": -1.0 / 3.0}


def dump_json(path: Path, value: object) -> None:
    def conv(x):
        if isinstance(x, np.ndarray):
            return x.tolist()
        if isinstance(x, (np.floating, np.integer)):
            return x.item()
        if isinstance(x, dict):
            return {str(k): conv(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [conv(v) for v in x]
        return x
    path.write_text(json.dumps(conv(value), indent=2) + "\n", encoding="utf-8")


def vectorized_profiles(x: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Normalized solutions of f'=-Bf for rows of B."""
    S = cumulative_trapezoid(B, x, axis=1, initial=0.0)
    logf = -S
    logf -= np.max(logf, axis=1)[:, None]
    f = np.exp(logf)
    norms = np.sqrt(np.trapezoid(f * f, x, axis=1))
    return f / norms[:, None]


def evaluate_q_kappa(basis, q: dict[str, np.ndarray],
                     kappa: dict[str, np.ndarray], target: np.ndarray):
    profiles = {}
    for sector in SECTORS:
        B = q[sector][:, None] * basis.O[None, :] + kappa[sector][:, None] * basis.E[None, :]
        profiles[sector] = vectorized_profiles(basis.x, B)
    H = higgs_profile(basis.x, s_H=0.514)
    Yu = overlap_matrix_trapezoid(basis.x, profiles["QL"], H, profiles["uR"])
    Yd = overlap_matrix_trapezoid(basis.x, profiles["QL"], H, profiles["dR"])
    result = flavor_observables(Yu, Yd)
    values = np.array([result["values"][key] for key in OBS], dtype=float)
    errors = 100.0 * (values / target - 1.0)
    return values, errors, result


def spurion_coefficients(params: np.ndarray):
    """Generation F=diag(-1,0,1) plus fixed SM-hypercharge spurion.

    q_Ai = exp[a_S Y_A + a_FS Y_A n_i]
    kappa_Ai = b_0 + b_F n_i + b_S Y_A + b_FS Y_A n_i

    The exact operator reduction then gives B_Ai=q_Ai O+kappa_Ai E.
    """
    aS, aFS, b0, bF, bS, bFS = map(float, params)
    q, kappa = {}, {}
    for sector in SECTORS:
        s = HYPERCHARGE[sector]
        q[sector] = np.exp(aS * s + aFS * s * NGEN)
        kappa[sector] = b0 + bF * NGEN + bS * s + bFS * s * NGEN
    return q, kappa


def residual_from_values(values: np.ndarray, target: np.ndarray) -> np.ndarray:
    return np.r_[np.log(values[:4] / target[:4]), (values[4:] - target[4:]) / target[4:]]


def fit_seeded_model(basis, target: np.ndarray, seed: np.ndarray,
                     bounds: tuple[float, float] | None):
    def fun(p):
        try:
            q, kappa = spurion_coefficients(p)
            values, _, _ = evaluate_q_kappa(basis, q, kappa, target)
            r = residual_from_values(values, target)
            return r if np.all(np.isfinite(r)) else np.ones(7) * 100.0
        except Exception:
            return np.ones(7) * 100.0

    kwargs = {}
    if bounds is not None:
        kwargs["bounds"] = bounds
    fit = least_squares(fun, seed, max_nfev=1200, xtol=1e-13, ftol=1e-13,
                        gtol=1e-13, x_scale="jac", **kwargs)
    q, kappa = spurion_coefficients(fit.x)
    values, errors, details = evaluate_q_kappa(basis, q, kappa, target)
    h = {sector: kappa[sector] / q[sector] for sector in SECTORS}
    flat_q = np.concatenate([q[s] for s in SECTORS])
    flat_k = np.concatenate([kappa[s] for s in SECTORS])
    flat_h = np.concatenate([h[s] for s in SECTORS])
    return {
        "parameters": fit.x,
        "values": values,
        "errors_pct": errors,
        "cost_sum_squares": float(np.sum(fit.fun**2)),
        "max_abs_percent_error": float(np.max(np.abs(errors))),
        "rms_percent_error": float(np.sqrt(np.mean(errors**2))),
        "q_min": float(np.min(flat_q)),
        "q_max": float(np.max(flat_q)),
        "kappa_min": float(np.min(flat_k)),
        "kappa_max": float(np.max(flat_k)),
        "h_min": float(np.min(flat_h)),
        "h_max": float(np.max(flat_h)),
        "condition_Yu": details["condition_Yu"],
        "condition_Yd": details["condition_Yd"],
        "q": q,
        "kappa": kappa,
        "h": h,
        "active_bound": bool(bounds is not None and np.any(
            np.isclose(np.abs(fit.x), max(abs(bounds[0]), abs(bounds[1])), rtol=0, atol=1e-6)
        )),
    }


def berry_covariance_test(x: np.ndarray, theta: np.ndarray) -> dict:
    """Verify that U^T U' is required by a local basis rotation.

    This is a covariance test, not a physical flavor fit. A two-channel mass
    matrix is solved in the fixed basis and in its rotating eigenbasis.
    """
    theta_s = CubicSpline(x, theta, bc_type="natural")
    theta_p = theta_s.derivative()
    # Smooth illustrative eigenvalues with a finite channel gap.
    m1 = 0.7 * np.tanh(0.7 * x) + 0.18 * np.exp(-x*x / 2.0)
    m2 = 2.1 + 0.35 * np.tanh(0.5 * x)
    m1s = CubicSpline(x, m1, bc_type="natural")
    m2s = CubicSpline(x, m2, bc_type="natural")

    def U(t):
        th = float(theta_s(t))
        return np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])

    def rhs_fixed(t, f):
        u = U(t)
        w = u @ np.diag([float(m1s(t)), float(m2s(t))]) @ u.T
        return -(w @ f)

    def rhs_rotating(t, g, include_connection=True):
        tp = float(theta_p(t))
        conn = np.array([[0.0, -tp], [tp, 0.0]]) if include_connection else np.zeros((2, 2))
        return -(np.diag([float(m1s(t)), float(m2s(t))]) + conn) @ g

    x0, x1 = float(x[0]), float(x[-1])
    f0 = np.array([1.0, 0.2])
    g0 = U(x0).T @ f0
    sample = np.linspace(x0, x1, 1601)
    sf = solve_ivp(rhs_fixed, (x0, x1), f0, t_eval=sample, method="DOP853",
                   rtol=1e-11, atol=1e-13)
    sg = solve_ivp(lambda t, y: rhs_rotating(t, y, True), (x0, x1), g0,
                   t_eval=sample, method="DOP853", rtol=1e-11, atol=1e-13)
    sg0 = solve_ivp(lambda t, y: rhs_rotating(t, y, False), (x0, x1), g0,
                    t_eval=sample, method="DOP853", rtol=1e-11, atol=1e-13)
    reconstructed = np.column_stack([U(t) @ sg.y[:, j] for j, t in enumerate(sample)])
    omitted = np.column_stack([U(t) @ sg0.y[:, j] for j, t in enumerate(sample)])
    scale = max(float(np.max(np.linalg.norm(sf.y, axis=0))), 1e-30)
    return {
        "rotating_basis_with_connection_relative_max_error": float(
            np.max(np.linalg.norm(reconstructed - sf.y, axis=0)) / scale
        ),
        "rotating_basis_without_connection_relative_max_error": float(
            np.max(np.linalg.norm(omitted - sf.y, axis=0)) / scale
        ),
        "connection_diagonal_is_zero": True,
        "interpretation": (
            "The real Hessian rotation generates an off-diagonal antisymmetric connection. "
            "It is required for basis covariance but is not an independent diagonal scalar potential."
        ),
    }


def main() -> None:
    out = ROOT / "results"
    profile = load_profile(ROOT / "data/background_profile_O4_regular_robin_full.csv")
    cmap = CoordinateMap()
    normalization = derive_baseline_normalization(profile, cmap)
    basis = build_wall_basis(profile, cmap, half_width=24.0, spacing=0.005,
                             normalization=normalization)
    target_json = json.loads((ROOT / "configs/targets_MZ.json").read_text())
    target = np.array([target_json["targets"][key] for key in OBS], dtype=float)

    # Exact canonical reduction.
    controls = json.loads((ROOT / "data/phaseB_formula_controls.json").read_text())["controls"]
    z = np.array([controls[k] for k in
                  ["h_Q", "h_u0", "h_u1", "a_d0", "a_d1", "h_d0", "h_d1"]])
    pars = route1_parameters(z)
    integ = cumulative_trapezoid(basis.E, basis.x, initial=0.0)
    i0 = int(np.argmin(np.abs(basis.x)))
    F_E = integ - integ[i0]
    embedding_rows = []
    max_reconstruction_error = 0.0
    for sector in SECTORS:
        for generation, (q, h) in enumerate(zip(pars[sector]["q"], pars[sector]["h"]), 1):
            kappa = float(q * h)
            A = np.exp(2.0 * kappa * F_E)
            W = A * q * basis.O
            dlogA_half = kappa * basis.E  # exact from definition of F_E
            B_reduced = W / A + dlogA_half
            B_target = q * (basis.O + h * basis.E)
            err = float(np.max(np.abs(B_reduced - B_target)))
            max_reconstruction_error = max(max_reconstruction_error, err)
            embedding_rows.append({
                "sector": sector,
                "generation": generation,
                "q": float(q),
                "h": float(h),
                "kappa_equals_qh": kappa,
                "A_min": float(np.min(A)),
                "A_max": float(np.max(A)),
                "A_max_over_min": float(np.max(A) / np.min(A)),
                "B_reconstruction_max_abs_error": err,
            })
    pd.DataFrame(embedding_rows).to_csv(out / "exact_anisotropic_embedding.csv", index=False)

    # Hessian angle and matrix connection diagnostics.
    hess = pd.read_csv(ROOT / "data/baseline_local_hessian_and_mixing.csv")
    r = hess["r"].to_numpy(float)
    two_theta = np.unwrap(np.arctan2(
        2.0 * hess["H_cross"].to_numpy(float),
        hess["H_PhiPhi"].to_numpy(float) - hess["H_phiphi"].to_numpy(float),
    ))
    theta_r = 0.5 * two_theta
    R_peak = float(profile["R_peak"].iloc[0])
    xh = cmap.alpha * (r - R_peak)
    theta = np.interp(basis.x, xh, theta_r, left=theta_r[0], right=theta_r[-1])
    theta_prime = np.gradient(theta, basis.x)
    soft = np.interp(basis.x, xh, hess["H_eigen_soft"].to_numpy(float))
    hard = np.interp(basis.x, xh, hess["H_eigen_hard"].to_numpy(float))
    gap = np.maximum(hard - soft, 1e-12)
    core = np.abs(basis.x) <= 5.0
    candidates = {
        "abs_theta_prime": np.abs(theta_prime),
        "theta_prime_squared_over_Hessian_gap": theta_prime**2 / gap,
        "abs_theta_prime_over_Hessian_gap": np.abs(theta_prime) / gap,
    }
    berry_rows = []
    for name, raw in candidates.items():
        normalized = raw / max(float(np.max(np.abs(raw[core]))), 1e-30)
        berry_rows.append({
            "candidate": name,
            "correlation_with_canonical_E": float(
                np.corrcoef(normalized[core], basis.E[core])[0, 1]
            ),
            "rms_difference_after_max_normalization": float(
                np.sqrt(np.mean((normalized[core] - basis.E[core])**2))
            ),
        })
    pd.DataFrame(berry_rows).to_csv(out / "berry_matrix_connection_candidates.csv", index=False)
    berry_cov = berry_covariance_test(basis.x[::8], theta[::8])
    dump_json(out / "berry_basis_covariance_test.json", berry_cov)

    # Predeclared two-spurion fits. Seeds are from a coarse global search and
    # are refined here on the full canonical grid.
    models = {
        "hypercharge_two_spurion_natural_bound_3": {
            "seed": np.array([-1.00281251, 1.25924324, 0.32058696,
                              -0.73367663, 0.59328660, 1.09033983]),
            "bounds": (-3.0, 3.0),
            "classification": "O(1)-coefficient test",
        },
        "hypercharge_two_spurion_relaxed_bound_5": {
            "seed": np.array([-2.15670916, -4.27766047, -5.0,
                              -3.64398350, -3.84964931, 5.0]),
            "bounds": (-5.0, 5.0),
            "classification": "relaxed test; coefficients at bounds are a warning",
        },
        "hypercharge_two_spurion_unbounded_seeded": {
            "seed": np.array([-1.81675383, -7.59086627, 0.81916273,
                              -1.28598453, 8.67623366, 6.80101173]),
            "bounds": None,
            "classification": "diagnostic only; rejects naturalness if charges become extreme",
        },
    }
    rows = []
    detailed = {}
    for name, spec in models.items():
        fit = fit_seeded_model(basis, target, spec["seed"], spec["bounds"])
        detailed[name] = fit
        row = {
            "model": name,
            "n_continuous_coefficients": 6,
            "generation_spurion": "diag(-1,0,+1)",
            "sector_spurion": "SM hypercharge (1/6,2/3,-1/3)",
            "classification": spec["classification"],
            "max_abs_percent_error": fit["max_abs_percent_error"],
            "rms_percent_error": fit["rms_percent_error"],
            "q_min": fit["q_min"], "q_max": fit["q_max"],
            "kappa_min": fit["kappa_min"], "kappa_max": fit["kappa_max"],
            "h_min": fit["h_min"], "h_max": fit["h_max"],
            "condition_Yu": fit["condition_Yu"],
            "condition_Yd": fit["condition_Yd"],
            "active_bound": fit["active_bound"],
            "parameters": json.dumps(fit["parameters"].tolist()),
        }
        row.update({f"error_{key}_pct": float(value)
                    for key, value in zip(OBS, fit["errors_pct"])})
        rows.append(row)
    pd.DataFrame(rows).to_csv(out / "two_spurion_hypercharge_results.csv", index=False)
    dump_json(out / "two_spurion_hypercharge_details.json", detailed)

    summary = {
        "version": "0.2.0",
        "status": "exact operator reduction complete; minimal two-spurion test fails naturalness/accuracy",
        "exact_zero_mode_equations": {
            "action": (
                "Hermitian action with 4D kinetic weight Z_parallel, transverse weight Z_perp, "
                "and mass kernel W. After Psi=Z_parallel^{-1/2} chi, define A=Z_perp/Z_parallel "
                "and M=Z_parallel^{-1/2} W Z_parallel^{-1/2}."
            ),
            "left": "f_L' = -A_y f_L -(1/2)(ln A)' f_L - A^{-1} M f_L",
            "right": "f_R' = -A_y f_R -(1/2)(ln A)' f_R + A^{-1} M f_R",
            "sign_note": (
                "Coordinate orientation or the sign assigned to the right-handed wall mass may be reversed; "
                "the invariant statement is that the mass term changes sign between chiralities."
            ),
        },
        "key_correction": (
            "A common isotropic position-dependent kinetic factor Z_parallel=Z_perp cancels under "
            "proper Hermitian canonical normalization. The earlier scalar kinetic-connection interpretation "
            "requires transverse/brane anisotropy or an explicit localized mass operator."
        ),
        "canonical_route1_embedding": {
            "formula": (
                "Let F_E'=E, A_Ai=exp(2 kappa_Ai F_E), W_Ai=A_Ai q_Ai O, "
                "kappa_Ai=q_Ai h_Ai. Then the left-handed reduced equation is exactly "
                "f'=-[q O+kappa E]f."
            ),
            "max_B_reconstruction_abs_error": max_reconstruction_error,
            "largest_transverse_anisotropy_ratio": float(
                max(row["A_max_over_min"] for row in embedding_rows)
            ),
            "warning": (
                "The exact embedding proves operator existence, not naturalness. The Q_L outer generations "
                "require an A_max/A_min ratio about 3.7e4 if the entire E term comes from kinetic anisotropy."
            ),
        },
        "berry_connection": {
            **berry_cov,
            "linear_abs_theta_prime_correlation_with_E": berry_rows[0]["correlation_with_canonical_E"],
            "adiabatic_second_order_theta2_over_gap_correlation_with_E": berry_rows[1]["correlation_with_canonical_E"],
            "conclusion": (
                "The Hessian connection is a real off-diagonal two-channel effect. Treating |theta'| as a "
                "standalone diagonal flavor potential is not justified. A physical contribution requires "
                "channel-dependent masses/couplings and, after elimination, is generally second order."
            ),
        },
        "minimal_two_spurion_test": {
            "model": (
                "F=diag(-1,0,+1), S=SM hypercharge; q=exp(a_S S+a_FS FS), "
                "kappa=b0+bF F+bS S+bFS FS"
            ),
            "n_coefficients": 6,
            "natural_bound_3_max_abs_percent_error": detailed[
                "hypercharge_two_spurion_natural_bound_3"
            ]["max_abs_percent_error"],
            "relaxed_bound_5_max_abs_percent_error": detailed[
                "hypercharge_two_spurion_relaxed_bound_5"
            ]["max_abs_percent_error"],
            "unbounded_seeded_max_abs_percent_error": detailed[
                "hypercharge_two_spurion_unbounded_seeded"
            ]["max_abs_percent_error"],
            "unbounded_seeded_h_max": detailed[
                "hypercharge_two_spurion_unbounded_seeded"
            ]["h_max"],
            "conclusion": (
                "The smallest generation-plus-hypercharge bilinear does not provide a natural closure. "
                "Relaxing it drives coefficients/charges to boundaries or extreme profile parameters."
            ),
        },
        "next_recommended_step": (
            "Keep the exact local operator reduction, but enlarge the flavor algebra minimally: test a second "
            "independent representation spurion or a non-Abelian/discrete flavor generator. In parallel, build "
            "the actual two-channel bulk/brane fermion mass matrix so the Hessian connection contributes "
            "through controlled adiabatic elimination rather than as a guessed scalar profile."
        ),
    }
    dump_json(out / "exact_operator_and_two_spurion_summary.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
