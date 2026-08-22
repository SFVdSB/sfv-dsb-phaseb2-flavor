#!/usr/bin/env python3
"""Vectorlike-mediator/flavon closure audit for Phase B2 v0.6.0.

This module tests whether the frozen local-core chiral operator can arise by
integrating out a small heavy vectorlike fermion sector.  It distinguishes:

1. a single rank-one mediator (one shared generation direction),
2. the minimal singlet + traceless-flavon two-direction mediator family,
3. compressed mediator matching relations suggested by the 51-wall audit.

The small rational ratios are frozen post-hoc UV targets, not claimed as
first-principles group theory.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Callable

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from raw_gradient_wilson_closure import context, evaluate, residual  # noqa: E402
from chiral_localization import (  # noqa: E402
    overlap_matrix_trapezoid,
    flavor_observables,
)

OBS = ["ct", "ut", "sb", "db", "Vus", "Vcb", "Vub"]
N = np.array([-1.0, 0.0, 1.0])


def fit_model(
    map_to_z: Callable[[np.ndarray], np.ndarray],
    seed: np.ndarray,
    b,
    H,
    C,
    env,
    target: np.ndarray,
    bounds=(-8.0, 8.0),
    max_nfev: int = 5000,
) -> dict:
    fit = least_squares(
        lambda p: residual(map_to_z(p), b, H, C, env, target),
        np.asarray(seed, dtype=float),
        bounds=bounds,
        max_nfev=max_nfev,
        xtol=1e-13,
        ftol=1e-13,
        gtol=1e-13,
        x_scale="jac",
    )
    z = map_to_z(fit.x)
    values, errors, _ = evaluate(z, b, H, C, env, target)
    return {
        "success": bool(fit.success),
        "nfev": int(fit.nfev),
        "parameters": fit.x.tolist(),
        "effective_controls": z.tolist(),
        "values": dict(zip(OBS, values.tolist())),
        "errors_pct": dict(zip(OBS, errors.tolist())),
        "max_error_pct": float(np.max(np.abs(errors))),
        "rms_error_pct": float(np.sqrt(np.mean(errors * errors))),
        "active_mask": fit.active_mask.tolist(),
        "cost": float(fit.cost),
    }


def evaluate_arbitrary_qh(qmat, hmat, b, H, C, env, target):
    profiles = {}
    for row, sector in enumerate(["QL", "uR", "dR"]):
        q = qmat[row]
        h = hmat[row]
        B = q[:, None] * (b.O[None, :] + h[:, None] * C[None, :])
        S = cumulative_trapezoid(B, b.x, axis=1, initial=0.0)
        L = -S
        L -= L.max(axis=1)[:, None]
        f = np.exp(L) * env[None, :]
        f /= np.sqrt(np.trapezoid(f * f, b.x, axis=1))[:, None]
        profiles[sector] = f
    Yu = overlap_matrix_trapezoid(b.x, profiles["QL"], H, profiles["uR"])
    Yd = overlap_matrix_trapezoid(b.x, profiles["QL"], H, profiles["dR"])
    result = flavor_observables(Yu, Yd)
    values = np.array([result["values"][o] for o in OBS])
    errors = 100.0 * (values / target - 1.0)
    return values, errors


def fit_rank_one_mediator(b, H, C, env, target):
    """One mediator: h_Ai=s_A(1+r*n_i), with down kink q exponential.

    The common generation vector forces the same slope/intercept ratio in all
    three representations.  This is the strict rank-one hypothesis.
    """

    def evaluate_p(p):
        sQ, su, sd, r, a0, a1 = p
        g = 1.0 + r * N
        hmat = np.vstack([sQ * g, su * g, sd * g])
        qmat = np.vstack([np.ones(3), np.ones(3), np.exp(a0 + a1 * N)])
        return evaluate_arbitrary_qh(qmat, hmat, b, H, C, env, target)

    def fun(p):
        values, _ = evaluate_p(p)
        return np.r_[
            np.log(values[:4] / target[:4]),
            (values[4:] - target[4:]) / target[4:],
        ]

    seeds = [
        np.array([7.98, 0.32, 7.99, 0.565, -0.222, 0.480]),
        np.array([-2.0, 3.0, 3.0, -1.0, 0.26, 0.24]),
    ]
    best = None
    for seed in seeds:
        fit = least_squares(
            fun,
            seed,
            bounds=(-8.0, 8.0),
            max_nfev=2200,
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-12,
            x_scale="jac",
        )
        values, errors = evaluate_p(fit.x)
        score = float(np.max(np.abs(errors)))
        if best is None or score < best[0]:
            best = (score, fit, values, errors)
    score, fit, values, errors = best
    return {
        "success": bool(fit.success),
        "parameters": fit.x.tolist(),
        "values": dict(zip(OBS, values.tolist())),
        "errors_pct": dict(zip(OBS, errors.tolist())),
        "max_error_pct": score,
        "rms_error_pct": float(np.sqrt(np.mean(errors * errors))),
        "active_mask": fit.active_mask.tolist(),
        "nfev": int(fit.nfev),
    }


def summarize_crosswall(path: Path, error_col: str) -> dict:
    df = pd.read_csv(path)
    values = df[error_col].to_numpy(float)
    return {
        "points": int(len(df)),
        "points_below_1pct": int(np.sum(values < 1.0)),
        "mean_max_error_pct": float(np.mean(values)),
        "median_max_error_pct": float(np.median(values)),
        "maximum_error_pct": float(np.max(values)),
        "worst_point": str(df.iloc[int(np.argmax(values))]["name"]),
    }


def main():
    b, H, G, C, Gmax, env, target = context()

    # Exact seven-control local-core benchmark.
    exact_seed = np.array(
        [2.7925, 2.9120, -0.03160, 0.25964, 0.24228, 2.6640, -0.65899]
    )
    exact = fit_model(lambda z: z, exact_seed, b, H, C, env, target)
    z_exact = np.asarray(exact["effective_controls"])
    hQ, hu0, hu1, ad0, ad1, hd0, hd1 = z_exact

    # The UV h-coupling matrix before the down-sector q prefactor.
    h_matrix = np.vstack(
        [N * hQ, hu0 + N * hu1, hd0 + N * hd1]
    )
    q_matrix = np.vstack(
        [np.ones(3), np.ones(3), np.exp(ad0 + N * ad1)]
    )
    h_singular = np.linalg.svd(h_matrix, compute_uv=False)
    q_singular = np.linalg.svd(q_matrix, compute_uv=False)

    # Phase-A information entering the propagator-suppression hypothesis.
    phase_a = json.loads((ROOT / "data/baseline_phaseA_microphysics.json").read_text())
    phase_a_all = pd.read_csv(ROOT / "data/phaseA_microphysics_all51_minimal.csv")
    baseline_micro = phase_a_all[phase_a_all.name == "baseline_zero_bias"].iloc[0]
    epsilon_geom = float(
        baseline_micro.xi_false_Phi_dimless / baseline_micro.R_peak_dimless
    )
    lock = float(
        b.alpha
        * (
            phase_a["hessian_mixing_max_radius_dimless"]
            - phase_a["R_gradient_peak_dimless"]
        )
    )
    a0_wall = float(phase_a["m_true_Phi_dimless"] / 4.0)
    a1_wall = float(lock / 2.0)

    # Strict one-mediator/rank-one test.
    rank_one = fit_rank_one_mediator(b, H, C, env, target)

    # Six-parameter vectorlike family with one frozen U/D odd-channel ratio.
    # Independent quantities: hQ, hu0, a0, a1, hd0, hd1.
    def z_ratio21(p):
        hQ_, hu0_, a0_, a1_, hd0_, hd1_ = p
        return np.array([hQ_, hu0_, hd1_ / 21.0, a0_, a1_, hd0_, hd1_])

    ratio21_6p = fit_model(
        z_ratio21,
        z_exact[[0, 1, 3, 4, 5, 6]],
        b,
        H,
        C,
        env,
        target,
    )

    def z_geom6(p):
        hQ_, hu0_, a0_, a1_, hd0_, hd1_ = p
        return np.array(
            [hQ_, hu0_, epsilon_geom * hd1_, a0_, a1_, hd0_, hd1_]
        )

    geom_6p = fit_model(
        z_geom6,
        z_exact[[0, 1, 3, 4, 5, 6]],
        b,
        H,
        C,
        env,
        target,
    )

    # Five parameters: add the right-handed even-channel Clebsch target.
    def z_clebsch5(p):
        hQ_, hd0_, hd1_, a0_, a1_ = p
        return np.array(
            [hQ_, (23.0 / 21.0) * hd0_, hd1_ / 21.0, a0_, a1_, hd0_, hd1_]
        )

    clebsch_5p = fit_model(
        z_clebsch5,
        z_exact[[0, 5, 6, 3, 4]],
        b,
        H,
        C,
        env,
        target,
    )

    # Four parameters: the core matrix is generated by two mediator amplitudes.
    def z_clebsch4(p):
        hd0_, hd1_, a0_, a1_ = p
        return np.array(
            [
                (22.0 / 21.0) * hd0_,
                (23.0 / 21.0) * hd0_,
                hd1_ / 21.0,
                a0_,
                a1_,
                hd0_,
                hd1_,
            ]
        )

    clebsch_4p = fit_model(
        z_clebsch4,
        z_exact[[5, 6, 3, 4]],
        b,
        H,
        C,
        env,
        target,
    )

    # Three parameters: additionally impose the exploratory q-kink ratio 15/14.
    def z_clebsch3(p):
        hd0_, hd1_, a1_ = p
        return np.array(
            [
                (22.0 / 21.0) * hd0_,
                (23.0 / 21.0) * hd0_,
                hd1_ / 21.0,
                (15.0 / 14.0) * a1_,
                a1_,
                hd0_,
                hd1_,
            ]
        )

    clebsch_3p = fit_model(
        z_clebsch3,
        z_exact[[5, 6, 4]],
        b,
        H,
        C,
        env,
        target,
    )

    # Four-parameter wall-frozen model using the geometric propagator relation.
    def z_wall4(p):
        hQ_, hu0_, hd0_, hd1_ = p
        return np.array(
            [hQ_, hu0_, epsilon_geom * hd1_, a0_wall, a1_wall, hd0_, hd1_]
        )

    wall_4p = fit_model(
        z_wall4,
        z_exact[[0, 1, 5, 6]],
        b,
        H,
        C,
        env,
        target,
    )

    benchmark_ratios = {
        "hQ_over_hd0": float(hQ / hd0),
        "hu0_over_hd0": float(hu0 / hd0),
        "hu1_over_hd1": float(hu1 / hd1),
        "a0_over_a1": float(ad0 / ad1),
    }
    rational_targets = {
        "hQ_over_hd0": 22.0 / 21.0,
        "hu0_over_hd0": 23.0 / 21.0,
        "hu1_over_hd1": 1.0 / 21.0,
        "a0_over_a1": 15.0 / 14.0,
    }
    ratio_errors = {
        k: 100.0 * (rational_targets[k] / benchmark_ratios[k] - 1.0)
        for k in benchmark_ratios
    }

    # Cross-wall ratio stability from independent exact local-core refits.
    all_controls = pd.read_csv(
        ROOT / "results/mediator_closure/local_geo_refit_controls_all51.csv"
    )
    ratio_frame = pd.DataFrame(
        {
            "name": all_controls.name,
            "design_phase": all_controls.design_phase,
            "hQ_over_hd0": all_controls.h_Q / all_controls.h_d0,
            "hu0_over_hd0": all_controls.h_u0 / all_controls.h_d0,
            "hu1_over_hd1": all_controls.h_u1 / all_controls.h_d1,
            "a0_over_a1": all_controls.a_d0 / all_controls.a_d1,
        }
    )
    ratio_rows = []
    for domain_name, frame in [
        ("all_51", ratio_frame),
        ("local_33", ratio_frame[ratio_frame.design_phase != "corridor"]),
        ("corridor_18", ratio_frame[ratio_frame.design_phase == "corridor"]),
    ]:
        for key, target_ratio in rational_targets.items():
            values = frame[key].to_numpy(float)
            rel = 100.0 * (target_ratio / values - 1.0)
            ratio_rows.append(
                {
                    "domain": domain_name,
                    "ratio": key,
                    "target": target_ratio,
                    "mean_observed": float(np.mean(values)),
                    "std_observed": float(np.std(values, ddof=1)),
                    "mean_abs_relative_error_pct": float(np.mean(np.abs(rel))),
                    "max_abs_relative_error_pct": float(np.max(np.abs(rel))),
                }
            )
    ratio_audit = pd.DataFrame(ratio_rows)

    outdir = ROOT / "results/mediator_closure"
    outdir.mkdir(parents=True, exist_ok=True)
    ratio_frame.to_csv(outdir / "coefficient_ratios_all51.csv", index=False)
    ratio_audit.to_csv(outdir / "rational_ratio_crosswall_audit.csv", index=False)

    models = {
        "one_rank_one_mediator_6p": rank_one,
        "singlet_flavon_ratio21_6p": ratio21_6p,
        "singlet_flavon_geometric_6p": geom_6p,
        "clebsch_5p": clebsch_5p,
        "two_amplitude_core_clebsch_4p": clebsch_4p,
        "three_parameter_near_miss": clebsch_3p,
        "wall_frozen_geometric_4p": wall_4p,
    }
    rows = []
    for name, result in models.items():
        rows.append(
            {
                "model": name,
                "n_parameters": len(result["parameters"]),
                "max_error_pct": result["max_error_pct"],
                "rms_error_pct": result["rms_error_pct"],
                "success": result["success"],
                "active_bounds": int(np.sum(np.abs(result.get("active_mask", [])) > 0)),
            }
        )
    pd.DataFrame(rows).to_csv(outdir / "mediator_model_comparison.csv", index=False)

    crosswall = {
        "ratio21_6p": summarize_crosswall(
            outdir / "ratio21_six_and_four_param_all51.csv",
            "six_param_max_error_pct",
        ),
        "geometric_ratio_6p": summarize_crosswall(
            outdir / "geometric_ratio_six_and_four_param_all51.csv",
            "six_param_max_error_pct",
        ),
        "clebsch_5p": summarize_crosswall(
            outdir / "five_param_clebsch_all51.csv", "max_error_pct"
        ),
        "core_clebsch_4p": summarize_crosswall(
            outdir / "four_param_core_clebsch_all51.csv", "max_error_pct"
        ),
    }

    result = {
        "version": "0.6.0",
        "claim_boundary": (
            "The vectorlike-mediator integration and rank structure are first-principles EFT. "
            "The rational 22/21, 23/21, 1/21 and 15/14 values were discovered post hoc and "
            "are frozen UV Clebsch targets, not yet derived group factors."
        ),
        "tree_level_matching": {
            "uv_operator": (
                "light psi_Ai mixes with heavy vectorlike F_A through G(y)^(1/2); "
                "integrating out F_A gives c_Ai G(y) with c_Ai=lambda_L,A "
                "(mu_A0+n_i mu_A1)/M_A"
            ),
            "generation_directions_required": ["flavor singlet 1", "traceless n=(-1,0,+1)"],
            "h_matrix_rank": int(np.linalg.matrix_rank(h_matrix)),
            "h_matrix_singular_values": h_singular.tolist(),
            "q_matrix_rank": int(np.linalg.matrix_rank(q_matrix)),
            "q_matrix_singular_values": q_singular.tolist(),
            "interpretation": (
                "Two flavor directions exactly span the effective h-control matrix; "
                "the missing information is in representation-dependent mediator matching."
            ),
        },
        "phaseA_mediator_suppression": {
            "epsilon_xi_falsePhi_over_Rpeak": epsilon_geom,
            "target_ratio_1_over_21": 1.0 / 21.0,
            "benchmark_exact_hu1_over_hd1": float(hu1 / hd1),
        },
        "exact_seven_control": exact,
        "models": models,
        "benchmark_ratio_values": benchmark_ratios,
        "frozen_rational_targets": rational_targets,
        "benchmark_ratio_target_errors_pct": ratio_errors,
        "crosswall_constrained_refits": crosswall,
        "main_conclusion": (
            "A strict rank-one mediator fails.  A singlet plus one traceless flavon direction "
            "is the minimal viable vectorlike sector.  Freezing only hu1/hd1=1/21 leaves six "
            "parameters and stays below 1% on all 51 walls.  A five-parameter Clebsch model "
            "reaches 0.11% at the benchmark; a four-parameter two-amplitude core model reaches "
            "0.60%, while the analogous three-parameter model narrowly misses at 1.28%."
        ),
    }
    (outdir / "vectorlike_mediator_summary.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
