#!/usr/bin/env python3
"""Residual-spurion derivation audit for Z0 and ZF (Phase B2 v1.3.0).

The explicit protected seed action fixes the 21 representation kernel, but the
absolute-amplitude checkpoint introduced
    Z0 = 1 + epsilon_c/21
    ZF = exp(-2 epsilon_c/15)
as a post-hoc benchmark closure.

This module asks what a local residual-relaxation spurion R can actually derive.
The minimal linked-response form is
    Z0 = exp(+lambda0 epsilon_c/21)
    ZF = exp(-2 lambdaF epsilon_c/15).
A common O(22)-normalized residual coupling sets lambda0=lambdaF=lambda_res.
Multiplicity fixes the projector traces 21 and 15, but does not fix the
remaining Wilson coefficient lambda_res.  The script quantifies the benchmark,
one- and two-coefficient fits, alternative resummations, operator counting, and
the inappropriate use of 51-wall compensation refits as a strict universality
condition.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from absolute_amplitude_closure import context, controls_from_amplitudes  # noqa: E402
from raw_gradient_wilson_closure import evaluate, residual  # noqa: E402

OBS = ["ct", "ut", "sb", "db", "Vus", "Vcb", "Vub"]


def amplitudes(eps: float, mphi: float, lock: float, gl: float, gr: float,
               lambda0: float = 1.0, lambdaf: float = 1.0,
               scheme: str = "exponential") -> dict[str, float]:
    weak = 0.5 * (gl ** -2 + gr ** -2)
    if scheme == "exponential":
        z0 = float(np.exp(lambda0 * eps / 21.0))
        zf = float(np.exp(-2.0 * lambdaf * eps / 15.0))
    elif scheme == "linear_exp":
        z0 = float(1.0 + lambda0 * eps / 21.0)
        zf = float(np.exp(-2.0 * lambdaf * eps / 15.0))
    elif scheme == "propagator":
        z0 = float(1.0 / (1.0 - lambda0 * eps / 21.0))
        zf = float((1.0 + lambdaf * eps / 15.0) ** -2.0)
    elif scheme == "determinant":
        z0 = float((1.0 + lambda0 * eps) ** (1.0 / 21.0))
        zf = float((1.0 + lambdaf * eps) ** (-2.0 / 15.0))
    else:
        raise ValueError(f"unknown scheme {scheme}")
    return {
        "Z0": z0,
        "Z_F": zf,
        "weak_inverse_even": weak,
        "c_d0": weak * z0,
        "c_d1": -zf / (4.0 * gl * gr),
        "a_d0": zf * mphi / 4.0,
        "a_d1": zf * lock / 2.0,
    }


def eval_amplitudes(amp: dict[str, float], gmax: float, b, H, C, env, target) -> dict:
    z = controls_from_amplitudes(amp, gmax)
    vals, errs, _ = evaluate(z, b, H, C, env, target)
    return {
        "controls": z,
        "values": vals,
        "errors": errs,
        "max_error_pct": float(np.max(np.abs(errs))),
        "rms_error_pct": float(np.sqrt(np.mean(errs * errs))),
    }


def fit_lambdas(common: bool, eps, mphi, lock, gl, gr, gmax, b, H, C, env, target):
    if common:
        def make(p):
            return amplitudes(eps, mphi, lock, gl, gr, float(p[0]), float(p[0]))
        x0 = np.array([1.0])
    else:
        def make(p):
            return amplitudes(eps, mphi, lock, gl, gr, float(p[0]), float(p[1]))
        x0 = np.array([1.0, 1.0])
    fit = least_squares(
        lambda p: residual(controls_from_amplitudes(make(p), gmax), b, H, C, env, target),
        x0, bounds=(-10.0, 10.0), max_nfev=3000,
        xtol=1e-14, ftol=1e-14, gtol=1e-14,
    )
    result = eval_amplitudes(make(fit.x), gmax, b, H, C, env, target)
    result["fit_parameters"] = fit.x
    result["amplitudes"] = make(fit.x)
    return result


def infer_coefficients(summary: dict) -> dict:
    eps = float(summary["inputs"]["center_energy_excess_fraction"])
    weak = float(summary["candidate_amplitudes"]["weak_inverse_even"])
    ref = summary["four_parameter_refit_amplitudes"]
    z0 = float(ref["c_d0"] / weak)
    zf_channels = {
        "from_c_d1": float(summary["shared_attenuation"]["from_c_d1"]),
        "from_a_d0": float(summary["shared_attenuation"]["from_a_d0"]),
        "from_a_d1": float(summary["shared_attenuation"]["from_a_d1"]),
    }
    zf_mean = float(np.mean(list(zf_channels.values())))
    return {
        "Z0_inferred": z0,
        "ZF_channels": zf_channels,
        "ZF_mean": zf_mean,
        "lambda0_exponential": float(21.0 * np.log(z0) / eps),
        "lambda0_linear": float(21.0 * (z0 - 1.0) / eps),
        "lambdaF_exponential": float(-15.0 * np.log(zf_mean) / (2.0 * eps)),
    }


def crosswall_coefficient_audit(summary: dict) -> tuple[pd.DataFrame, dict]:
    track = pd.read_csv(ROOT / "results/absolute_amplitude_closure/crosswall_amplitude_tracking.csv")
    micro = pd.read_csv(ROOT / "data/phaseA_amplitude_invariants_all51.csv")
    df = track.merge(micro, on=["name", "design_phase", "x", "y", "z"])
    gl = float(summary["inputs"]["g_L"]); gr = float(summary["inputs"]["g_R"])
    weak = 0.5 * (gl ** -2 + gr ** -2)
    alpha = 0.6909375570964031
    df["lock"] = alpha * (df.hessian_mixing_max_radius_dimless - df.R_gradient_peak_dimless)
    eps = df.center_energy_excess_fraction
    df["Z0_inferred"] = df.c_d0_fit / weak
    df["ZF_from_c_d1"] = df.c_d1_fit / (-1.0 / (4.0 * gl * gr))
    df["ZF_from_a_d0"] = df.a_d0_fit / (df.m_true_Phi_dimless / 4.0)
    df["ZF_from_a_d1"] = df.a_d1_fit / (df.lock / 2.0)
    df["lambda0"] = 21.0 * np.log(df.Z0_inferred) / eps
    for key in ["c_d1", "a_d0", "a_d1"]:
        df[f"lambdaF_{key}"] = -15.0 * np.log(df[f"ZF_from_{key}"]) / (2.0 * eps)
    cols = [
        "name", "design_phase", "x", "y", "z", "center_energy_excess_fraction",
        "Z0_inferred", "ZF_from_c_d1", "ZF_from_a_d0", "ZF_from_a_d1",
        "lambda0", "lambdaF_c_d1", "lambdaF_a_d0", "lambdaF_a_d1",
    ]
    out = df[cols].copy()
    stats = {}
    for label, frame in [("local33", out[out.design_phase != "corridor"]), ("all51", out)]:
        stats[label] = {}
        for col in ["lambda0", "lambdaF_c_d1", "lambdaF_a_d0", "lambdaF_a_d1"]:
            vals = frame[col].to_numpy(float)
            stats[label][col] = {
                "mean": float(np.mean(vals)), "std": float(np.std(vals, ddof=1)),
                "min": float(np.min(vals)), "max": float(np.max(vals)),
            }
    return out, stats


def gauge_normalizer_exploration() -> pd.DataFrame:
    """Exploratory only: compare simple wall scalars with compensation-fit c_d0.

    This is not promoted to a derivation because the candidates were compared after
    seeing the compensation refits and several contain bubble-radius dependence that
    is not expected in a local 4D gauge kinetic coefficient.
    """
    full_path = Path("/mnt/data/work_phaseA/sfv-dsb-microphysical-dictionary-phaseA-v0.1.0/results/microphysics_dictionary_all_points.csv")
    if not full_path.exists():
        return pd.DataFrame()
    full = pd.read_csv(full_path)
    fits = pd.read_csv(ROOT / "results/mediator_closure/four_param_core_clebsch_all51.csv")
    gm = pd.read_csv(ROOT / "results/raw_gradient_wilson/raw_gradient_wilson_coefficients_all51.csv")[["name", "Gmax_y"]]
    df = full.merge(fits, on=["name", "design_phase", "x", "y", "z"]).merge(gm, on="name")
    df["c_d0_fit"] = df.hd0 / df.Gmax_y
    base = df[df.name == "baseline_zero_bias"].iloc[0]
    summary = json.loads((ROOT / "results/absolute_amplitude_closure/absolute_amplitude_closure_summary.json").read_text())
    weak = float(summary["candidate_amplitudes"]["weak_inverse_even"])
    candidates = {
        "R_action_equivalent": df.R_action_equivalent_dimless,
        "R_peak": df.R_peak_dimless,
        "action_FWHM": df.action_FWHM_dimless,
        "gradient_FWHM": df.gradient_FWHM_dimless,
        "sigma_gradient": df.sigma_gradient_dimless,
        "old_form_xiPhiT2_sigmaGrad_over_wGrad": (
            df.xi_true_Phi_dimless**2 * df.sigma_gradient_dimless / df.gradient_FWHM_dimless
        ),
        "old_form_xiPhiT2_sigmaAction_over_wAction": (
            df.xi_true_Phi_dimless**2 * df.sigma_action_equivalent_dimless / df.action_FWHM_dimless
        ),
    }
    rows = []
    for name, val in candidates.items():
        rel = val / float(val[df.name == "baseline_zero_bias"].iloc[0])
        pred = weak * rel * np.exp(df.center_energy_excess_fraction / 21.0)
        err = 100.0 * (pred / df.c_d0_fit - 1.0)
        loc = df.design_phase != "corridor"
        rows.append({
            "candidate": name,
            "local33_mean_abs_error_pct": float(np.mean(np.abs(err[loc]))),
            "local33_max_abs_error_pct": float(np.max(np.abs(err[loc]))),
            "all51_mean_abs_error_pct": float(np.mean(np.abs(err))),
            "all51_max_abs_error_pct": float(np.max(np.abs(err))),
            "promoted_to_derivation": False,
        })
    return pd.DataFrame(rows).sort_values("local33_mean_abs_error_pct")


def main() -> None:
    out = ROOT / "results/residual_spurion_normalization"
    out.mkdir(parents=True, exist_ok=True)
    old = json.loads((ROOT / "results/absolute_amplitude_closure/absolute_amplitude_closure_summary.json").read_text())
    b, H, G, C, gmax, env, target = context()
    eps = float(old["inputs"]["center_energy_excess_fraction"])
    mphi = float(old["inputs"]["m_true_Phi_dimless"])
    lock = float(old["inputs"]["L_lock"])
    gl = float(old["inputs"]["g_L"]); gr = float(old["inputs"]["g_R"])

    model_rows = []
    model_details = {}
    for scheme in ["linear_exp", "exponential", "propagator", "determinant"]:
        amp = amplitudes(eps, mphi, lock, gl, gr, scheme=scheme)
        res = eval_amplitudes(amp, gmax, b, H, C, env, target)
        model_rows.append({
            "scheme": scheme, "lambda0": 1.0, "lambdaF": 1.0,
            "Z0": amp["Z0"], "Z_F": amp["Z_F"],
            "max_error_pct": res["max_error_pct"], "rms_error_pct": res["rms_error_pct"],
        })
        model_details[scheme] = {
            "amplitudes": amp, "controls": res["controls"].tolist(),
            "errors_pct": dict(zip(OBS, res["errors"].tolist())),
            "max_error_pct": res["max_error_pct"], "rms_error_pct": res["rms_error_pct"],
        }
    pd.DataFrame(model_rows).to_csv(out / "normalization_resummation_comparison.csv", index=False)

    one = fit_lambdas(True, eps, mphi, lock, gl, gr, gmax, b, H, C, env, target)
    two = fit_lambdas(False, eps, mphi, lock, gl, gr, gmax, b, H, C, env, target)
    inferred = infer_coefficients(old)
    cw, cwstats = crosswall_coefficient_audit(old)
    cw.to_csv(out / "crosswall_inferred_residual_coefficients.csv", index=False)
    gauge = gauge_normalizer_exploration()
    gauge.to_csv(out / "exploratory_wall_gauge_normalizer_audit.csv", index=False)

    operator_table = pd.DataFrame([
        {"operator": "R X0^2", "symmetry": "PS singlet", "coefficient": "lambda_0", "fixed_by_multiplicity": False},
        {"operator": "R Tr(X4^2)", "symmetry": "SU(4) adjoint norm", "coefficient": "lambda_4", "fixed_by_multiplicity": False},
        {"operator": "R Tr(XL^2)", "symmetry": "SU(2)L adjoint norm", "coefficient": "lambda_L", "fixed_by_multiplicity": False},
        {"operator": "R Tr(XR^2)", "symmetry": "SU(2)R adjoint norm", "coefficient": "lambda_R", "fixed_by_multiplicity": False},
        {"operator": "R J_wall S", "symmetry": "wall-source dressing", "coefficient": "eta_0", "fixed_by_multiplicity": False},
        {"operator": "R J_F S", "symmetry": "fermion-source dressing", "coefficient": "eta_F", "fixed_by_multiplicity": False},
    ])
    operator_table.to_csv(out / "residual_spurion_operator_ledger.csv", index=False)

    exp_unit = model_details["exponential"]
    tests = {
        "unit_common_coefficient_below_1pct": exp_unit["max_error_pct"] < 1.0,
        "best_common_coefficient_near_one": abs(float(one["fit_parameters"][0]) - 1.0) < 0.1,
        "one_coefficient_model_below_1pct": one["max_error_pct"] < 1.0,
        "two_coefficient_model_below_1pct": two["max_error_pct"] < 1.0,
        "multiplicity_does_not_fix_wilson_coefficients": bool((~operator_table.fixed_by_multiplicity).all()),
        "crosswall_compensation_lambdas_not_universal": cwstats["local33"]["lambda0"]["std"] > 0.5,
    }
    (out / "verification_tests.json").write_text(json.dumps(tests, indent=2) + "\n")

    summary = {
        "version": "1.3.0",
        "derivation": {
            "minimal_residual_spurion_action": (
                "A dimensionless residual-relaxation spurion R with <R>=epsilon_c dresses the normalized "
                "PS-singlet collective source and the two fermion-source vertices. Linked local insertions give "
                "Z0=exp(lambda0 epsilon_c/21) and ZF=exp(-2 lambdaF epsilon_c/15)."
            ),
            "common_O22_normalized_limit": "lambda0=lambdaF=lambda_res",
            "multiplicity_result": (
                "21 and 15 identify normalized PS and SU(4) trace projectors, but trace normalization alone "
                "does not fix the associated Wilson coefficient."
            ),
        },
        "benchmark": {
            "epsilon_c": eps,
            "unit_common_coefficient": model_details["exponential"],
            "best_one_common_coefficient": {
                "lambda_res": float(one["fit_parameters"][0]),
                "amplitudes": one["amplitudes"],
                "errors_pct": dict(zip(OBS, one["errors"].tolist())),
                "max_error_pct": one["max_error_pct"], "rms_error_pct": one["rms_error_pct"],
            },
            "best_two_coefficients": {
                "lambda0": float(two["fit_parameters"][0]),
                "lambdaF": float(two["fit_parameters"][1]),
                "max_error_pct": two["max_error_pct"], "rms_error_pct": two["rms_error_pct"],
            },
            "coefficients_inferred_from_four_amplitude_reference": inferred,
            "resummation_comparison": model_details,
        },
        "crosswall_audit": {
            "interpretation": (
                "The 51-wall amplitudes are compensation refits chosen to hold the same observed flavor targets "
                "while changing the universe's wall. They are useful diagnostics but are not independent measurements "
                "of UV Wilson coefficients. Their non-universality cannot by itself falsify a fixed UV action."
            ),
            "inferred_coefficient_statistics": cwstats,
        },
        "gauge_normalization_audit": {
            "status": "not independently closed",
            "reason": (
                "The existing gauge-normalization work supplies baseline gL,gR, but its group weights/kernel "
                "coefficients were calibrated. No predeclared local wall kernel in the present repository uniquely "
                "predicts gL,gR across X,Y,Z. Exploratory correlations were not promoted to laws."
            ),
        },
        "claim_boundary": (
            "The explicit seed action reduces the four absolute amplitudes to one common order-one residual-response "
            "coefficient while preserving sub-percent flavor. The value lambda_res=1 gives a zero-continuous-fit "
            "0.632% result, and the best common value is 1.076. However O(22)/Pati-Salam multiplicity fixes the "
            "projector dimensions, not the Wilson coefficient itself. A microscopic residual-seed interaction or "
            "independent observable is still required to prove lambda_res=1 rather than assume unit normalization."
        ),
        "verification_tests": tests,
        "main_conclusion": (
            "Z0 and ZF can be generated coherently by one local residual spurion, so the four-amplitude problem "
            "collapses to one natural coefficient. The current symmetries do not uniquely determine that coefficient; "
            "therefore a full parameter-free first-principles prediction has not yet been reached."
        ),
    }
    (out / "residual_spurion_normalization_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
