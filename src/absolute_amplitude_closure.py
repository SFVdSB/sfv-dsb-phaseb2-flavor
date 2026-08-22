#!/usr/bin/env python3
"""Absolute-amplitude closure for Phase B2 v1.2.0.

This module tests whether the four remaining amplitudes of the protected
Pati-Salam mediator model can be replaced by a multiplicity-normalized
cross-sector matching law.

The frozen candidate is
    c_d0 = 0.5(g_L^-2+g_R^-2) (1+epsilon_c/21)
    zeta_F = exp(-2 epsilon_c/15)
    c_d1 = -zeta_F/(4 g_L g_R)
    a_d0 = zeta_F m_Phi,T / 4
    a_d1 = zeta_F L_lock / 2
where L_lock=alpha(R_mix-R_grad).  The protected 21-kernel then gives
    h_Q=(22/21) h_d0, h_u0=(23/21)h_d0, h_u1=h_d1/21,
with h_d=c_d Gmax.

The multiplicity formulas were found after inspecting the benchmark.  They
are a frozen leading-order matching hypothesis, not yet a derived loop result.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from raw_gradient_wilson_closure import context, evaluate, residual  # noqa: E402
from chiral_localization import route1_parameters  # noqa: E402

OBS = ["ct", "ut", "sb", "db", "Vus", "Vcb", "Vub"]
ALPHA = 0.6909375570964031


def candidate_amplitudes(epsilon_c: float, m_true_phi: float, lock: float,
                         g_l: float, g_r: float) -> dict[str, float]:
    weak_inverse_even = 0.5 * (g_l ** -2 + g_r ** -2)
    z0 = 1.0 + epsilon_c / 21.0
    zf = float(np.exp(-2.0 * epsilon_c / 15.0))
    return {
        "weak_inverse_even": float(weak_inverse_even),
        "Z0": float(z0),
        "Z_F": zf,
        "c_d0": float(weak_inverse_even * z0),
        "c_d1": float(-zf / (4.0 * g_l * g_r)),
        "a_d0": float(zf * m_true_phi / 4.0),
        "a_d1": float(zf * lock / 2.0),
    }


def controls_from_amplitudes(amps: dict[str, float], gmax: float) -> np.ndarray:
    hd0 = amps["c_d0"] * gmax
    hd1 = amps["c_d1"] * gmax
    return np.array([
        (22.0 / 21.0) * hd0,
        (23.0 / 21.0) * hd0,
        hd1 / 21.0,
        amps["a_d0"],
        amps["a_d1"],
        hd0,
        hd1,
    ])


def fit_four_parameter(b, H, C, env, target) -> dict:
    def map4(p):
        hd0, hd1, a0, a1 = p
        return np.array([
            (22.0 / 21.0) * hd0,
            (23.0 / 21.0) * hd0,
            hd1 / 21.0,
            a0, a1, hd0, hd1,
        ])

    fit = least_squares(
        lambda p: residual(map4(p), b, H, C, env, target),
        np.array([2.68427, -0.66266, 0.26134, 0.24186]),
        bounds=(-8.0, 8.0), max_nfev=2500,
        xtol=1e-13, ftol=1e-13, gtol=1e-13, x_scale="jac",
    )
    z = map4(fit.x)
    values, errors, _ = evaluate(z, b, H, C, env, target)
    return {
        "parameters": fit.x,
        "controls": z,
        "values": values,
        "errors": errors,
        "max_error_pct": float(np.max(np.abs(errors))),
        "rms_error_pct": float(np.sqrt(np.mean(errors * errors))),
    }


def fit_shared_zeta(c0: float, gmax: float, m_true_phi: float, lock: float,
                    g_l: float, g_r: float, b, H, C, env, target) -> dict:
    hd0 = c0 * gmax

    def controls(zeta: float) -> np.ndarray:
        c1 = -zeta / (4.0 * g_l * g_r)
        hd1 = c1 * gmax
        return np.array([
            (22.0 / 21.0) * hd0,
            (23.0 / 21.0) * hd0,
            hd1 / 21.0,
            zeta * m_true_phi / 4.0,
            zeta * lock / 2.0,
            hd0,
            hd1,
        ])

    fit = least_squares(
        lambda p: residual(controls(float(p[0])), b, H, C, env, target),
        np.array([0.9915]), bounds=(0.8, 1.2), max_nfev=1000,
        xtol=1e-14, ftol=1e-14, gtol=1e-14,
    )
    zeta = float(fit.x[0])
    z = controls(zeta)
    values, errors, _ = evaluate(z, b, H, C, env, target)
    return {
        "zeta": zeta,
        "controls": z,
        "values": values,
        "errors": errors,
        "max_error_pct": float(np.max(np.abs(errors))),
        "rms_error_pct": float(np.sqrt(np.mean(errors * errors))),
    }


def spectrum_for_controls(z: np.ndarray, b, env: np.ndarray) -> tuple[pd.DataFrame, dict]:
    loge = np.log(np.maximum(env, 1e-300))
    bgeo = -np.gradient(loge, b.x, edge_order=2)
    x = np.linspace(-24.0, 24.0, 801)
    dx = x[1] - x[0]
    pars = route1_parameters(z)
    rows = []
    for sector in ["QL", "uR", "dR"]:
        for j in range(3):
            b0 = pars[sector]["q"][j] * (
                b.O + pars[sector]["h"][j] * (
                    b.dPhi_dy**2 + b.dphi_dy**2
                ) / np.max(b.dPhi_dy**2 + b.dphi_dy**2)
            ) + bgeo
            bb = np.interp(x, b.x, b0)
            dbb = np.gradient(bb, dx, edge_order=2)
            for op, vall in [("desired", bb * bb - dbb), ("opposite", bb * bb + dbb)]:
                v = vall[1:-1]
                main = 2.0 / dx**2 + v
                off = np.full(len(v) - 1, -1.0 / dx**2)
                ham = diags([off, main, off], [-1, 0, 1], format="csr")
                ev = np.sort(eigsh(ham, k=4, which="SA", return_eigenvectors=False,
                                   tol=1e-7, maxiter=50000))
                threshold = min(bb[0] ** 2, bb[-1] ** 2)
                row = {
                    "sector": sector,
                    "generation": j + 1,
                    "operator": op,
                    "near_zero_count_abs_lt_1e3": int(np.sum(np.abs(ev) < 1e-3)),
                    "below_threshold_count": int(np.sum(ev < threshold - 1e-4)),
                    "asymptotic_threshold": float(threshold),
                }
                row.update({f"eigenvalue_{i}": float(vv) for i, vv in enumerate(ev)})
                rows.append(row)
    frame = pd.DataFrame(rows)
    desired = frame[frame.operator == "desired"]
    opposite = frame[frame.operator == "opposite"]
    summary = {
        "desired_profiles_with_exactly_one_near_zero": int(np.sum(desired.near_zero_count_abs_lt_1e3 == 1)),
        "opposite_profiles_with_zero_near_zero": int(np.sum(opposite.near_zero_count_abs_lt_1e3 == 0)),
        "minimum_opposite_eigenvalue": float(opposite.eigenvalue_0.min()),
        "minimum_nonzero_desired_eigenvalue": float(desired.eigenvalue_1.min()),
        "total_profiles": 9,
    }
    return frame, summary


def central_axis_response(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    rows = []
    for coord in ["x", "y", "z"]:
        plus_name = f"axis_{coord}_p2p0pct"
        minus_name = f"axis_{coord}_m2p0pct"
        p = df[df.name == plus_name]
        m = df[df.name == minus_name]
        if len(p) != 1 or len(m) != 1:
            continue
        p = p.iloc[0]; m = m.iloc[0]
        delta = float(p[coord] - m[coord])
        for col in columns:
            # signed log derivative for positive amplitudes, log magnitude otherwise.
            deriv = (np.log(abs(float(p[col]))) - np.log(abs(float(m[col])))) / delta
            rows.append({"coordinate": coord, "amplitude": col, "dln_amplitude_dfraction": float(deriv)})
    return pd.DataFrame(rows)


def main() -> None:
    out = ROOT / "results/absolute_amplitude_closure"
    out.mkdir(parents=True, exist_ok=True)
    b, H, G, C, gmax, env, target = context()

    micro = pd.read_csv(ROOT / "data/phaseA_amplitude_invariants_all51.csv")
    base = micro[micro.name == "baseline_zero_bias"].iloc[0]
    ps = pd.read_csv(ROOT / "results/o22_seed_protection/pati_salam_block_inputs.csv").set_index("block")
    g_l = float(ps.loc["SU2L", "g_match"])
    g_r = float(ps.loc["SU2R", "g_match"])
    lock = ALPHA * (float(base.hessian_mixing_max_radius_dimless) - float(base.R_gradient_peak_dimless))

    exact = fit_four_parameter(b, H, C, env, target)
    hd0, hd1, a0, a1 = exact["parameters"]
    exact_raw = {
        "c_d0": float(hd0 / gmax),
        "c_d1": float(hd1 / gmax),
        "a_d0": float(a0),
        "a_d1": float(a1),
    }

    amps = candidate_amplitudes(float(base.center_energy_excess_fraction),
                                float(base.m_true_Phi_dimless), lock, g_l, g_r)
    z0 = controls_from_amplitudes(amps, gmax)
    values, errors, _ = evaluate(z0, b, H, C, env, target)

    # Infer the common attenuation separately from each of the three channels.
    zetas = {
        "from_c_d1": float(exact_raw["c_d1"] / (-1.0 / (4.0 * g_l * g_r))),
        "from_a_d0": float(exact_raw["a_d0"] / (float(base.m_true_Phi_dimless) / 4.0)),
        "from_a_d1": float(exact_raw["a_d1"] / (lock / 2.0)),
    }
    zeta_values = np.array(list(zetas.values()))
    zeta_diag = {
        **zetas,
        "mean": float(np.mean(zeta_values)),
        "std": float(np.std(zeta_values, ddof=1)),
        "relative_span_pct": float(100.0 * (np.max(zeta_values) - np.min(zeta_values)) / np.mean(zeta_values)),
        "predicted_exp_minus_2eps_over_15": float(amps["Z_F"]),
    }

    shared_fit = fit_shared_zeta(amps["c_d0"], gmax, float(base.m_true_Phi_dimless),
                                 lock, g_l, g_r, b, H, C, env, target)

    comparison_rows = []
    for name in ["c_d0", "c_d1", "a_d0", "a_d1"]:
        comparison_rows.append({
            "amplitude": name,
            "four_parameter_refit": exact_raw[name],
            "zero_fit_prediction": amps[name],
            "relative_error_pct": 100.0 * (amps[name] / exact_raw[name] - 1.0),
        })
    pd.DataFrame(comparison_rows).to_csv(out / "baseline_amplitude_comparison.csv", index=False)
    pd.DataFrame([
        {"source": k, "zeta": v} for k, v in zetas.items()
    ] + [{"source": "mean", "zeta": zeta_diag["mean"]},
         {"source": "multiplicity_prediction", "zeta": amps["Z_F"]}]
    ).to_csv(out / "shared_attenuation_inference.csv", index=False)

    pd.DataFrame({
        "observable": OBS,
        "target": target,
        "zero_fit_value": values,
        "zero_fit_error_pct": errors,
        "four_parameter_value": exact["values"],
        "four_parameter_error_pct": exact["errors"],
    }).to_csv(out / "baseline_zero_fit_observables.csv", index=False)

    # Chiral spectrum.
    spec, spec_summary = spectrum_for_controls(z0, b, env)
    spec.to_csv(out / "zero_fit_partner_spectrum.csv", index=False)
    (out / "zero_fit_partner_spectrum_summary.json").write_text(json.dumps(spec_summary, indent=2) + "\n")

    # Cross-wall tracking against independently refitted four-amplitude values.
    fits = pd.read_csv(ROOT / "results/mediator_closure/four_param_core_clebsch_all51.csv")
    gmax_df = pd.read_csv(ROOT / "results/raw_gradient_wilson/raw_gradient_wilson_coefficients_all51.csv")[["name", "Gmax_y"]]
    cw = fits.merge(gmax_df, on="name").merge(micro, on="name", suffixes=("", "_micro"))
    cw["lock"] = ALPHA * (cw.hessian_mixing_max_radius_dimless - cw.R_gradient_peak_dimless)
    cw["c_d0_fit"] = cw.hd0 / cw.Gmax_y
    cw["c_d1_fit"] = cw.hd1 / cw.Gmax_y
    pred_rows = []
    for row in cw.itertuples(index=False):
        aa = candidate_amplitudes(float(row.center_energy_excess_fraction),
                                  float(row.m_true_Phi_dimless), float(row.lock), g_l, g_r)
        rr = row._asdict()
        for key in ["c_d0", "c_d1", "a_d0", "a_d1"]:
            fit_key = key + "_fit" if key.startswith("c_") else key.replace("d", "").replace("_", "")
        pred_rows.append({
            "name": row.name, "design_phase": row.design_phase,
            "x": row.x, "y": row.y, "z": row.z,
            "c_d0_fit": row.c_d0_fit, "c_d0_pred": aa["c_d0"],
            "c_d1_fit": row.c_d1_fit, "c_d1_pred": aa["c_d1"],
            "a_d0_fit": row.a0, "a_d0_pred": aa["a_d0"],
            "a_d1_fit": row.a1, "a_d1_pred": aa["a_d1"],
            "Z_F_pred": aa["Z_F"],
        })
    track = pd.DataFrame(pred_rows)
    for key in ["c_d0", "c_d1", "a_d0", "a_d1"]:
        track[key + "_relative_error_pct"] = 100.0 * (track[key + "_pred"] / track[key + "_fit"] - 1.0)
    track.to_csv(out / "crosswall_amplitude_tracking.csv", index=False)

    tracking_summary = {}
    for label, frame in [("local33", track[track.design_phase != "corridor"]), ("all51", track)]:
        tracking_summary[label] = {}
        for key in ["c_d0", "c_d1", "a_d0", "a_d1"]:
            ee = frame[key + "_relative_error_pct"].to_numpy(float)
            tracking_summary[label][key] = {
                "mean_abs_error_pct": float(np.mean(np.abs(ee))),
                "maximum_abs_error_pct": float(np.max(np.abs(ee))),
                "mean_signed_error_pct": float(np.mean(ee)),
            }
    (out / "crosswall_tracking_summary.json").write_text(json.dumps(tracking_summary, indent=2) + "\n")

    # Compare local central response of refitted and predicted amplitudes.
    response_input = track.copy()
    response_input["c_d0_fit_abs"] = response_input.c_d0_fit
    response_input["c_d1_fit_abs"] = response_input.c_d1_fit
    response_input["a_d0_fit_abs"] = response_input.a_d0_fit
    response_input["a_d1_fit_abs"] = response_input.a_d1_fit
    response_input["c_d0_pred_abs"] = response_input.c_d0_pred
    response_input["c_d1_pred_abs"] = response_input.c_d1_pred
    response_input["a_d0_pred_abs"] = response_input.a_d0_pred
    response_input["a_d1_pred_abs"] = response_input.a_d1_pred
    resp = central_axis_response(response_input, [
        "c_d0_fit_abs", "c_d1_fit_abs", "a_d0_fit_abs", "a_d1_fit_abs",
        "c_d0_pred_abs", "c_d1_pred_abs", "a_d0_pred_abs", "a_d1_pred_abs",
    ])
    resp.to_csv(out / "amplitude_response_matrix.csv", index=False)

    # Verification tests.
    tests = {
        "zero_fit_max_error_below_1pct": bool(np.max(np.abs(errors)) < 1.0),
        "shared_zeta_relative_span_below_0p05pct": bool(zeta_diag["relative_span_pct"] < 0.05),
        "candidate_zeta_within_0p05pct_of_inferred_mean": bool(abs(amps["Z_F"] / zeta_diag["mean"] - 1.0) * 100.0 < 0.05),
        "all_four_amplitudes_within_0p1pct_at_baseline": bool(max(abs(r["relative_error_pct"]) for r in comparison_rows) < 0.1),
        "desired_9_of_9": spec_summary["desired_profiles_with_exactly_one_near_zero"] == 9,
        "opposite_9_of_9": spec_summary["opposite_profiles_with_zero_near_zero"] == 9,
        "protected_representation_ratios_exact": bool(np.allclose([
            z0[0] / z0[5], z0[1] / z0[5], z0[2] / z0[6]
        ], [22/21, 23/21, 1/21], rtol=0, atol=1e-13)),
    }
    (out / "verification_tests.json").write_text(json.dumps(tests, indent=2) + "\n")

    summary = {
        "version": "1.2.0",
        "claim_boundary": (
            "The protected 21 representation ratios are inherited from the explicit auxiliary derivation. "
            "The four absolute-amplitude formulas are a frozen multiplicity-normalized leading-order matching "
            "hypothesis identified after inspecting the benchmark. Their numerical success is not yet an independent "
            "loop derivation or blind prediction."
        ),
        "inputs": {
            "g_L": g_l, "g_R": g_r,
            "center_energy_excess_fraction": float(base.center_energy_excess_fraction),
            "m_true_Phi_dimless": float(base.m_true_Phi_dimless),
            "L_lock": float(lock), "Gmax_y": float(gmax),
        },
        "candidate_formulas": {
            "c_d0": "0.5*(g_L^-2+g_R^-2)*(1+epsilon_c/21)",
            "Z_F": "exp(-2*epsilon_c/15)",
            "c_d1": "-Z_F/(4*g_L*g_R)",
            "a_d0": "Z_F*m_true_Phi/4",
            "a_d1": "Z_F*L_lock/2",
            "representation_relations": "h_Q=(22/21)h_d0, h_u0=(23/21)h_d0, h_u1=h_d1/21",
        },
        "candidate_amplitudes": amps,
        "four_parameter_refit_amplitudes": exact_raw,
        "baseline_amplitude_relative_errors_pct": {
            row["amplitude"]: row["relative_error_pct"] for row in comparison_rows
        },
        "shared_attenuation": zeta_diag,
        "zero_fit_flavor": {
            "controls": z0.tolist(),
            "values": dict(zip(OBS, values.tolist())),
            "errors_pct": dict(zip(OBS, errors.tolist())),
            "max_error_pct": float(np.max(np.abs(errors))),
            "rms_error_pct": float(np.sqrt(np.mean(errors * errors))),
        },
        "one_fitted_common_attenuation_diagnostic": {
            "zeta": shared_fit["zeta"],
            "max_error_pct": shared_fit["max_error_pct"],
            "rms_error_pct": shared_fit["rms_error_pct"],
        },
        "four_parameter_reference": {
            "max_error_pct": exact["max_error_pct"],
            "rms_error_pct": exact["rms_error_pct"],
        },
        "chiral_spectrum": spec_summary,
        "crosswall_tracking": tracking_summary,
        "verification_tests": tests,
        "main_conclusion": (
            "At the benchmark, the four free amplitudes collapse to two physically organized normalizations and, "
            "under the frozen multiplicity attenuation Z_F=exp(-2 epsilon_c/15), to a zero-continuous-fit flavor "
            "realization with sub-percent accuracy. The formulas track the local 33-wall design only at the 0.5-3% "
            "amplitude level, so the baseline closure is promising but not yet a universal first-principles law."
        ),
    }
    (out / "absolute_amplitude_closure_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
