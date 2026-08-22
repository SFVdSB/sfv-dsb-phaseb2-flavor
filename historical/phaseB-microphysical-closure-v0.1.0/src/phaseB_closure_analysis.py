#!/usr/bin/env python3
"""Phase-B microphysical closure analysis for SFV/dSB chiral localization.

This script distinguishes three levels:
1. the exact seven-control Route-I benchmark;
2. physically constrained fits in which selected controls are fixed by Phase-A wall invariants;
3. an exploratory, post-hoc zero-continuous-fit formula map.

The formula map is a hypothesis-generation device, not a first-principles derivation.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from chiral_localization import (  # noqa: E402
    CoordinateMap,
    build_wall_basis,
    derive_baseline_normalization,
    evaluate_route1,
    load_profile,
)

OBS = ["ct", "ut", "sb", "db", "Vus", "Vcb", "Vub"]
CONTROLS = ["h_Q", "h_u0", "h_u1", "a_d0", "a_d1", "h_d0", "h_d1"]
CANONICAL = np.array([
    1.9243768271959274,
    0.37591840296783147,
    -0.033295237698478884,
    0.31709198515056736,
    0.48708450408784837,
    0.10405120245019506,
    -0.14850198624291774,
])


def dump_json(path: Path, obj: object) -> None:
    def convert(v):
        if isinstance(v, np.ndarray):
            return v.tolist()
        if isinstance(v, (np.floating, np.integer)):
            return v.item()
        if isinstance(v, dict):
            return {str(k): convert(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)):
            return [convert(x) for x in v]
        return v
    path.write_text(json.dumps(convert(obj), indent=2) + "\n", encoding="utf-8")


def values(basis, controls: np.ndarray) -> np.ndarray:
    obs = evaluate_route1(basis, controls)["observables"]["values"]
    return np.array([obs[k] for k in OBS], dtype=float)


def residual(model: np.ndarray, target: np.ndarray) -> np.ndarray:
    return np.r_[np.log(model[:4] / target[:4]), (model[4:] - target[4:]) / target[4:]]


def phase_a_formula_controls(d: dict, cmap: CoordinateMap) -> tuple[np.ndarray, dict]:
    """Exploratory post-hoc formula map built only from Phase-A quantities."""
    locking_interval = cmap.alpha * (
        d["hessian_mixing_max_radius_dimless"] - d["R_gradient_peak_dimless"]
    )
    z = np.array([
        1.5 / d["Phi_gradient_integral"],
        np.sqrt(np.pi) / (2.0 * d["rho"]),
        -1.0 / (3.0 * d["tachyonic_soft_end_dimless"]),
        d["xi_true_Phi_dimless"] / 3.0,
        locking_interval,
        (5.0 / 3.0) * d["center_energy_excess_fraction"],
        -0.5 * d["phi_gradient_fraction"],
    ])
    definitions = {
        "h_Q": "3/(2 I_Phi), I_Phi=integral (dPhi/dr)^2 dr",
        "h_u0": "sqrt(pi)/(2 rho)",
        "h_u1": "-1/(3 r_tachyonic_end)",
        "a_d0": "xi_true_Phi/3",
        "a_d1": "alpha (R_max_mixing - R_total_gradient_peak)",
        "h_d0": "(5/3) times center energy-excess fraction",
        "h_d1": "-1/2 times brane-gradient fraction",
    }
    return z, definitions


def fit_with_fixed(basis, target: np.ndarray, predicted: np.ndarray,
                   fixed_indices: list[int]) -> tuple[np.ndarray, np.ndarray]:
    fixed = np.array(fixed_indices, dtype=int)
    active = np.array([i for i in range(7) if i not in fixed], dtype=int)

    def unpack(a: np.ndarray) -> np.ndarray:
        z = predicted.copy()
        z[active] = a
        return z

    if len(active):
        result = least_squares(
            lambda a: residual(values(basis, unpack(a)), target),
            CANONICAL[active],
            max_nfev=2500,
            xtol=1e-13,
            ftol=1e-13,
            gtol=1e-13,
            x_scale="jac",
        )
        z = unpack(result.x)
    else:
        z = predicted.copy()
    return z, values(basis, z)


def main() -> None:
    out = ROOT / "results"
    profile = load_profile(ROOT / "data/background_profile_O4_regular_robin_full.csv")
    cmap = CoordinateMap()
    norm = derive_baseline_normalization(profile, cmap)
    basis = build_wall_basis(profile, cmap, half_width=24.0, spacing=0.005,
                             normalization=norm)
    phase_a = json.loads((ROOT / "data/baseline_microphysics_dictionary.json").read_text())
    targets_json = json.loads((ROOT / "configs/targets_MZ.json").read_text())
    target = np.array([targets_json["targets"][k] for k in OBS], dtype=float)

    predicted, definitions = phase_a_formula_controls(phase_a, cmap)
    formula_values = values(basis, predicted)
    formula_errors = 100.0 * (formula_values / target - 1.0)

    control_map = pd.DataFrame({
        "control": CONTROLS,
        "canonical_fitted": CANONICAL,
        "phaseA_formula": predicted,
        "relative_difference": predicted / CANONICAL - 1.0,
        "formula_definition": [definitions[k] for k in CONTROLS],
        "status": [
            "post-hoc candidate", "post-hoc candidate", "post-hoc candidate",
            "post-hoc candidate", "51-wall supported candidate",
            "post-hoc candidate", "post-hoc candidate",
        ],
    })
    control_map.to_csv(out / "baseline_microphysical_control_map.csv", index=False)

    pd.DataFrame({
        "observable": OBS,
        "target": target,
        "posthoc_formula_model": formula_values,
        "percent_error": formula_errors,
    }).to_csv(out / "posthoc_formula_observables.csv", index=False)

    models = {
        "canonical_7_fit": [],
        "locking_interval_fixed": [4],
        "locking_interval_and_gradient_fraction_fixed": [4, 6],
        "four_formula_controls_fixed": [0, 1, 4, 6],
        "six_formula_controls_fixed": [0, 1, 2, 3, 4, 6],
        "all_seven_formula_controls_fixed": list(range(7)),
    }
    rows = []
    for name, fixed in models.items():
        z, model = fit_with_fixed(basis, target, predicted, fixed)
        err = 100.0 * (model / target - 1.0)
        row = {
            "model": name,
            "n_fitted": 7 - len(fixed),
            "n_fixed_by_formula": len(fixed),
            "fixed_controls": "+".join(CONTROLS[i] for i in fixed),
            "max_abs_percent_error": float(np.max(np.abs(err))),
            "rms_percent_error": float(np.sqrt(np.mean(err**2))),
        }
        row.update({f"control_{k}": float(z[i]) for i, k in enumerate(CONTROLS)})
        row.update({f"error_{k}_pct": float(err[i]) for i, k in enumerate(OBS)})
        rows.append(row)
    pd.DataFrame(rows).to_csv(out / "constrained_model_results.csv", index=False)

    refits = pd.read_csv(out / "refit_controls_all_walls.csv")
    relation_rows = []
    for control in CONTROLS:
        fitted = refits[f"fit_{control}"].to_numpy(float)
        pred = refits[f"pred_{control}"].to_numpy(float)
        rel = pred / fitted - 1.0
        relation_rows.append({
            "control": control,
            "correlation": float(np.corrcoef(fitted, pred)[0, 1]),
            "mean_abs_relative_error_pct": float(100 * np.mean(np.abs(rel))),
            "median_abs_relative_error_pct": float(100 * np.median(np.abs(rel))),
            "max_abs_relative_error_pct": float(100 * np.max(np.abs(rel))),
            "baseline_relative_error_pct": float(100 * rel[0]),
        })
    pd.DataFrame(relation_rows).to_csv(out / "all_wall_control_relation_summary.csv", index=False)

    blind = pd.read_csv(out / "formula4_blind_tests.csv")
    up = {"ct", "ut"}; down = {"sb", "db"}; mixing = {"Vus", "Vcb", "Vub"}
    balanced_mask = blind["calibration"].map(
        lambda s: len(set(s.split(",")) & up) == 1
        and len(set(s.split(",")) & down) == 1
        and len(set(s.split(",")) & mixing) == 1
    )
    balanced = blind[balanced_mask]

    spectrum = pd.read_csv(out / "posthoc_formula_partner_spectrum.csv")
    desired = spectrum[spectrum.operator == "desired"]
    opposite = spectrum[spectrum.operator == "opposite"]

    summary = {
        "version": "0.1.0",
        "status": "partial microphysical closure with an exploratory post-hoc benchmark formula",
        "accepted_empirical_relation": {
            "control": "a_d1",
            "formula": definitions["a_d1"],
            "baseline_relative_error_pct": float(
                100 * (predicted[4] / CANONICAL[4] - 1.0)
            ),
            "all_51_wall_correlation": float(
                np.corrcoef(refits.fit_a_d1, refits.pred_a_d1)[0, 1]
            ),
            "all_51_wall_mean_abs_relative_error_pct": float(
                100 * np.mean(np.abs(refits.pred_a_d1 / refits.fit_a_d1 - 1.0))
            ),
            "all_51_wall_max_abs_relative_error_pct": float(
                100 * np.max(np.abs(refits.pred_a_d1 / refits.fit_a_d1 - 1.0))
            ),
        },
        "benchmark_formula": {
            "continuous_fit_count": 0,
            "controls": dict(zip(CONTROLS, predicted)),
            "max_abs_percent_error": float(np.max(np.abs(formula_errors))),
            "rms_percent_error": float(np.sqrt(np.mean(formula_errors**2))),
            "classification": "hypothesis generated after inspecting the seven-control solution; not a prediction",
        },
        "three_fit_four_fixed_model": {
            "fixed_controls": ["h_Q", "h_u0", "a_d1", "h_d1"],
            "fitted_controls": ["h_u1", "a_d0", "h_d0"],
            "max_abs_percent_error_all_seven": float(
                pd.DataFrame(rows).set_index("model").loc[
                    "four_formula_controls_fixed", "max_abs_percent_error"
                ]
            ),
        },
        "internal_balanced_holdout_test": {
            "number_of_calibration_sets": int(len(balanced)),
            "sets_with_heldout_max_below_1pct": int(np.sum(balanced.heldout_max_pct < 1.0)),
            "median_heldout_max_percent_error": float(balanced.heldout_max_pct.median()),
            "worst_heldout_max_percent_error": float(balanced.heldout_max_pct.max()),
            "warning": "not scientifically blind because formulas were discovered using the fitted benchmark",
        },
        "spectral_audit_formula_model": {
            "desired_with_one_near_zero": int(np.sum(desired.near_zero == 1)),
            "opposite_with_zero_near_zero": int(np.sum(opposite.near_zero == 0)),
            "minimum_opposite_eigenvalue": float(opposite.eigenvalue_0.min()),
        },
        "structural_conclusion": (
            "The scalar wall determines O(x), E(x), and several geometric combinations, "
            "but it does not determine generation charges or fermion Wilson coefficients. "
            "A flavor symmetry or explicit UV fermion sector remains necessary for a first-principles Phase-B derivation."
        ),
    }
    dump_json(out / "phaseB_summary.json", summary)
    dump_json(out / "posthoc_formula_controls.json", {
        "controls": dict(zip(CONTROLS, predicted)),
        "definitions": definitions,
        "claim_boundary": "post-hoc candidate; freeze for future external tests before calling predictive",
    })
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
