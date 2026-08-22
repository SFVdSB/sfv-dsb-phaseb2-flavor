#!/usr/bin/env python3
"""Extended 51-wall audit of candidate control relations.

Requires the full Part-I profile repository and the Phase-A dictionary.
This script re-fits all seven Route-I controls at each wall and compares them
with candidate Phase-A formulas. It is not run by the compact `run_all.sh`.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from chiral_localization import (  # noqa: E402
    CoordinateMap, build_wall_basis, derive_baseline_normalization,
    evaluate_route1, load_profile,
)

OBS = ["ct", "ut", "sb", "db", "Vus", "Vcb", "Vub"]
CONTROLS = ["h_Q", "h_u0", "h_u1", "a_d0", "a_d1", "h_d0", "h_d1"]
CANONICAL = np.array([
    1.9243768271959274, 0.37591840296783147, -0.033295237698478884,
    0.31709198515056736, 0.48708450408784837, 0.10405120245019506,
    -0.14850198624291774,
])


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--part1-root", type=Path, required=True)
    p.add_argument("--phaseA-root", type=Path, required=True)
    p.add_argument("--output", type=Path,
                   default=ROOT / "results/refit_controls_all_walls.csv")
    return p.parse_args()


def main() -> None:
    a = parse_args()
    targets_json = json.loads((ROOT / "configs/targets_MZ.json").read_text())
    target = np.array([targets_json["targets"][k] for k in OBS])
    cmap = CoordinateMap()
    base_profile = load_profile(a.part1_root / "results/profiles/baseline_zero_bias.csv")
    norm = derive_baseline_normalization(base_profile, cmap)
    dictionary = pd.read_csv(
        a.phaseA_root / "results/microphysics_dictionary_all_points.csv"
    ).set_index("name")

    def values(basis, z):
        o = evaluate_route1(basis, z)["observables"]["values"]
        return np.array([o[k] for k in OBS])

    def residual(v):
        return np.r_[np.log(v[:4] / target[:4]), (v[4:] - target[4:]) / target[4:]]

    rows = []
    warm = CANONICAL.copy()
    for name, d in dictionary.iterrows():
        if d.design_phase == "corridor":
            path = a.part1_root / "results/corridor/profiles" / f"{name}.csv"
        else:
            path = a.part1_root / "results/profiles" / f"{name}.csv"
        basis = build_wall_basis(load_profile(path), cmap, 24.0, 0.01, norm)
        fit = least_squares(
            lambda z: residual(values(basis, z)), warm, max_nfev=1000,
            xtol=2e-11, ftol=2e-11, gtol=2e-11, x_scale="jac",
        )
        z = fit.x
        v = values(basis, z)
        errors = 100 * (v / target - 1)
        predicted = np.array([
            1.5 / d.Phi_gradient_integral,
            np.sqrt(np.pi) / (2 * d.rho),
            -1 / (3 * d.tachyonic_soft_end_dimless),
            d.xi_true_Phi_dimless / 3,
            cmap.alpha * (d.hessian_mixing_max_radius_dimless - d.R_gradient_peak_dimless),
            (5 / 3) * d.center_energy_excess_fraction,
            -0.5 * d.phi_gradient_fraction,
        ])
        row = {
            "name": name, "design_phase": d.design_phase,
            "x": d.x, "y": d.y, "z": d.z,
            "success": fit.success, "nfev": fit.nfev,
            "max_pct": np.max(np.abs(errors)),
        }
        row.update({f"fit_{n}": z[i] for i, n in enumerate(CONTROLS)})
        row.update({f"pred_{n}": predicted[i] for i, n in enumerate(CONTROLS)})
        row.update({f"err_{n}": errors[i] for i, n in enumerate(OBS)})
        rows.append(row)
        warm = z
    a.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(a.output, index=False)
    print(f"Wrote {len(rows)} rows to {a.output}")


if __name__ == "__main__":
    main()
