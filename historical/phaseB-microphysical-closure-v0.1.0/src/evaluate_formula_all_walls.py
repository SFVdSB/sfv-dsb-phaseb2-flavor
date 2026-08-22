#!/usr/bin/env python3
"""Apply the frozen exploratory Phase-A formula map to all 51 walls.

Requires the full Part-I profile repository and the Phase-A dictionary.
No flavor controls are re-fitted.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from chiral_localization import (  # noqa: E402
    CoordinateMap, build_wall_basis, derive_baseline_normalization,
    evaluate_route1, load_profile,
)

OBS = ["ct", "ut", "sb", "db", "Vus", "Vcb", "Vub"]
CONTROLS = ["h_Q", "h_u0", "h_u1", "a_d0", "a_d1", "h_d0", "h_d1"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--part1-root", type=Path, required=True)
    p.add_argument("--phaseA-root", type=Path, required=True)
    p.add_argument("--output", type=Path,
                   default=ROOT / "results/formula_closure_all_walls.csv")
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
    )
    rows = []
    for _, d in dictionary.iterrows():
        name = d["name"]
        if d.design_phase == "corridor":
            path = a.part1_root / "results/corridor/profiles" / f"{name}.csv"
        else:
            path = a.part1_root / "results/profiles" / f"{name}.csv"
        basis = build_wall_basis(load_profile(path), cmap, 24.0, 0.01, norm)
        controls = np.array([
            1.5 / d.Phi_gradient_integral,
            np.sqrt(np.pi) / (2 * d.rho),
            -1 / (3 * d.tachyonic_soft_end_dimless),
            d.xi_true_Phi_dimless / 3,
            cmap.alpha * (d.hessian_mixing_max_radius_dimless - d.R_gradient_peak_dimless),
            (5 / 3) * d.center_energy_excess_fraction,
            -0.5 * d.phi_gradient_fraction,
        ])
        o = evaluate_route1(basis, controls)["observables"]["values"]
        values = np.array([o[k] for k in OBS])
        errors = 100 * (values / target - 1)
        row = {
            "name": name, "design_phase": d.design_phase,
            "x": d.x, "y": d.y, "z": d.z,
            "max_abs_percent_error": np.max(np.abs(errors)),
            "rms_percent_error": np.sqrt(np.mean(errors**2)),
        }
        row.update({CONTROLS[i]: controls[i] for i in range(7)})
        row.update({f"model_{OBS[i]}": values[i] for i in range(7)})
        row.update({f"error_{OBS[i]}_pct": errors[i] for i in range(7)})
        rows.append(row)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(a.output, index=False)
    print(f"Wrote {len(rows)} rows to {a.output}")


if __name__ == "__main__":
    main()
