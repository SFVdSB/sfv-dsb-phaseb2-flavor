#!/usr/bin/env python3
"""One-loop seed-sector beta-function audit for the SFV/dSB Pati-Salam construction.

Model:
  X0: real Pati-Salam singlet seed
  X4: real adjoint of SU(4)_C (15 components)
  XL: real adjoint of SU(2)_L (3 components)
  XR: real adjoint of SU(2)_R (3 components)
  S2: real gauge-singlet mediator coordinate with cubic source -mu_a S2 X_a^2

At the O(22) boundary all seed masses M_a^2 and S2 couplings mu_a are equal.
The exact one-loop gauge-breaking pieces are
  (16 pi^2) beta(M_a^2)|g = -6 C2(a) g_a^2 M_a^2
  (16 pi^2) beta(mu_a)|g  = -6 C2(a) g_a^2 mu_a
for a=4,L,R, while the singlet has no gauge term. Common O(22)-symmetric
quartic/trilinear pieces cancel in the block-to-singlet ratios.

The script integrates these RGEs, derives several possible induced-kernel
scalings, and refits the four remaining flavor amplitudes for each N_eff.
"""
from __future__ import annotations

import json
from math import exp, log, pi, sqrt
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from explicit_ps_seed_lagrangian import fit_continuous_N_grid  # noqa: E402

LAMBDA_SEED = 2.4100e14
M_MATCH = 1.6280e13
G_UV = {"SU4": 0.580379389023194, "SU2L": 0.541853061524504, "SU2R": 0.5230584114817481}

# Minimal Pati-Salam running content used only for the numerical integration:
# 3 x [(4,2,1)+(4bar,1,2)], one complex (1,2,2) bidoublet, and the real seed adjoints.
# Extra Pati-Salam-breaking multiplets would change b_a but not the analytic seed beta structure.
BLOCKS = [
    {"block": "SU4", "multiplicity": 15, "C2_adj": 4.0, "b_minimal": -10.0},
    {"block": "SU2L", "multiplicity": 3, "C2_adj": 2.0, "b_minimal": -8.0 / 3.0},
    {"block": "SU2R", "multiplicity": 3, "C2_adj": 2.0, "b_minimal": -8.0 / 3.0},
]


def integrate_block(block: dict, lambda_seed: float = LAMBDA_SEED, m_match: float = M_MATCH) -> dict:
    """Integrate g, M^2/M0^2, and mu/mu0 from UV to IR."""
    g0 = G_UV[block["block"]]
    C = block["C2_adj"]
    b = block["b_minimal"]
    t0, t1 = log(lambda_seed), log(m_match)

    def rhs(_t, y):
        g, rM, rmu = y
        loop = 16.0 * pi * pi
        return [b * g**3 / loop, -6.0 * C * g * g * rM / loop, -6.0 * C * g * g * rmu / loop]

    sol = solve_ivp(rhs, (t0, t1), [g0, 1.0, 1.0], rtol=2e-11, atol=2e-13)
    g_ir, rM, rmu = sol.y[:, -1]
    return {
        **block,
        "g_UV": g0,
        "g_IR": float(g_ir),
        "M2_ratio_IR_over_UV_singlet": float(rM),
        "mu_ratio_IR_over_UV_singlet": float(rmu),
        "mu_over_M2_ratio": float(rmu / rM),
        "kernel_ratio_massloop_mu2": float(rmu**2),
        "kernel_ratio_derivative_mu2_over_M2": float(rmu**2 / rM),
        "kernel_ratio_inverse_mass_1_over_M2": float(1.0 / rM),
        "kernel_ratio_mu2_over_M4": float(rmu**2 / rM**2),
    }


def neff(rows: list[dict], key: str) -> float:
    return float(sum(r["multiplicity"] * r[key] for r in rows))


def flavor_at_neff(values: list[float]) -> pd.DataFrame:
    return fit_continuous_N_grid(np.array(values, dtype=float))


def neff_for_scale_ratio(ratio: float, scenario: str) -> float:
    rows = [integrate_block(b, lambda_seed=ratio * M_MATCH, m_match=M_MATCH) for b in BLOCKS]
    key = {
        "derivative": "kernel_ratio_derivative_mu2_over_M2",
        "inverse_mass": "kernel_ratio_inverse_mass_1_over_M2",
        "protected": "kernel_ratio_mu2_over_M4",
    }[scenario]
    return neff(rows, key)


def main() -> None:
    out = ROOT / "results" / "seed_sector_beta"
    out.mkdir(parents=True, exist_ok=True)

    rows = [integrate_block(b) for b in BLOCKS]
    blocks_df = pd.DataFrame(rows)
    blocks_df.to_csv(out / "one_loop_seed_running_by_block.csv", index=False)

    scenarios = [
        {
            "scenario": "O22 boundary imposed at matching scale",
            "kernel_scaling": "k_a/k_0=1",
            "status": "protected boundary",
            "N_eff": 21.0,
            "interpretation": "No Pati-Salam running interval between equality boundary and mediator matching.",
        },
        {
            "scenario": "source-to-mass-squared protected ratio",
            "kernel_scaling": "mu_a/M_a^2",
            "status": "one-loop gauge protected",
            "N_eff": 21.0,
            "interpretation": "Because beta_ln(mu)=beta_ln(M^2), the ratio mu/M^2 is unchanged by gauge running at one loop.",
        },
        {
            "scenario": "composite derivative susceptibility",
            "kernel_scaling": "mu_a^2/M_a^2",
            "status": "enhancing",
            "N_eff": neff(rows, "kernel_ratio_derivative_mu2_over_M2"),
            "interpretation": "The one-loop derivative expansion of a scalar bubble has this scaling; charged adjoint blocks are enhanced in the IR.",
        },
        {
            "scenario": "inverse-mass tree/auxiliary marker",
            "kernel_scaling": "1/M_a^2",
            "status": "suppressive",
            "N_eff": neff(rows, "kernel_ratio_inverse_mass_1_over_M2"),
            "interpretation": "A tree-level inverse-mass kernel is suppressed because gauge running raises the charged seed masses toward the IR.",
        },
        {
            "scenario": "zero-momentum loop mass marker",
            "kernel_scaling": "mu_a^2 times renormalized B0(0)",
            "status": "scheme/boundary dependent",
            "N_eff": neff(rows, "kernel_ratio_massloop_mu2"),
            "interpretation": "The local S mass term is logarithmically divergent and requires a counterterm; its finite sign and relative threshold are not fixed by the beta function alone.",
        },
        {
            "scenario": "dimension-six normalized kernel",
            "kernel_scaling": "mu_a^2/M_a^4",
            "status": "one-loop gauge protected",
            "N_eff": neff(rows, "kernel_ratio_mu2_over_M4"),
            "interpretation": "Gauge factors cancel exactly at this order because M_a^2 and mu_a run identically.",
        },
    ]

    # Flavor fits at every numerically meaningful N_eff, plus exact 21.
    nvals = sorted(set(round(float(x["N_eff"]), 9) for x in scenarios if isinstance(x["N_eff"], (int, float))))
    fdf = flavor_at_neff(nvals)
    fdf.to_csv(out / "predicted_Neff_flavor_refits.csv", index=False)
    fmap = {round(float(r.N_eff), 9): r for _, r in fdf.iterrows()}
    for s in scenarios:
        key = round(float(s["N_eff"]), 9)
        if key in fmap:
            s["max_flavor_error_pct"] = float(fmap[key].max_error_pct)
            s["rms_flavor_error_pct"] = float(fmap[key].rms_error_pct)
    pd.DataFrame(scenarios).to_csv(out / "effective_kernel_scenarios.csv", index=False)

    # The earlier continuous audit gives the exact sub-1% interval N_eff=[20.49,21.16].
    upper_ratio = brentq(lambda r: neff_for_scale_ratio(r, "derivative") - 21.16, 1.0, 5.0)
    lower_ratio = brentq(lambda r: neff_for_scale_ratio(r, "inverse_mass") - 20.49, 1.0, 5.0)

    summary = {
        "version": "1.0.0",
        "model": {
            "fields": "X0 singlet plus X4,XL,XR real adjoint seed blocks; S2 gauge-singlet source",
            "boundary": "M_a^2=M_0^2 and mu_a=mu_0 at Lambda_seed under O(22)_seed",
            "minimal_gauge_beta_coefficients": {b["block"]: b["b_minimal"] for b in BLOCKS},
            "claim_boundary": "The analytic gauge-breaking beta terms are model independent. Numerical g running uses a minimal Pati-Salam matter content; extra breaking multiplets change the numerical size.",
        },
        "one_loop_gauge_beta_functions": {
            "mass": "(16 pi^2) beta(M_a^2)|g = -6 C2(adj_a) g_a^2 M_a^2",
            "trilinear_source": "(16 pi^2) beta(mu_a)|g = -6 C2(adj_a) g_a^2 mu_a",
            "singlet": "no gauge contribution",
            "protected_relation": "d ln(mu_a/M_a^2)/d ln Q = 0 at one-loop gauge order",
            "common_terms": "O(22)-symmetric quartic and scalar terms cancel in adjoint/singlet ratios; non-O(22) scalar couplings would add model-dependent splitting.",
        },
        "scales": {
            "Lambda_seed_GeV": LAMBDA_SEED,
            "M_match_GeV": M_MATCH,
            "ratio": LAMBDA_SEED / M_MATCH,
            "log_interval": log(LAMBDA_SEED / M_MATCH),
        },
        "blocks": rows,
        "scenario_results": scenarios,
        "sub1pct_interval_inherited": {"N_eff_min": 20.49, "N_eff_max": 21.16},
        "maximum_scale_ratio_for_sub1pct": {
            "derivative_enhancing_mu2_over_M2": upper_ratio,
            "inverse_mass_suppressive_1_over_M2": lower_ratio,
        },
        "main_conclusion": (
            "The exact one-loop seed RGEs do not select a universal suppressive or enhancing sign by themselves; "
            "the sign belongs to the definition of the induced mediator kernel. M_a^2 and mu_a receive identical "
            "gauge factors, so mu/M^2 and mu^2/M^4 are protected. The concrete loop-induced derivative kernel "
            "mu^2/M^2 is enhancing and gives N_eff about 23.73 across the current scale interval, producing a "
            "roughly 6.5% worst flavor error. Therefore the derivative-bubble realization is ruled out unless the "
            "O(22) equality is imposed within about 19% of the mediator matching scale or additional dynamics cancels "
            "the running. A protected auxiliary/composite realization based on mu/M^2 remains viable at one loop."
        ),
    }
    (out / "seed_sector_beta_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
