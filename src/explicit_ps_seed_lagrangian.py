#!/usr/bin/env python3
"""Explicit Pati-Salam seed Lagrangian and N=21 threshold audit.

The construction uses two real scalar mediator coordinates S=(S1,S2) and a
heavy seed sector 1 + Adj(G_PS).  The singlet seed couples to S1+S2; every
adjoint seed component couples to S2.  At quadratic order, any common seed
threshold produces the Gram matrix

    C^T C = [[1,1],[1,1+dim Adj(G_PS)]] = [[1,1],[1,22]].

This script derives the matrix, the T3R contractions, the continuous N_eff
flavor tolerance, and the amount of singlet/adjoint threshold equality needed.
"""
from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from raw_gradient_wilson_closure import context, evaluate, residual  # noqa: E402

OBS = ["ct", "ut", "sb", "db", "Vus", "Vcb", "Vub"]


def ps_group_data() -> dict:
    dim_su4 = 4**2 - 1
    dim_su2l = 2**2 - 1
    dim_su2r = 2**2 - 1
    dim_adj = dim_su4 + dim_su2l + dim_su2r
    c2_4 = Fraction(15, 8)
    c2_2 = Fraction(3, 4)
    return {
        "group": "SU(4)_C x SU(2)_L x SU(2)_R",
        "adjoint_dimensions": {"SU4": dim_su4, "SU2L": dim_su2l, "SU2R": dim_su2r},
        "dim_adjoint_total": dim_adj,
        "seed_space_dimension_1_plus_adj": 1 + dim_adj,
        "component_level_C2_for_(4,2)": str(c2_4 + c2_2),
        "component_level_verdict": (
            "An elementary (4,2) component receives C2=21/8 from gauge dressing, not 21. "
            "The exact integer 21 therefore requires a multiplicity/trace threshold, not a bare component Casimir."
        ),
    }


def incidence_and_gram(dim_adj: int = 21) -> tuple[np.ndarray, np.ndarray]:
    # row 0: singlet seed couples to S1+S2; rows 1..dim_adj: adjoint components couple to S2.
    C = np.zeros((1 + dim_adj, 2), float)
    C[0] = [1.0, 1.0]
    C[1:, 1] = 1.0
    return C, C.T @ C


def exact_contractions(N: float) -> dict:
    M = np.array([[1.0, 1.0], [1.0, 1.0 + N]], float)
    Minv = np.linalg.inv(M)
    e1 = np.array([1.0, 0.0])
    e2 = np.array([0.0, 1.0])
    sectors = {"Q_L": 0.0, "u_R": 0.5, "d_R": -0.5}
    factors = {}
    for name, t3r in sectors.items():
        r = np.array([1.0, -2.0 * t3r])
        factors[name] = float(e1 @ Minv @ r)
    return {
        "N_eff": float(N),
        "matrix": M.tolist(),
        "inverse": Minv.tolist(),
        "sector_source_vector": "r_A=(1,-2 T3R(A))",
        "sector_factors": factors,
        "ratios_to_down": {k: float(v / factors["d_R"]) for k, v in factors.items()},
        "odd_up_path_e2_Minv_e2": float(e2 @ Minv @ e2),
        "odd_down_reference": float(factors["d_R"]),
        "odd_up_over_down": float((e2 @ Minv @ e2) / factors["d_R"]),
        "determinant": float(np.linalg.det(M)),
        "eigenvalues": np.linalg.eigvalsh(M).tolist(),
    }


def fit_continuous_N_grid(N_values: np.ndarray) -> pd.DataFrame:
    b, H, G, C, Gmax, env, target = context()

    def mapz(p: np.ndarray, N: float) -> np.ndarray:
        hd0, hd1, a0, a1 = p
        return np.array(
            [
                (N + 1.0) / N * hd0,
                (N + 2.0) / N * hd0,
                hd1 / N,
                a0,
                a1,
                hd0,
                hd1,
            ]
        )

    # Start from the known N=21 optimum.
    seed = np.array([2.684270573214709, -0.6626623294011261, 0.2613395807102711, 0.241864927584543])
    order = sorted(range(len(N_values)), key=lambda i: abs(float(N_values[i]) - 21.0))
    solved: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    p_seed = seed.copy()
    for idx in order:
        N = float(N_values[idx])
        fit = least_squares(
            lambda p: residual(mapz(p, N), b, H, C, env, target),
            p_seed,
            bounds=(-8.0, 8.0),
            max_nfev=1800,
            xtol=2e-11,
            ftol=2e-11,
            gtol=2e-11,
            x_scale="jac",
        )
        vals, errs, _ = evaluate(mapz(fit.x, N), b, H, C, env, target)
        solved[idx] = (fit.x.copy(), vals, errs)
        p_seed = fit.x.copy()

    rows = []
    for idx, N in enumerate(N_values):
        p, vals, errs = solved[idx]
        row = {
            "N_eff": float(N),
            "threshold_ratio_kAdj_over_k0": float(N / 21.0),
            "threshold_mismatch_pct": float(100.0 * (N / 21.0 - 1.0)),
            "max_error_pct": float(np.max(np.abs(errs))),
            "rms_error_pct": float(np.sqrt(np.mean(errs * errs))),
            "h_d0": float(p[0]),
            "h_d1": float(p[1]),
            "a_d0": float(p[2]),
            "a_d1": float(p[3]),
        }
        row.update({f"error_{k}_pct": float(v) for k, v in zip(OBS, errs)})
        rows.append(row)
    return pd.DataFrame(rows)


def interval_from_grid(df: pd.DataFrame, threshold: float = 1.0) -> dict:
    good = df[df.max_error_pct < threshold].copy()
    if good.empty:
        return {"threshold_pct": threshold, "exists": False}
    best = df.loc[df.max_error_pct.idxmin()]
    return {
        "threshold_pct": threshold,
        "exists": True,
        "N_eff_min": float(good.N_eff.min()),
        "N_eff_max": float(good.N_eff.max()),
        "kAdj_over_k0_min": float(good.threshold_ratio_kAdj_over_k0.min()),
        "kAdj_over_k0_max": float(good.threshold_ratio_kAdj_over_k0.max()),
        "mismatch_pct_min": float(good.threshold_mismatch_pct.min()),
        "mismatch_pct_max": float(good.threshold_mismatch_pct.max()),
        "best_continuous_N_eff": float(best.N_eff),
        "best_continuous_max_error_pct": float(best.max_error_pct),
        "grid_step": float(np.diff(np.sort(df.N_eff.unique())).min()),
    }



def radiative_splitting_estimate() -> pd.DataFrame:
    """Order-marker for O(22)-breaking from Pati-Salam gauge running.

    Uses delta k/k ~ c_g g^2 Cbar/(16 pi^2) ln(Lambda/M), with c_g=1.
    This is not a full beta-function calculation; it only gauges plausibility.
    """
    cbar=(15*4+3*2+3*2)/21
    rows=[]
    for g in (0.45,0.55,0.65,0.75):
        for hierarchy in (2,3,10,30,100):
            delta=(g*g*cbar/(16*np.pi*np.pi))*np.log(hierarchy)
            rows.append({
                'g_PS':g,
                'Lambda_over_M':hierarchy,
                'weighted_adjoint_C2':cbar,
                'estimated_abs_threshold_splitting_pct_cg1':100*delta,
                'within_broad_1pct_flavor_tolerance_if_favorable_sign':bool(delta<=0.024285714285714355),
            })
    return pd.DataFrame(rows)

def seed_threshold_formula() -> dict:
    return {
        "seed_fields": {
            "X0": "one heavy Pati-Salam singlet seed",
            "XA": "one complete real adjoint seed X^A, A=1,...,21",
            "S": "two heavy real scalar mediator coordinates S1,S2",
        },
        "renormalizable_seed_lagrangian": (
            "L_seed = 1/2 (D X_A)^2 + 1/2 (d X_0)^2 "
            "-1/2 [M_X^2+2 mu S2] X_A X_A "
            "-1/2 [M_X^2+2 mu (S1+S2)] X_0^2."
        ),
        "quadratic_threshold": (
            "Delta V^(2)=1/2 kappa_X S^T [v0 v0^T + dim(G_PS) vA vA^T] S, "
            "with v0=(1,1), vA=(0,1)."
        ),
        "general_nonuniversal_threshold": (
            "M_S^2/k0 = [[1,1],[1,1+N_eff]], N_eff=21 kAdj/k0."
        ),
        "equal_threshold_limit": "kAdj=k0 => N_eff=dim(G_PS)=21 and M_S^2 proportional to [[1,1],[1,22]].",
        "technical_naturalness": (
            "Equality of the common S2 coupling and seed mass can be imposed by a seed-space O(22)-symmetric "
            "boundary condition before Pati-Salam gauging.  Pati-Salam gauge interactions break that equality "
            "radiatively, so the continuous tolerance audit quantifies the matching accuracy required."
        ),
    }


def effective_lagrangian() -> dict:
    return {
        "wall_source": "J_W = mu_W I_G(y) e1, e1=(1,0)",
        "sector_source": "J_A = y_A (bar Psi_A Psi_A) r_A, r_A=(1,-2 T3R(A))",
        "integrated_operator": (
            "Delta L_eff = mu_W y_A I_G(y) (bar Psi_A Psi_A) e1^T (M_S^2)^(-1) r_A"
        ),
        "exact_sector_pattern_at_N21": {
            "d_R": "1",
            "Q_L": "22/21",
            "u_R": "23/21",
        },
        "family_direction": "A real family spurion F=diag(-1,0,+1) supplies n_i.",
        "odd_path": (
            "Use the same two mediator coordinates with P_u=(1+2T3R)/2 and P_d=(1-2T3R)/2: "
            "the down odd path uses e1^T M^-1 (1,1)=1, while the up odd path uses e2^T M^-1 e2=1/21."
        ),
        "claim_boundary": (
            "This derives the rational ratios from an explicit local threshold Lagrangian.  It does not yet derive "
            "the four absolute amplitudes or prove the O(22)-symmetric seed boundary condition from the SFV action."
        ),
    }


def main() -> None:
    out = ROOT / "results/explicit_ps_seed_lagrangian"
    out.mkdir(parents=True, exist_ok=True)

    group = ps_group_data()
    C, gram = incidence_and_gram(group["dim_adjoint_total"])
    pd.DataFrame(C, columns=["coupling_to_S1", "coupling_to_S2"]).to_csv(out / "seed_incidence_matrix_22x2.csv", index_label="seed_index")
    pd.DataFrame(gram, index=["S1", "S2"], columns=["S1", "S2"]).to_csv(out / "seed_gram_mass_matrix.csv")

    # Fine continuous audit around the selected integer.
    N_values = np.round(np.linspace(19.5, 22.5, 301), 10)
    scan = fit_continuous_N_grid(N_values)
    scan.to_csv(out / "continuous_Neff_flavor_tolerance.csv", index=False)
    interval_1 = interval_from_grid(scan, 1.0)
    interval_06 = interval_from_grid(scan, 0.6)
    split = radiative_splitting_estimate()
    split.to_csv(out / "radiative_seed_splitting_order_marker.csv", index=False)

    result = {
        "version": "0.8.0",
        "group_data": group,
        "seed_threshold_lagrangian": seed_threshold_formula(),
        "incidence_matrix_shape": list(C.shape),
        "gram_matrix": gram.tolist(),
        "gram_determinant": float(np.linalg.det(gram)),
        "exact_N21_contractions": exact_contractions(21.0),
        "effective_wall_fermion_lagrangian": effective_lagrangian(),
        "continuous_flavor_tolerance_below_1pct": interval_1,
        "continuous_flavor_tolerance_below_0p6pct": interval_06,
        "radiative_splitting_order_marker": {
            "weighted_adjoint_C2": float((15*4+3*2+3*2)/21),
            "formula": "|delta k/k| ~ c_g g_PS^2 Cbar/(16 pi^2) ln(Lambda/M), shown with c_g=1",
            "claim_boundary": "order-of-magnitude marker only; a complete seed-sector RGE is not yet derived"
        },
        "best_grid_row": scan.loc[scan.max_error_pct.idxmin()].to_dict(),
        "main_conclusion": (
            "A fully local renormalizable seed threshold with one Pati-Salam singlet plus a complete 21-component "
            "adjoint generates the mediator Gram matrix [[1,1],[1,22]] at quadratic order.  T3R then generates "
            "the 21:22:23 sector offsets, and the same inverse matrix generates the 1:21 odd-slope ratio.  The "
            "remaining nontrivial assumption is equality of the singlet and per-adjoint-component threshold kernels; "
            "its required accuracy is quantified by the continuous flavor scan."
        ),
    }
    (out / "explicit_ps_seed_lagrangian_summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
