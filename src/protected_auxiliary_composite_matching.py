#!/usr/bin/env python3
"""Protected auxiliary/composite matching for the SFV/dSB Pati-Salam seed construction.

This checkpoint distinguishes an algebraic/Legendre collective response from the excluded
propagating scalar-bubble derivative susceptibility.

For the gauge-invariant singlet and adjoint-norm collective coordinates Y0 and YA, use the
zero-momentum first-order response action

  Gamma_aux = 1/2 a0 Y0^2 - Y0 (S1+S2)
            + 21 [1/2 aA YA^2 - YA S2],

where a_i=M_i^2/mu_i and k_i=a_i^{-1}=mu_i/M_i^2.  Eliminating Y0,YA gives the pullback
kernel

  K_S = k0 v0 v0^T + 21 kA vA vA^T,
  v0=(1,1), vA=(0,1).

At the O(22) boundary kA=k0.  The one-loop gauge beta functions of M^2 and mu are equal,
so k=mu/M^2 remains equal under the gauge running already audited in v1.0.0.

The script verifies basis/canonical-normalization invariance, reconstructs the flavor ratios,
and quantifies finite-momentum contamination if the collective channels acquire kinetic terms.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from explicit_ps_seed_lagrangian import fit_continuous_N_grid  # noqa: E402

MULT = {"SU4": 15, "SU2L": 3, "SU2R": 3}
V0 = np.array([1.0, 1.0])
VA = np.array([0.0, 1.0])
E1 = np.array([1.0, 0.0])
E2 = np.array([0.0, 1.0])
T3R = {"Q_L": 0.0, "u_R": 0.5, "d_R": -0.5}
NMIN, NMAX = 20.49, 21.16


def kernel_from_weights(k0: float, k4: float, kL: float, kR: float) -> np.ndarray:
    return k0 * np.outer(V0, V0) + (15*k4 + 3*kL + 3*kR) * np.outer(VA, VA)


def sector_contractions(K: np.ndarray) -> dict:
    Ki = np.linalg.inv(K)
    factors = {}
    for sec, t in T3R.items():
        r = np.array([1.0, -2.0*t])
        factors[sec] = float(E1 @ Ki @ r)
    d = factors["d_R"]
    return {
        "inverse": Ki.tolist(),
        "sector_factors": factors,
        "ratios_to_down": {k: float(v/d) for k, v in factors.items()},
        "odd_up_e2_Kinv_e2": float(E2 @ Ki @ E2),
        "odd_up_over_down": float((E2 @ Ki @ E2)/d),
    }


def canonical_basis_invariance(K: np.ndarray, ntrial: int = 250, seed: int = 21022) -> dict:
    rng = np.random.default_rng(seed)
    Ki = np.linalg.inv(K)
    max_err = 0.0
    for _ in range(ntrial):
        while True:
            R = rng.normal(size=(2,2))
            if abs(np.linalg.det(R)) > 0.2:
                break
        Jw = rng.normal(size=2)
        Jf = rng.normal(size=2)
        base = float(Jw @ Ki @ Jf)
        Kp = R.T @ K @ R
        Jwp = R.T @ Jw
        Jfp = R.T @ Jf
        transformed = float(Jwp @ np.linalg.inv(Kp) @ Jfp)
        max_err = max(max_err, abs(base-transformed))
    return {"trials": ntrial, "max_abs_cross_term_error": max_err}


def load_beta_rows() -> pd.DataFrame:
    p = ROOT / "results" / "seed_sector_beta" / "one_loop_seed_running_by_block.csv"
    return pd.read_csv(p)


def protected_running_kernel() -> tuple[np.ndarray, dict]:
    df = load_beta_rows()
    ratios = {r.block: float(r.mu_over_M2_ratio) for _, r in df.iterrows()}
    K = kernel_from_weights(1.0, ratios["SU4"], ratios["SU2L"], ratios["SU2R"])
    neff = 15*ratios["SU4"] + 3*ratios["SU2L"] + 3*ratios["SU2R"]
    return K, {"block_k_over_k0": ratios, "N_eff": float(neff)}


def flavor_row_at_N(N: float) -> dict:
    df = fit_continuous_N_grid(np.array([float(N)]))
    return df.iloc[0].to_dict()


def wall_momentum_estimates() -> dict:
    phaseA = json.loads((Path("/mnt/data/work_phaseA/sfv-dsb-microphysical-dictionary-phaseA-v0.1.0/results/phaseA_summary.json")).read_text())
    v = float(phaseA["canonical_scale"]["v_brane_GeV"])
    wa = float(phaseA["wall"]["action_FWHM"])
    wg = float(phaseA["wall"]["gradient_FWHM"])
    Mmatch = 1.6280e13
    Ea = v/wa
    Eg = v/wg
    return {
        "v_brane_GeV": v,
        "action_FWHM_dimless": wa,
        "gradient_FWHM_dimless": wg,
        "inverse_action_width_GeV": Ea,
        "inverse_gradient_width_GeV": Eg,
        "M_match_GeV": Mmatch,
        "q_action_if_composite_mass_equals_Mmatch": Ea/Mmatch,
        "q_gradient_if_composite_mass_equals_Mmatch": Eg/Mmatch,
    }


def finite_momentum_neff(q: float, rho4: float, rhoL: float, rhoR: float) -> float:
    """N_eff for k_a(q)/k0(q)=(1+q^2)/(1+rho_a q^2).

    rho_a=(Z_a/M_a^2)/(Z_0/M_0^2).  rho=1 preserves exact equality at every momentum.
    """
    def rr(rho: float) -> float:
        return (1.0 + q*q)/(1.0 + rho*q*q)
    return 15*rr(rho4) + 3*rr(rhoL) + 3*rr(rhoR)


def finite_momentum_scan() -> pd.DataFrame:
    rows=[]
    qvals=[0.0,0.25,0.5,1.0,2.0,3.0,3.6]
    deltas=[-0.03,-0.02,-0.01,0.0,0.01,0.02,0.03,0.05]
    # use a common relative derivative-normalization mismatch for a transparent envelope
    for q in qvals:
        for d in deltas:
            rho=1.0+d
            N=finite_momentum_neff(q,rho,rho,rho)
            rows.append({
                "q_p_over_Mcomp":q,
                "rho_ZoverM2_adj_over_singlet":rho,
                "derivative_normalization_mismatch_pct":100*d,
                "N_eff":N,
                "within_inherited_sub1pct_N_interval":bool(NMIN <= N <= NMAX),
            })
    return pd.DataFrame(rows)


def mismatch_bounds(q: float) -> dict:
    """Common rho-1 bounds implied by N interval."""
    if q == 0:
        return {"q":0.0,"delta_min":None,"delta_max":None,"note":"no momentum sensitivity"}
    # r=N/21=(1+q2)/(1+(1+delta)q2)
    def delta_for_N(N: float) -> float:
        r=N/21.0
        return ((1+q*q)/r - 1.0)/(q*q) - 1.0
    # N increases when delta decreases
    return {
        "q":q,
        "delta_min":delta_for_N(NMAX),
        "delta_max":delta_for_N(NMIN),
        "delta_min_pct":100*delta_for_N(NMAX),
        "delta_max_pct":100*delta_for_N(NMIN),
    }


def main() -> None:
    out = ROOT / "results" / "protected_auxiliary_matching"
    out.mkdir(parents=True, exist_ok=True)

    K_exact = kernel_from_weights(1.0,1.0,1.0,1.0)
    exact = sector_contractions(K_exact)
    basis = canonical_basis_invariance(K_exact)
    K_run, run = protected_running_kernel()
    run_contractions = sector_contractions(K_run)
    flavor = flavor_row_at_N(run["N_eff"])
    momenta = wall_momentum_estimates()
    fm = finite_momentum_scan()
    fm.to_csv(out / "finite_momentum_composite_contamination_scan.csv", index=False)
    bounds = pd.DataFrame([mismatch_bounds(q) for q in [0.25,0.5,1.0,2.0,3.0,3.6]])
    bounds.to_csv(out / "finite_momentum_allowed_derivative_mismatch.csv", index=False)

    # Auxiliary incidence ledger: one singlet collective plus one normalized adjoint norm carrying multiplicity 21.
    ledger = pd.DataFrame([
        {"collective":"Y0","gauge_object":"X0^2","multiplicity":1,"incidence_S1":1,"incidence_S2":1,"response":"k0=mu0/M0^2"},
        {"collective":"YAdj","gauge_object":"(1/21) Tr_Adj X^2","multiplicity":21,"incidence_S1":0,"incidence_S2":1,"response":"kAdj=muAdj/MAdj^2"},
    ])
    ledger.to_csv(out / "gauge_invariant_collective_channel_ledger.csv", index=False)
    pd.DataFrame(K_exact,index=["S1","S2"],columns=["S1","S2"]).to_csv(out / "protected_auxiliary_gram_matrix.csv")

    summary={
        "version":"1.1.0",
        "matching_type":"zero-momentum auxiliary/Legendre collective response",
        "gauge_invariant_collectives":{
            "singlet":"Y0 conjugate to X0^2",
            "adjoint_norm":"YAdj conjugate to (1/21) Tr_Adj X^2; multiplicity 21 is explicit",
        },
        "first_order_response_action":(
            "Gamma_aux = 1/2 (M0^2/mu0) Y0^2 - Y0(S1+S2) "
            "+ 21[1/2 (MAdj^2/muAdj) YAdj^2 - YAdj S2]."
        ),
        "eliminated_kernel":(
            "K_S = (mu0/M0^2) v0 v0^T + 21(muAdj/MAdj^2) vAdj vAdj^T, "
            "v0=(1,1), vAdj=(0,1)."
        ),
        "exact_equal_response_kernel":K_exact.tolist(),
        "exact_determinant":float(np.linalg.det(K_exact)),
        "exact_contractions":exact,
        "canonical_basis_invariance":basis,
        "one_loop_running":run,
        "one_loop_running_kernel":K_run.tolist(),
        "one_loop_contractions":run_contractions,
        "flavor_at_protected_Neff":{
            "N_eff":float(flavor["N_eff"]),
            "max_error_pct":float(flavor["max_error_pct"]),
            "rms_error_pct":float(flavor["rms_error_pct"]),
        },
        "finite_momentum":{
            "response":"k_a(p)=mu_a/(M_a^2+Z_a p^2)",
            "exact_condition_all_p":"Z_a/M_a^2 equal across singlet and adjoint blocks",
            "wall_momentum_estimates":momenta,
            "interpretation":(
                "The pure auxiliary limit has no derivative contamination. A propagating composite is also safe if its "
                "dimensionless derivative ratio Z/M^2 is O(22)-universal. At wall momenta of a few times M_match, "
                "percent-level nonuniversality in Z/M^2 can move N_eff outside the flavor tolerance."
            ),
        },
        "claim_boundary":(
            "This is an exact local auxiliary/Legendre EFT matching and demonstrates that canonical normalization does not "
            "convert mu/M^2 into the excluded mu^2/M^2 kernel. It is not yet a fully dynamical non-supersymmetric composite "
            "bound-state calculation; such a completion must reproduce the same zero-momentum response and control Z/M^2."
        ),
        "main_conclusion":(
            "The protected coefficient survives the explicit auxiliary matching. The 21 Gram matrix, its inverse rational "
            "relations, and the 0.596% flavor result are unchanged at one-loop gauge order. The remaining risk is not hidden "
            "canonical normalization but derivative dynamics of a propagating composite; the strictly auxiliary realization "
            "avoids that risk, while a composite realization needs approximately universal Z/M^2 at the wall scale."
        ),
    }
    (out / "protected_auxiliary_matching_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__ == "__main__":
    main()
