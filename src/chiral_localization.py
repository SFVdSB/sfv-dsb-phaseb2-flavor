#!/usr/bin/env python3
"""Core routines for SFV/dSB chiral Dirac localization.

The canonical scalar Route-I construction is
    B_Ai(x) = q_Ai [O(x) + h_Ai E(x)]
with x=y-y_H.  The physical mass is +B for Q_L and -B for right-handed
singlets.  In either case the desired zero-mode amplitude solves f'=-B f.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Tuple

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, quad, solve_ivp
from scipy.interpolate import CubicSpline
from scipy.special import gamma


@dataclass(frozen=True)
class CoordinateMap:
    y_H: float = 1.02
    s_H: float = 0.514
    R_peak: float = 5.860202508437851
    w_FWHM: float = 1.7517899927739728
    alpha: float = 0.6909375570964031


@dataclass(frozen=True)
class BasisNormalization:
    Phi_left: float
    Phi_right: float
    E_scale: float


@dataclass
class WallBasis:
    x: np.ndarray
    y: np.ndarray
    O: np.ndarray
    E: np.ndarray
    Phi: np.ndarray
    phi: np.ndarray
    dPhi_dy: np.ndarray
    dphi_dy: np.ndarray
    normalization: BasisNormalization
    R_peak: float
    alpha: float


def load_profile(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"r", "Phi", "Phi_prime", "phi", "phi_prime"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Profile missing columns: {sorted(missing)}")
    return df


def _interp_extended(x_data: np.ndarray, values: np.ndarray, x: np.ndarray,
                     left: float | None = None, right: float | None = None) -> np.ndarray:
    if left is None:
        left = float(values[0])
    if right is None:
        right = float(values[-1])
    return np.interp(x, x_data, values, left=left, right=right)


def derive_baseline_normalization(profile: pd.DataFrame, cmap: CoordinateMap,
                                  probe_half_width: float = 32.0,
                                  probe_points: int = 12801) -> BasisNormalization:
    x = np.linspace(-probe_half_width, probe_half_width, probe_points)
    R = float(profile["R_peak"].iloc[0]) if "R_peak" in profile else cmap.R_peak
    x_data = cmap.alpha * (profile["r"].to_numpy(float) - R)
    dPhi_dy = profile["Phi_prime"].to_numpy(float) / cmap.alpha
    dphi_dy = profile["phi_prime"].to_numpy(float) / cmap.alpha
    g = dPhi_dy**2 + dphi_dy**2
    gx = _interp_extended(x_data, g, x, left=0.0, right=0.0)
    gmx = _interp_extended(x_data, g, -x, left=0.0, right=0.0)
    E_scale = float(np.max(0.5 * (gx + gmx)))
    if not np.isfinite(E_scale) or E_scale <= 0:
        raise ValueError("Invalid gradient-mode normalization")
    return BasisNormalization(
        Phi_left=float(profile["Phi"].iloc[0]),
        Phi_right=float(profile["Phi"].iloc[-1]),
        E_scale=E_scale,
    )


def build_wall_basis(profile: pd.DataFrame, cmap: CoordinateMap,
                     half_width: float = 24.0, spacing: float = 0.005,
                     normalization: BasisNormalization | None = None) -> WallBasis:
    n = int(round(2 * half_width / spacing)) + 1
    x = np.linspace(-half_width, half_width, n)
    R = float(profile["R_peak"].iloc[0]) if "R_peak" in profile else cmap.R_peak
    x_data = cmap.alpha * (profile["r"].to_numpy(float) - R)

    Phi_data = profile["Phi"].to_numpy(float)
    phi_data = profile["phi"].to_numpy(float)
    dPhi_data = profile["Phi_prime"].to_numpy(float) / cmap.alpha
    dphi_data = profile["phi_prime"].to_numpy(float) / cmap.alpha

    if normalization is None:
        normalization = derive_baseline_normalization(profile, cmap)

    Phi = _interp_extended(x_data, Phi_data, x)
    phi = _interp_extended(x_data, phi_data, x)
    dPhi = _interp_extended(x_data, dPhi_data, x, left=0.0, right=0.0)
    dphi = _interp_extended(x_data, dphi_data, x, left=0.0, right=0.0)

    Phi_m = _interp_extended(x_data, Phi_data, -x)
    denom = normalization.Phi_left - normalization.Phi_right
    if abs(denom) < 1e-14:
        raise ValueError("Degenerate Phi vacuum normalization")
    T = 1.0 - 2.0 * (Phi - normalization.Phi_right) / denom
    Tm = 1.0 - 2.0 * (Phi_m - normalization.Phi_right) / denom
    O = 0.5 * (T - Tm)

    grad_data = dPhi_data**2 + dphi_data**2
    grad = _interp_extended(x_data, grad_data, x, left=0.0, right=0.0)
    grad_m = _interp_extended(x_data, grad_data, -x, left=0.0, right=0.0)
    E = 0.5 * (grad + grad_m) / normalization.E_scale

    return WallBasis(
        x=x,
        y=x + cmap.y_H,
        O=O,
        E=E,
        Phi=Phi,
        phi=phi,
        dPhi_dy=dPhi,
        dphi_dy=dphi,
        normalization=normalization,
        R_peak=R,
        alpha=cmap.alpha,
    )


def basis_diagnostics(basis: WallBasis) -> Dict[str, float]:
    O_ref = np.interp(-basis.x, basis.x, basis.O)
    E_ref = np.interp(-basis.x, basis.x, basis.E)
    odd_error = float(np.max(np.abs(basis.O + O_ref)))
    even_error = float(np.max(np.abs(basis.E - E_ref)))
    core = np.abs(basis.x) <= 5.0
    mat = np.column_stack([basis.O[core], basis.E[core]])
    s = np.linalg.svd(mat, compute_uv=False)
    condition = float(s[0] / s[-1])
    corr = float(np.corrcoef(basis.O[core], basis.E[core])[0, 1])
    return {
        "odd_parity_max_abs_error": odd_error,
        "even_parity_max_abs_error": even_error,
        "O_left": float(basis.O[0]),
        "O_right": float(basis.O[-1]),
        "E_left": float(basis.E[0]),
        "E_right": float(basis.E[-1]),
        "E_max": float(np.max(basis.E)),
        "basis_condition_number_core": condition,
        "basis_correlation_core": corr,
    }


def higgs_profile(x: np.ndarray, s_H: float = 0.514, center_x: float = 0.0) -> np.ndarray:
    return np.pi**(-0.25) * s_H**(-0.5) * np.exp(-(x - center_x)**2 / (2 * s_H**2))


def cumulative_integral_from_left(x: np.ndarray, values: np.ndarray) -> np.ndarray:
    return cumulative_trapezoid(values, x, initial=0.0)


def normalized_zero_mode_from_B(x: np.ndarray, B: np.ndarray) -> Tuple[np.ndarray, float]:
    """Return desired mode f solving f'=-B f and log normalization offset."""
    S = cumulative_integral_from_left(x, B)
    logf = -S
    shift = float(np.max(logf))
    f = np.exp(logf - shift)
    norm = float(np.sqrt(np.trapezoid(f * f, x)))
    if not np.isfinite(norm) or norm <= 0:
        raise FloatingPointError("Zero-mode normalization failed")
    return f / norm, shift + np.log(norm)


def raw_mode_with_center_value_one(x: np.ndarray, B: np.ndarray, sign: int = -1) -> np.ndarray:
    """Return exp(sign*integral_0^x B) with f(0)=1, stabilized for finite domains."""
    S_left = cumulative_integral_from_left(x, B)
    i0 = int(np.argmin(np.abs(x)))
    S0 = S_left[i0]
    exponent = sign * (S_left - S0)
    # Retain physical finite-domain growth while avoiding overflow.
    return np.exp(np.clip(exponent, -700.0, 700.0))


def solve_zero_mode_ode(x: np.ndarray, B: np.ndarray,
                        rtol: float = 1e-10, atol: float = 1e-12) -> np.ndarray:
    """Independent solve of f'=-B(x)f from x=0 in both directions."""
    spline = CubicSpline(x, B, bc_type="natural")
    i0 = int(np.argmin(np.abs(x)))
    x0 = float(x[i0])
    xr = x[i0:]
    xl_desc = x[:i0 + 1][::-1]

    def rhs(t, y):
        return [-float(spline(t)) * y[0]]

    sol_r = solve_ivp(rhs, (x0, float(x[-1])), [1.0], t_eval=xr,
                      rtol=rtol, atol=atol, method="DOP853")
    sol_l = solve_ivp(rhs, (x0, float(x[0])), [1.0], t_eval=xl_desc,
                      rtol=rtol, atol=atol, method="DOP853")
    if not sol_r.success or not sol_l.success:
        raise RuntimeError(f"ODE solve failed: {sol_l.message}; {sol_r.message}")
    left = sol_l.y[0][::-1]
    right = sol_r.y[0][1:]
    f = np.concatenate([left, right])
    norm = float(np.sqrt(np.trapezoid(f * f, x)))
    return f / norm


def density_moments(x: np.ndarray, f: np.ndarray) -> Dict[str, float]:
    p = f * f
    norm = float(np.trapezoid(p, x))
    p = p / norm
    mean = float(np.trapezoid(x * p, x))
    var = float(np.trapezoid((x - mean)**2 * p, x))
    sigma = float(np.sqrt(max(var, 0.0)))
    skew = float(np.trapezoid(((x - mean) / sigma)**3 * p, x)) if sigma > 0 else 0.0
    peak = float(x[int(np.argmax(f))])
    half = 0.5 * float(np.max(f))
    mask = f >= half
    amp_fwhm = float(x[mask][-1] - x[mask][0]) if np.any(mask) else float("nan")
    halfp = 0.5 * float(np.max(p))
    maskp = p >= halfp
    density_fwhm = float(x[maskp][-1] - x[maskp][0]) if np.any(maskp) else float("nan")
    return {
        "norm": norm,
        "mean_x": mean,
        "sigma_density": sigma,
        "equivalent_gaussian_amplitude_s": np.sqrt(2.0) * sigma,
        "skewness_density": skew,
        "peak_x": peak,
        "amplitude_FWHM": amp_fwhm,
        "density_FWHM": density_fwhm,
    }


def count_zero_crossings(x: np.ndarray, values: np.ndarray, zero_tol: float = 1e-10) -> Tuple[int, list[float]]:
    roots: list[float] = []
    for i in range(len(x) - 1):
        a, b = values[i], values[i + 1]
        if abs(a) <= zero_tol:
            roots.append(float(x[i]))
        elif a * b < 0:
            t = abs(a) / (abs(a) + abs(b))
            roots.append(float(x[i] * (1 - t) + x[i + 1] * t))
    if abs(values[-1]) <= zero_tol:
        roots.append(float(x[-1]))
    unique: list[float] = []
    for r in roots:
        if not unique or abs(r - unique[-1]) > 2 * (x[1] - x[0]):
            unique.append(r)
    return len(unique), unique


def route1_parameters(z: np.ndarray) -> Dict[str, Dict[str, np.ndarray]]:
    """Expand the canonical 7-control vector into q_i and h_i arrays."""
    if len(z) != 7:
        raise ValueError("Canonical Route-I vector must have seven controls")
    hQ, hu0, hu1, ad0, ad1, hd0, hd1 = map(float, z)
    n = np.array([-1.0, 0.0, 1.0])
    return {
        "QL": {"q": np.ones(3), "h": n * hQ},
        "uR": {"q": np.ones(3), "h": hu0 + n * hu1},
        "dR": {"q": np.exp(ad0 + n * ad1), "h": hd0 + n * hd1},
    }


def build_route1_profiles(basis: WallBasis, z: np.ndarray) -> Dict[str, Dict[str, np.ndarray]]:
    pars = route1_parameters(z)
    out: Dict[str, Dict[str, np.ndarray]] = {}
    for sector, values in pars.items():
        q, h = values["q"], values["h"]
        B = q[:, None] * (basis.O[None, :] + h[:, None] * basis.E[None, :])
        F = np.vstack([normalized_zero_mode_from_B(basis.x, row)[0] for row in B])
        out[sector] = {"q": q, "h": h, "B": B, "profiles": F}
    return out


def overlap_matrix_trapezoid(x: np.ndarray, left: np.ndarray, H: np.ndarray,
                              right: np.ndarray) -> np.ndarray:
    return np.trapezoid(left[:, None, :] * H[None, None, :] * right[None, :, :], x, axis=2)


def overlap_matrix_adaptive(x: np.ndarray, left: np.ndarray, H: np.ndarray,
                            right: np.ndarray, epsabs: float = 1e-11,
                            epsrel: float = 1e-10) -> np.ndarray:
    spl_L = [CubicSpline(x, row, bc_type="natural") for row in left]
    spl_R = [CubicSpline(x, row, bc_type="natural") for row in right]
    spl_H = CubicSpline(x, H, bc_type="natural")
    Y = np.empty((left.shape[0], right.shape[0]))
    for i in range(left.shape[0]):
        for j in range(right.shape[0]):
            val, _ = quad(lambda t: float(spl_L[i](t) * spl_H(t) * spl_R[j](t)),
                          float(x[0]), float(x[-1]), epsabs=epsabs, epsrel=epsrel,
                          limit=400)
            Y[i, j] = val
    return Y


def flavor_observables(Yu: np.ndarray, Yd: np.ndarray) -> Dict[str, object]:
    Uu, su, Vtu = np.linalg.svd(Yu)
    Ud, sd, Vtd = np.linalg.svd(Yd)
    iu = np.argsort(su)
    id_ = np.argsort(sd)
    su = su[iu]
    sd = sd[id_]
    Uu = Uu[:, iu]
    Ud = Ud[:, id_]
    Vtu = Vtu[iu, :]
    Vtd = Vtd[id_, :]
    Vckm_signed = Uu.T @ Ud
    Vckm = np.abs(Vckm_signed)
    values = {
        "ct": float(su[1] / su[2]),
        "ut": float(su[0] / su[2]),
        "sb": float(sd[1] / sd[2]),
        "db": float(sd[0] / sd[2]),
        "Vus": float(Vckm[0, 1]),
        "Vcb": float(Vckm[1, 2]),
        "Vub": float(Vckm[0, 2]),
    }
    return {
        "values": values,
        "up_singular_values": su,
        "down_singular_values": sd,
        "Uu": Uu,
        "Ud": Ud,
        "Vtu": Vtu,
        "Vtd": Vtd,
        "V_CKM_signed": Vckm_signed,
        "V_CKM_abs": Vckm,
        "condition_Yu": float(np.linalg.cond(Yu)),
        "condition_Yd": float(np.linalg.cond(Yd)),
    }


def evaluate_route1(basis: WallBasis, z: np.ndarray, s_H: float = 0.514,
                    center_x: float = 0.0, adaptive: bool = False) -> Dict[str, object]:
    model = build_route1_profiles(basis, z)
    H = higgs_profile(basis.x, s_H=s_H, center_x=center_x)
    overlap = overlap_matrix_adaptive if adaptive else overlap_matrix_trapezoid
    Yu = overlap(basis.x, model["QL"]["profiles"], H, model["uR"]["profiles"])
    Yd = overlap(basis.x, model["QL"]["profiles"], H, model["dR"]["profiles"])
    obs = flavor_observables(Yu, Yd)
    return {"model": model, "H": H, "Yu": Yu, "Yd": Yd, "observables": obs}


def normalized_gaussian(x: np.ndarray, center: float, s: float) -> np.ndarray:
    return np.pi**(-0.25) * s**(-0.5) * np.exp(-(x - center)**2 / (2 * s**2))


def gaussian_benchmark_observables(x: np.ndarray, geometry: dict) -> Dict[str, object]:
    yH = float(geometry["H"]["yH"])
    sH = float(geometry["H"]["sH"])
    H = normalized_gaussian(x, 0.0, sH)
    arrays = {}
    for sector in ("QL", "uR", "dR"):
        arrays[sector] = np.vstack([
            normalized_gaussian(x, float(y - yH), float(s))
            for y, s in geometry[sector]["ys"]
        ])
    Yu = overlap_matrix_trapezoid(x, arrays["QL"], H, arrays["uR"])
    Yd = overlap_matrix_trapezoid(x, arrays["QL"], H, arrays["dR"])
    return {"Yu": Yu, "Yd": Yd, "profiles": arrays, "H": H,
            "observables": flavor_observables(Yu, Yd)}


def analytic_kink_profile(x: np.ndarray, amplitude: float, inverse_width: float) -> np.ndarray:
    p = amplitude / inverse_width
    norm_sq = np.sqrt(np.pi) * gamma(p) / (inverse_width * gamma(p + 0.5))
    return np.cosh(inverse_width * x)**(-p) / np.sqrt(norm_sq)


def analytic_kink_validation(half_width: float = 20.0, spacing: float = 0.001,
                             amplitude: float = 1.3, inverse_width: float = 0.8) -> Dict[str, float]:
    n = int(round(2 * half_width / spacing)) + 1
    x = np.linspace(-half_width, half_width, n)
    B = amplitude * np.tanh(inverse_width * x)
    fq, _ = normalized_zero_mode_from_B(x, B)
    fo = solve_zero_mode_ode(x, B)
    fe = analytic_kink_profile(x, amplitude, inverse_width)
    fe = fe / np.sqrt(np.trapezoid(fe * fe, x))
    return {
        "quadrature_L2_error": float(np.sqrt(np.trapezoid((fq - fe)**2, x))),
        "ode_L2_error": float(np.sqrt(np.trapezoid((fo - fe)**2, x))),
        "quadrature_vs_ode_L2_error": float(np.sqrt(np.trapezoid((fq - fo)**2, x))),
        "quadrature_max_abs_error": float(np.max(np.abs(fq - fe))),
        "ode_max_abs_error": float(np.max(np.abs(fo - fe))),
        "analytic_norm_on_domain": float(np.trapezoid(fe * fe, x)),
        "amplitude": amplitude,
        "inverse_width": inverse_width,
        "half_width": half_width,
        "spacing": spacing,
    }
