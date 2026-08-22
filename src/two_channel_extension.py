#!/usr/bin/env python3
"""Independent extension 1: Hessian-rotating two-channel fermion model.

The local eigenbasis is fixed by the scalar Hessian angle.  The fermion mass
matrix is
    M_f = U_H diag(q O, q O + beta Delta_H) U_H^T,
with Delta_H=sqrt(lambda_hard-lambda_soft)/alpha.  The script solves the exact
first-order two-channel zero-mode equation and reports the induced light-channel
core correction.  The 51-wall audit is optional because the compact release does
not bundle all full-resolution wall profiles.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from chiral_localization import CoordinateMap, load_profile, derive_baseline_normalization, build_wall_basis, normalized_zero_mode_from_B


def hessian(P,q,p):
    a=p['lambda_Phi']*(3*P**2-p['rho']**2)+2*p['bias']+2*p['g']*q**2
    d=p['lambda_brane']*(3*q**2-1)-2*p['mu2_tilde']+2*p['g']*P**2
    b=4*p['g']*P*q
    disc=np.sqrt((a-d)**2+4*b*b)
    return a,d,b,.5*(a+d-disc),.5*(a+d+disc)


def solve_case(profile: pd.DataFrame, meta: dict, normalization, beta=1.0, spacing=.005):
    cmap=CoordinateMap()
    basis=build_wall_basis(profile,cmap,24.0,spacing,normalization)
    r=profile.r.to_numpy(float); P=profile.Phi.to_numpy(float); ph=profile.phi.to_numpy(float)
    p={k:float(meta[k]) for k in ['lambda_Phi','rho','bias','g','lambda_brane','mu2_tilde']}
    a,d,b,lo,hi=hessian(P,ph,p)
    theta=.5*np.unwrap(np.arctan2(2*b,a-d))
    xd=cmap.alpha*(r-float(meta['R_peak']))
    th=np.interp(basis.x,xd,theta,left=theta[0],right=theta[-1])
    tp=np.gradient(th,basis.x)
    los=np.interp(basis.x,xd,lo,left=lo[0],right=lo[-1])
    his=np.interp(basis.x,xd,hi,left=hi[0],right=hi[-1])
    gap=beta*np.sqrt(np.maximum(his-los,1e-12))/cmap.alpha
    Osp=CubicSpline(basis.x,basis.O); tps=CubicSpline(basis.x,tp); gps=CubicSpline(basis.x,gap)
    def rhs(x,y):
        t=float(tps(x)); ml=float(Osp(x)); mh=ml+float(gps(x))
        return [-ml*y[0]+t*y[1],-t*y[0]-mh*y[1]]
    sol=solve_ivp(rhs,(basis.x[0],basis.x[-1]),[1.,0.],t_eval=basis.x,method='DOP853',rtol=1e-10,atol=1e-12)
    if not sol.success: raise RuntimeError(sol.message)
    light,heavy=sol.y
    ratio=np.divide(heavy,light,out=np.zeros_like(heavy),where=np.abs(light)>1e-250)
    correction=-tp*ratio
    core=np.abs(basis.x)<=6
    k=float(np.dot(correction[core],basis.E[core])/np.dot(basis.E[core],basis.E[core]))
    corr=float(np.corrcoef(correction[core],basis.E[core])[0,1])
    amp=np.sqrt(light*light+heavy*heavy); amp/=np.sqrt(np.trapezoid(amp*amp,basis.x))
    scalar,_=normalized_zero_mode_from_B(basis.x,basis.O+k*basis.E)
    return {
        'beta':beta,'k_geo':k,'correlation_with_E':corr,
        'relative_shape_residual':float(np.linalg.norm(correction[core]-k*basis.E[core])/np.linalg.norm(k*basis.E[core])),
        'heavy_norm_fraction':float(np.trapezoid(heavy*heavy,basis.x)/np.trapezoid(light*light+heavy*heavy,basis.x)),
        'scalar_profile_overlap':float(np.trapezoid(amp*scalar,basis.x)),
        'gap_min':float(gap.min()),'gap_at_core':float(gap[np.argmin(abs(basis.x))]),
    }


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',type=Path,default=ROOT); ap.add_argument('--output',type=Path,default=None)
    a=ap.parse_args(); root=a.root.resolve(); out=a.output or root/'results/independent_extensions/two_channel_beta_scan.json'
    profile=load_profile(root/'data/background_profile_O4_regular_robin_full.csv'); cmap=CoordinateMap(); norm=derive_baseline_normalization(profile,cmap)
    meta={'lambda_Phi':.1,'rho':2.357142857,'bias':0.,'g':2.313019,'lambda_brane':1e-8,'mu2_tilde':1e-8,'R_peak':float(profile.R_peak.iloc[0])}
    rows=[solve_case(profile,meta,norm,b) for b in [.25,.5,1.,2.,4.]]
    out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(rows,indent=2)+'\n'); print(json.dumps(rows,indent=2))
if __name__=='__main__': main()
