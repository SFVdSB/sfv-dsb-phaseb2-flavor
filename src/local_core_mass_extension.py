#!/usr/bin/env python3
"""Independent extension 3: local wall-core mass invariants.

Fits the unchanged seven-control Route-I sector structure after replacing the
symmetrized reporting mode E(x) by genuinely local wall invariants.  This tests
operator naturalness, not flavor-parameter compression.
"""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from chiral_localization import CoordinateMap,load_profile,derive_baseline_normalization,build_wall_basis,higgs_profile,overlap_matrix_trapezoid,flavor_observables,route1_parameters
OBS=['ct','ut','sb','db','Vus','Vcb','Vub'];SECTORS=['QL','uR','dR']

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT);ap.add_argument('--output',type=Path,default=None);a=ap.parse_args();root=a.root.resolve();out=a.output or root/'results/independent_extensions/local_core_primary_refit.json'
 profile=load_profile(root/'data/background_profile_O4_regular_robin_full.csv');cmap=CoordinateMap();norm=derive_baseline_normalization(profile,cmap);basis=build_wall_basis(profile,cmap,24,.005,norm)
 targetj=json.loads((root/'configs/targets_MZ.json').read_text());target=np.array([targetj['targets'][o] for o in OBS]);H=higgs_profile(basis.x,.514)
 C=basis.dPhi_dy**2+basis.dphi_dy**2;C/=C.max()
 z0=np.array([2.9748097239308926,3.191115303149601,-.03098638280151201,.25265488283253024,.23413864783043356,2.8951808354446698,-.6850396342903919])
 def evaluate(z):
  pars=route1_parameters(z);f={}
  for s in SECTORS:
   B=pars[s]['q'][:,None]*(basis.O[None,:]+pars[s]['h'][:,None]*C[None,:]);S=cumulative_trapezoid(B,basis.x,axis=1,initial=0);L=-S;L-=L.max(axis=1)[:,None];ff=np.exp(L);ff/=np.sqrt(np.trapezoid(ff*ff,basis.x,axis=1))[:,None];f[s]=ff
  Yu=overlap_matrix_trapezoid(basis.x,f['QL'],H,f['uR']);Yd=overlap_matrix_trapezoid(basis.x,f['QL'],H,f['dR']);r=flavor_observables(Yu,Yd);v=np.array([r['values'][o] for o in OBS]);return v,r
 def fun(z):
  v,_=evaluate(z);return np.r_[np.log(v[:4]/target[:4]),(v[4:]-target[4:])/target[4:]]
 fit=least_squares(fun,z0,bounds=(-5,5),max_nfev=1600,xtol=1e-13,ftol=1e-13,gtol=1e-13,x_scale='jac');v,r=evaluate(fit.x);err=100*(v/target-1);pars=route1_parameters(fit.x);k=np.concatenate([pars[s]['q']*pars[s]['h'] for s in SECTORS])
 result={'core':'raw_local_gradient_density','controls':fit.x.tolist(),'errors_pct':err.tolist(),'max_error_pct':float(max(abs(err))),'kappa_min':float(k.min()),'kappa_max':float(k.max()),'condition_Yu':r['condition_Yu'],'condition_Yd':r['condition_Yd']}
 out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
