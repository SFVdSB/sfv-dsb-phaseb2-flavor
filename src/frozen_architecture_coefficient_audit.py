#!/usr/bin/env python3
"""Freeze the successful local operator and audit coefficient closure.

This script does not claim a first-principles flavor derivation. It tests a
predeclared exploratory four-fixed/three-fitted relation and computes the local
identifiability matrix around the exact seven-control solution.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import least_squares

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from chiral_localization import (CoordinateMap,load_profile,derive_baseline_normalization,
 build_wall_basis,higgs_profile,overlap_matrix_trapezoid,flavor_observables,route1_parameters)

OBS=['ct','ut','sb','db','Vus','Vcb','Vub']
NAMES=['h_Q','h_u0','h_u1','a_d0','a_d1','h_d0','h_d1']
EXACT=np.array([2.9748097239308926,3.191115303149601,-0.03098638280151201,
                0.25265488283253024,0.23413864783043356,2.8951808354446698,
                -0.6850396342903919])

def build_context():
 p=load_profile(ROOT/'data/background_profile_O4_regular_robin_full.csv')
 cmap=CoordinateMap(); norm=derive_baseline_normalization(p,cmap)
 b=build_wall_basis(p,cmap,24,.005,norm); H=higgs_profile(b.x,.514)
 C=b.dPhi_dy**2+b.dphi_dy**2; C/=C.max()
 tj=json.loads((ROOT/'configs/targets_MZ.json').read_text())['targets']
 target=np.array([tj[o] for o in OBS])
 return b,H,C,target

def values(z,b,H,C):
 pars=route1_parameters(np.asarray(z,float)); f={}
 for s in ['QL','uR','dR']:
  B=pars[s]['q'][:,None]*b.O[None,:]+(pars[s]['q']*pars[s]['h'])[:,None]*C[None,:]
  S=cumulative_trapezoid(B,b.x,axis=1,initial=0); L=-S; L-=L.max(axis=1)[:,None]
  ff=np.exp(L); ff/=np.sqrt(np.trapezoid(ff*ff,b.x,axis=1))[:,None]; f[s]=ff
 Yu=overlap_matrix_trapezoid(b.x,f['QL'],H,f['uR'])
 Yd=overlap_matrix_trapezoid(b.x,f['QL'],H,f['dR'])
 r=flavor_observables(Yu,Yd)
 return np.array([r['values'][o] for o in OBS])

def residual(z,b,H,C,target):
 v=values(z,b,H,C)
 return np.r_[np.log(v[:4]/target[:4]),(v[4:]-target[4:])/target[4:]]

def main():
 b,H,C,target=build_context()
 m=json.loads((ROOT/'data/baseline_phaseA_microphysics.json').read_text())
 eps=m['center_energy_excess_fraction']; k=m['k_geo']
 lock=b.alpha*(m['hessian_mixing_max_radius_dimless']-m['R_gradient_peak_dimless'])
 # Exploratory rational/geometric hypothesis. It is explicitly post-hoc.
 hQ=3-2*eps/5
 fixed={0:hQ,1:hQ+k+eps/4,4:lock*(.5-eps/3),5:hQ-2*k/5}
 free=[2,3,6]
 def zfrom(p):
  z=EXACT.copy()
  for i,v in fixed.items(): z[i]=v
  z[free]=p
  return z
 fit=least_squares(lambda p:residual(zfrom(p),b,H,C,target),EXACT[free],bounds=(-5,5),
                   max_nfev=1000,xtol=1e-13,ftol=1e-13,gtol=1e-13,x_scale='jac')
 z3=zfrom(fit.x); e3=100*(values(z3,b,H,C)/target-1)
 # Complete seven-formula exploratory map for a zero-fit diagnostic.
 z0=np.array([hQ,hQ+k+eps/4,-(1-eps)/(3*m['tachyonic_soft_end_dimless']),
              .25+eps/24,lock*(.5-eps/3),hQ-2*k/5,
              -(m['Phi_gradient_fraction']-eps/3)])
 e0=100*(values(z0,b,H,C)/target-1)
 # Logarithmic observable Jacobian and singular spectrum.
 v=values(EXACT,b,H,C); J=np.zeros((7,7))
 for j in range(7):
  h=1e-5*max(abs(EXACT[j]),1.0); zp=EXACT.copy(); zm=EXACT.copy(); zp[j]+=h; zm[j]-=h
  J[:,j]=(np.log(values(zp,b,H,C))-np.log(values(zm,b,H,C)))/(2*h)
 U,s,Vt=np.linalg.svd(J)
 outdir=ROOT/'results/coefficient_origin'; outdir.mkdir(parents=True,exist_ok=True)
 np.savetxt(outdir/'local_log_observable_jacobian.csv',J,delimiter=',',header=','.join(NAMES),comments='')
 result={
  'claim_boundary':'exploratory post-hoc coefficient hypothesis; not a first-principles derivation',
  'frozen_operator':'q*O + kappa*I_G plus separately established universal Hessian correction',
  'exact_controls':dict(zip(NAMES,EXACT.tolist())),
  'four_fixed':{NAMES[i]:float(v) for i,v in fixed.items()},
  'three_fitted':{NAMES[i]:float(v) for i,v in zip(free,fit.x)},
  'four_fixed_three_fit_errors_pct':dict(zip(OBS,e3.tolist())),
  'four_fixed_three_fit_max_error_pct':float(np.max(abs(e3))),
  'zero_fit_candidate_controls':dict(zip(NAMES,z0.tolist())),
  'zero_fit_candidate_errors_pct':dict(zip(OBS,e0.tolist())),
  'zero_fit_candidate_max_error_pct':float(np.max(abs(e0))),
  'jacobian_singular_values':s.tolist(),
  'jacobian_condition_number':float(s[0]/s[-1]),
  'weakest_parameter_direction':dict(zip(NAMES,Vt[-1].tolist())),
  'interpretation':'all seven directions are locally identifiable, but one combination is weak; the remaining unknown is a flavor/Wilson-coefficient matching law, not one missing scalar coefficient'
 }
 (outdir/'frozen_architecture_coefficient_summary.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__': main()
