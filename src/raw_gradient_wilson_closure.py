#!/usr/bin/env python3
"""Raw-gradient Wilson-coefficient closure for Phase B2 v0.5.0.

The local core used in v0.3-v0.4 was numerically normalized as I_G=G/Gmax.
This script restores the raw local invariant
    G=(dPhi/dy)^2+(dphi/dy)^2
and reports the underlying Wilson coefficients c=kappa/Gmax.  It then tests a
frozen rational/geometric benchmark formula and audits how many coefficients
can be fixed while retaining <1% flavor accuracy.
"""
from __future__ import annotations
import itertools,json,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy.integrate import cumulative_trapezoid,solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from chiral_localization import (CoordinateMap,load_profile,derive_baseline_normalization,
 build_wall_basis,higgs_profile,overlap_matrix_trapezoid,flavor_observables,route1_parameters)
OBS=['ct','ut','sb','db','Vus','Vcb','Vub']
NAMES=['h_Q','h_u0','h_u1','a_d0','a_d1','h_d0','h_d1']
SECTORS=['QL','uR','dR']

def context():
 p=load_profile(ROOT/'data/background_profile_O4_regular_robin_full.csv')
 cmap=CoordinateMap();norm=derive_baseline_normalization(p,cmap)
 b=build_wall_basis(p,cmap,24,.01,norm);H=higgs_profile(b.x,.514)
 G=b.dPhi_dy**2+b.dphi_dy**2;Gmax=float(G.max());C=G/Gmax
 # exact universal two-channel envelope
 h=pd.read_csv(ROOT/'data/baseline_local_hessian_and_mixing.csv')
 r=h.r.to_numpy(float);R=float(p.R_peak.iloc[0]);xh=cmap.alpha*(r-R)
 th=.5*np.unwrap(np.arctan2(2*h.H_cross.to_numpy(float),h.H_PhiPhi.to_numpy(float)-h.H_phiphi.to_numpy(float)))
 th=np.interp(b.x,xh,th,left=th[0],right=th[-1]);tp=np.gradient(th,b.x)
 lo=np.interp(b.x,xh,h.H_eigen_soft);hi=np.interp(b.x,xh,h.H_eigen_hard)
 gap=np.sqrt(np.maximum(hi-lo,1e-12))/cmap.alpha
 tps=CubicSpline(b.x,tp);gps=CubicSpline(b.x,gap)
 def rhs(x,y):
  t=float(tps(x));d=float(gps(x));return [t*y[1],-t*y[0]-d*y[1]]
 sol=solve_ivp(rhs,(b.x[0],b.x[-1]),[1.,0.],t_eval=b.x,method='DOP853',rtol=1e-11,atol=1e-13)
 if not sol.success:raise RuntimeError(sol.message)
 env=np.sqrt(sol.y[0]**2+sol.y[1]**2);env/=env.max()
 targetj=json.loads((ROOT/'configs/targets_MZ.json').read_text())['targets']
 target=np.array([targetj[o] for o in OBS])
 return b,H,G,C,Gmax,env,target

def evaluate(z,b,H,C,env,target):
 pars=route1_parameters(np.asarray(z,float));f={}
 for s in SECTORS:
  q=pars[s]['q'];hh=pars[s]['h']
  B=q[:,None]*(b.O[None,:]+hh[:,None]*C[None,:])
  S=cumulative_trapezoid(B,b.x,axis=1,initial=0);L=-S;L-=L.max(axis=1)[:,None]
  ff=np.exp(L)*env[None,:];ff/=np.sqrt(np.trapezoid(ff*ff,b.x,axis=1))[:,None];f[s]=ff
 Yu=overlap_matrix_trapezoid(b.x,f['QL'],H,f['uR'])
 Yd=overlap_matrix_trapezoid(b.x,f['QL'],H,f['dR'])
 rr=flavor_observables(Yu,Yd);v=np.array([rr['values'][o] for o in OBS])
 e=100*(v/target-1)
 return v,e,rr

def residual(z,b,H,C,env,target):
 v,_,_=evaluate(z,b,H,C,env,target)
 return np.r_[np.log(v[:4]/target[:4]),(v[4:]-target[4:])/target[4:]]

def main():
 b,H,G,C,Gmax,env,target=context()
 # Refit exact local+geometric benchmark.
 seed=np.array([2.7925,2.9120,-.03160,.25964,.24228,2.6640,-.65899])
 fit=least_squares(lambda z:residual(z,b,H,C,env,target),seed,bounds=(-8,8),max_nfev=1800,
                   xtol=1e-13,ftol=1e-13,gtol=1e-13,x_scale='jac')
 exact=fit.x;_,exact_err,_=evaluate(exact,b,H,C,env,target)
 # Restore raw-gradient Wilson coefficients.
 raw={
  'c_Q':exact[0]/Gmax,
  'c_u0':exact[1]/Gmax,
  'c_u1':exact[2]/Gmax,
  'a_d0':exact[3],
  'a_d1':exact[4],
  'c_d0':exact[5]/Gmax,
  'c_d1':exact[6]/Gmax,
 }
 m=json.loads((ROOT/'data/baseline_phaseA_microphysics.json').read_text())
 lock=b.alpha*(m['hessian_mixing_max_radius_dimless']-m['R_gradient_peak_dimless'])
 formula_raw={'c_Q':11/3,'c_u0':23/6,'c_u1':-1/24,
              'a_d0':m['m_true_Phi_dimless']/4,'a_d1':lock/2,
              'c_d0':7/2,'c_d1':-7/8}
 z0=np.array([formula_raw['c_Q']*Gmax,formula_raw['c_u0']*Gmax,
              formula_raw['c_u1']*Gmax,formula_raw['a_d0'],formula_raw['a_d1'],
              formula_raw['c_d0']*Gmax,formula_raw['c_d1']*Gmax])
 _,zero_err,_=evaluate(z0,b,H,C,env,target)
 # Audit every fixed subset; formula values are held fixed, remaining controls fitted.
 rows=[]
 labels=['c_Q','c_u0','c_u1','a_d0','a_d1','c_d0','c_d1']
 for nf in range(8):
  for fixed in itertools.combinations(range(7),nf):
   free=[i for i in range(7) if i not in fixed]
   if free:
    def fun(u):
     z=z0.copy();z[free]=u;return residual(z,b,H,C,env,target)
    ft=least_squares(fun,exact[free],bounds=(-8,8),max_nfev=1200,
                     xtol=2e-12,ftol=2e-12,gtol=2e-12,x_scale='jac')
    z=z0.copy();z[free]=ft.x
   else:z=z0.copy()
   _,e,_=evaluate(z,b,H,C,env,target)
   rows.append({'n_fixed':nf,'fixed':'+'.join(labels[i] for i in fixed),
                'n_fitted':len(free),'max_error_pct':float(np.max(abs(e))),
                'rms_error_pct':float(np.sqrt(np.mean(e*e))),
                **dict(zip(NAMES,z.tolist()))})
 sdf=pd.DataFrame(rows)
 outdir=ROOT/'results/raw_gradient_wilson';outdir.mkdir(parents=True,exist_ok=True)
 sdf.to_csv(outdir/'raw_rational_formula_subset_audit.csv',index=False)
 best={}
 for nf in range(8):
  r=sdf[sdf.n_fixed==nf].sort_values(['max_error_pct','rms_error_pct']).iloc[0]
  zz=np.array([r[n] for n in NAMES],float)
  best[str(nf)]={'fixed':r.fixed,'n_fitted':int(r.n_fitted),
                 'max_error_pct':float(r.max_error_pct),'rms_error_pct':float(r.rms_error_pct),
                 'normalized_controls':dict(zip(NAMES,zz.tolist())),
                 'raw_parameters':{'c_Q':float(zz[0]/Gmax),'c_u0':float(zz[1]/Gmax),
                   'c_u1':float(zz[2]/Gmax),'a_d0':float(zz[3]),'a_d1':float(zz[4]),
                   'c_d0':float(zz[5]/Gmax),'c_d1':float(zz[6]/Gmax)}}
 # local observable Jacobian
 J=np.empty((7,7));eps=1e-5
 for j in range(7):
  dz=eps*max(1,abs(exact[j]));zp=exact.copy();zm=exact.copy();zp[j]+=dz;zm[j]-=dz
  J[:,j]=(residual(zp,b,H,C,env,target)-residual(zm,b,H,C,env,target))/(2*dz)
 _,sing,vt=np.linalg.svd(J,full_matrices=False)
 pd.DataFrame(J,index=OBS,columns=NAMES).to_csv(outdir/'local_geo_residual_jacobian.csv')
 result={
  'claim_boundary':'post-hoc rational/geometric hypothesis generated after inspecting the local solution; not yet first-principles or blind',
  'operator':'q O + c G_raw plus the separately calculated universal two-channel envelope',
  'Gmax_y':Gmax,
  'exact_normalized_controls':dict(zip(NAMES,exact.tolist())),
  'exact_errors_pct':dict(zip(OBS,exact_err.tolist())),
  'exact_raw_wilson_coefficients':raw,
  'frozen_raw_formula':formula_raw,
  'formula_relative_errors_vs_exact_raw_pct':{k:100*(formula_raw[k]/raw[k]-1) for k in raw},
  'zero_fit_errors_pct':dict(zip(OBS,zero_err.tolist())),
  'zero_fit_max_error_pct':float(np.max(abs(zero_err))),
  'best_by_number_fixed':best,
  'five_fixed_two_fitted':best['5'],
  'four_fixed_three_fitted':best['4'],
  'six_fixed_one_fitted':best['6'],
  'jacobian_singular_values':sing.tolist(),
  'jacobian_condition_number':float(sing[0]/sing[-1]),
  'weakest_normalized_control_direction':dict(zip(NAMES,vt[-1].tolist())),
  'interpretation':'normalizing the local gradient core hid nearly constant raw Wilson coefficients; the benchmark admits five frozen relations with two fitted controls below one percent, but the relations remain post-hoc'
 }
 (outdir/'raw_gradient_wilson_summary.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
