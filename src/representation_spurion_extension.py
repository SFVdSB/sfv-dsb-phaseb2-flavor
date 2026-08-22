from __future__ import annotations
import sys,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from exact_operator_and_two_spurion import evaluate_q_kappa,residual_from_values
from chiral_localization import CoordinateMap,load_profile,derive_baseline_normalization,build_wall_basis
OBS=['ct','ut','sb','db','Vus','Vcb','Vub']; S=['QL','uR','dR']; n=np.array([-1.,0.,1.]);Y={'QL':1/6,'uR':2/3,'dR':-1/3};D={'QL':0.,'uR':0.,'dR':1.};U={'QL':0.,'uR':1.,'dR':0.};C={'QL':0.,'uR':1.,'dR':-1.}
profile=load_profile(ROOT/'data/background_profile_O4_regular_robin_full.csv');cmap=CoordinateMap();norm=derive_baseline_normalization(profile,cmap);basis=build_wall_basis(profile,cmap,24,.01,norm)
tj=json.load(open(ROOT/'configs/targets_MZ.json'));target=np.array([tj['targets'][o] for o in OBS])

def make(model,p):
 q={};k={}
 if model=='Dq_Yk6':
  aD,aFD,b0,bF,bY,bFY=p
  for s in S:q[s]=np.exp(D[s]*(aD+aFD*n));k[s]=b0+bF*n+bY*Y[s]+bFY*Y[s]*n
 elif model=='Dq_Yk_Dint7':
  aD,aFD,b0,bF,bY,bFY,bD=p
  for s in S:q[s]=np.exp(D[s]*(aD+aFD*n));k[s]=b0+bF*n+bY*Y[s]+bFY*Y[s]*n+bD*D[s]
 elif model=='Dq_Yk_Dslope7':
  aD,aFD,b0,bF,bY,bFY,bFD=p
  for s in S:q[s]=np.exp(D[s]*(aD+aFD*n));k[s]=b0+bF*n+bY*Y[s]+bFY*Y[s]*n+bFD*D[s]*n
 elif model=='Dq_Yk_Cint7':
  aD,aFD,b0,bF,bY,bFY,bC=p
  for s in S:q[s]=np.exp(D[s]*(aD+aFD*n));k[s]=b0+bF*n+bY*Y[s]+bFY*Y[s]*n+bC*C[s]
 elif model=='Dq_Yk_Cslope7':
  aD,aFD,b0,bF,bY,bFY,bFC=p
  for s in S:q[s]=np.exp(D[s]*(aD+aFD*n));k[s]=b0+bF*n+bY*Y[s]+bFY*Y[s]*n+bFC*C[s]*n
 elif model=='projectors7':
  aD,aFD,bF,bU,bFU,bD,bFD=p
  for s in S:q[s]=np.exp(D[s]*(aD+aFD*n));k[s]=bF*n+bU*U[s]+bFU*U[s]*n+bD*D[s]+bFD*D[s]*n
 else:raise KeyError(model)
 return q,k

def fit(model,dim,bound):
 def fun(p):
  try:
   q,k=make(model,p);v,_,_=evaluate_q_kappa(basis,q,k,target);r=residual_from_values(v,target);return r if np.all(np.isfinite(r)) else np.ones(7)*100
  except Exception:return np.ones(7)*100
 rng=np.random.default_rng(20260718+dim+int(bound*10))
 seeds=[np.zeros(dim)]
 seeds += [rng.uniform(-1,1,dim) for _ in range(8)]
 best=None
 for seed in seeds:
  ft=least_squares(fun,seed,bounds=(-bound,bound),max_nfev=900,xtol=2e-12,ftol=2e-12,gtol=2e-12,x_scale='jac')
  q,k=make(model,ft.x);v,e,res=evaluate_q_kappa(basis,q,k,target);cost=np.sum(ft.fun**2)
  if best is None or cost<best[0]:best=(cost,ft.x,v,e,q,k,res)
 cost,p,v,e,q,k,res=best
 flatq=np.concatenate([q[s] for s in S]);flatk=np.concatenate([k[s] for s in S]);h=np.concatenate([k[s]/q[s] for s in S])
 return dict(model=model,parameters=dim,bound=bound,cost=cost,max_error_pct=float(np.max(abs(e))),rms_error_pct=float(np.sqrt(np.mean(e*e))),parameter_values=p.tolist(),errors=e.tolist(),q_min=float(flatq.min()),q_max=float(flatq.max()),kappa_min=float(flatk.min()),kappa_max=float(flatk.max()),h_min=float(h.min()),h_max=float(h.max()),active_bound=bool(np.any(np.isclose(abs(p),bound,atol=1e-5))))
rows=[]
mods=[('Dq_Yk6',6),('Dq_Yk_Dint7',7),('Dq_Yk_Dslope7',7),('Dq_Yk_Cint7',7),('Dq_Yk_Cslope7',7),('projectors7',7)]
for m,d in mods:
 for b in [3,5]:
  r=fit(m,d,b);rows.append(r);print(m,b,r['max_error_pct'],r['rms_error_pct'],r['active_bound'],r['parameter_values'])
open(ROOT/'results/independent_extensions/representation_spurion_refit_full.json','w').write(json.dumps(rows,indent=2))
pd.DataFrame([{k:v for k,v in r.items() if k not in ['parameter_values','errors']} for r in rows]).to_csv(ROOT/'results/independent_extensions/representation_spurion_refit_full.csv',index=False)
