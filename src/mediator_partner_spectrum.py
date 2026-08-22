#!/usr/bin/env python3
"""Partner-operator spectrum for compressed vectorlike-mediator models."""
from __future__ import annotations
import json,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from raw_gradient_wilson_closure import context
from chiral_localization import route1_parameters

def main():
 b,H,G,C,Gmax,env,target=context()
 summary=json.loads((ROOT/'results/mediator_closure/vectorlike_mediator_summary.json').read_text())
 chosen=['singlet_flavon_ratio21_6p','clebsch_5p','two_amplitude_core_clebsch_4p']
 loge=np.log(np.maximum(env,1e-300));Bgeo=-np.gradient(loge,b.x,edge_order=2)
 x=np.linspace(-24,24,801);dx=x[1]-x[0]
 rows=[];summaries={}
 for model in chosen:
  z=np.array(summary['models'][model]['effective_controls'],float);pars=route1_parameters(z)
  start=len(rows)
  for sector in ['QL','uR','dR']:
   for j in range(3):
    B0=pars[sector]['q'][j]*(b.O+pars[sector]['h'][j]*C)+Bgeo
    B=np.interp(x,b.x,B0);dB=np.gradient(B,dx,edge_order=2)
    for op,Vall in [('desired',B*B-dB),('opposite',B*B+dB)]:
     V=Vall[1:-1];main=2/dx**2+V;off=np.full(len(V)-1,-1/dx**2)
     ham=diags([off,main,off],[-1,0,1],format='csr')
     ev=np.sort(eigsh(ham,k=4,which='SA',return_eigenvectors=False,tol=1e-7,maxiter=50000))
     threshold=min(B[0]**2,B[-1]**2)
     row={'model':model,'sector':sector,'generation':j+1,'operator':op,
          'near_zero_count_abs_lt_1e3':int(np.sum(abs(ev)<1e-3)),
          'below_threshold_count':int(np.sum(ev<threshold-1e-4)),
          'asymptotic_threshold':float(threshold)}
     row.update({f'eigenvalue_{i}':float(v) for i,v in enumerate(ev)});rows.append(row)
  frame=pd.DataFrame(rows[start:]);des=frame[frame.operator=='desired'];opp=frame[frame.operator=='opposite']
  summaries[model]={'desired_profiles_with_exactly_one_near_zero':int(np.sum(des.near_zero_count_abs_lt_1e3==1)),
                    'opposite_profiles_with_zero_near_zero':int(np.sum(opp.near_zero_count_abs_lt_1e3==0)),
                    'minimum_opposite_eigenvalue':float(opp.eigenvalue_0.min()),
                    'minimum_nonzero_desired_eigenvalue':float(des.eigenvalue_1.min()),'total_profiles':9}
 out=ROOT/'results/mediator_closure';pd.DataFrame(rows).to_csv(out/'mediator_partner_spectrum.csv',index=False)
 (out/'mediator_partner_spectrum_summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
 print(json.dumps(summaries,indent=2))
if __name__=='__main__':main()
