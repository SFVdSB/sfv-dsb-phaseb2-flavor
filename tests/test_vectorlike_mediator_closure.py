import json
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]

def summary():
 return json.loads((ROOT/'results/mediator_closure/vectorlike_mediator_summary.json').read_text())

def test_rank_and_rank_one_failure():
 d=summary(); assert d['tree_level_matching']['h_matrix_rank']==2
 assert d['models']['one_rank_one_mediator_6p']['max_error_pct']>20

def test_compressed_benchmark_models():
 d=summary();m=d['models']
 assert m['singlet_flavon_ratio21_6p']['max_error_pct']<0.2
 assert m['clebsch_5p']['max_error_pct']<0.2
 assert m['two_amplitude_core_clebsch_4p']['max_error_pct']<1.0
 assert m['three_parameter_near_miss']['max_error_pct']>1.0

def test_crosswall_relation():
 d=summary()['crosswall_constrained_refits']
 assert d['ratio21_6p']['points']==51
 assert d['ratio21_6p']['points_below_1pct']==51
 assert d['ratio21_6p']['maximum_error_pct']<1.0

def test_partner_spectrum():
 d=json.loads((ROOT/'results/mediator_closure/mediator_partner_spectrum_summary.json').read_text())
 for v in d.values():
  assert v['desired_profiles_with_exactly_one_near_zero']==9
  assert v['opposite_profiles_with_zero_near_zero']==9
  assert v['minimum_opposite_eigenvalue']>0.5

def test_ratio_audit_written():
 df=pd.read_csv(ROOT/'results/mediator_closure/rational_ratio_crosswall_audit.csv')
 assert set(df.domain)=={'all_51','local_33','corridor_18'}
 assert len(df)==12
