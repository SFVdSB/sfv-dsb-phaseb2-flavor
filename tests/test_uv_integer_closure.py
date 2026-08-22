import json
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
SUMMARY=json.loads((ROOT/'results/uv_integer_closure/uv_integer_closure_summary.json').read_text())


def test_integer_matrix_exact_ratios():
    r=SUMMARY['route_1_integer_mediator']['generated_ratios']
    assert r['hQ_over_hd0']=='22/21'
    assert r['hu0_over_hd0']=='23/21'
    assert r['hu1_over_hd1']=='1/21'


def test_integer_matrix_is_positive():
    m=SUMMARY['route_1_integer_mediator']
    assert m['determinant']==21
    assert m['positive_definite']
    assert min(m['eigenvalues'])>0


def test_N21_unique_under_one_percent():
    s=SUMMARY['integer_scan_summary']
    assert s['best_integer_N']==21
    assert s['integers_below_1pct']==[21]
    assert s['best_max_error_pct']<1.0
    assert s['neighbor_N20_error_pct']>1.0
    assert s['neighbor_N22_error_pct']>1.0


def test_family_anomalies_cancel():
    a=SUMMARY['route_2_family_symmetry_and_anomalies']
    assert a['all_local_anomalies_zero']
    assert all(v=='0' for v in a['total_coefficients'].values())
    assert a['witten_anomalies_absent']


def test_pati_salam_pattern():
    p=SUMMARY['route_3_pati_salam_contractions']
    assert p['group_dimension']==21
    assert p['eight_times_C2_sum']==21
    assert p['exact_pattern']=={'Q_L':'22','u_R':'23','d_R':'21'}


def test_cross_project_consistency():
    c=SUMMARY['route_5_cross_project_consistency']
    assert c['Lorentz_architecture']['status']=='pass'
    assert c['strong_CP_architecture']['status']=='pass_conditionally'
    assert c['gauge_invariance']['status']=='pass_conditionally'


def test_scan_file_complete():
    df=pd.read_csv(ROOT/'results/uv_integer_closure/integer_N_flavor_scan_2_80.csv')
    assert len(df)==79
    assert df.N.min()==2 and df.N.max()==80
