import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/explicit_ps_seed_lagrangian'
SUMMARY=json.loads((OUT/'explicit_ps_seed_lagrangian_summary.json').read_text())

def test_ps_adjoint_dimension_is_21():
    assert SUMMARY['group_data']['dim_adjoint_total']==21
    assert SUMMARY['group_data']['seed_space_dimension_1_plus_adj']==22

def test_component_casimir_is_not_mistaken_for_21():
    assert SUMMARY['group_data']['component_level_C2_for_(4,2)']=='21/8'

def test_seed_gram_matrix_exact():
    g=np.array(SUMMARY['gram_matrix'],float)
    assert np.allclose(g,[[1,1],[1,22]],rtol=0,atol=1e-12)
    assert abs(np.linalg.det(g)-21)<1e-10

def test_exact_sector_ratios():
    d=SUMMARY['exact_N21_contractions']
    r=d['ratios_to_down']
    assert abs(r['Q_L']-22/21)<1e-12
    assert abs(r['u_R']-23/21)<1e-12
    assert abs(r['d_R']-1)<1e-12
    assert abs(d['odd_up_over_down']-1/21)<1e-12

def test_continuous_tolerance_contains_equal_threshold():
    q=SUMMARY['continuous_flavor_tolerance_below_1pct']
    assert q['exists']
    assert q['N_eff_min'] <= 21 <= q['N_eff_max']
    assert q['kAdj_over_k0_min'] <= 1 <= q['kAdj_over_k0_max']

def test_continuous_scan_reaches_subpercent():
    df=pd.read_csv(OUT/'continuous_Neff_flavor_tolerance.csv')
    row=df.iloc[(df.N_eff-21).abs().argmin()]
    assert row.max_error_pct < 0.61
    assert df.max_error_pct.min() < 0.56

def test_incidence_has_one_singlet_and_21_adjoint_rows():
    df=pd.read_csv(OUT/'seed_incidence_matrix_22x2.csv')
    assert len(df)==22
    vals=df[['coupling_to_S1','coupling_to_S2']].to_numpy(float)
    assert np.allclose(vals[0],[1,1])
    assert np.allclose(vals[1:],[0,1])
