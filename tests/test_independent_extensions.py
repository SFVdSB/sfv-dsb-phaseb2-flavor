import json
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
RES=ROOT/'results/independent_extensions'

def test_two_channel_baseline_and_crosswall():
    s=json.loads((RES/'two_channel_extension_summary.json').read_text())
    b=s['baseline']; c=s['crosswall_51']
    assert 0.18 < b['k_geo'] < 0.22
    assert 0.04 < b['heavy_norm_fraction'] < 0.07
    assert b['scalar_profile_overlap'] > 0.998
    assert c['count']==51
    assert c['k_geo_relative_span'] < 0.09

def test_one_extra_spurion_fails_but_projectors_succeed():
    d=pd.read_csv(RES/'representation_spurion_summary.csv')
    one=d[(d.model=='Dq_Yk_Dslope7') & (d.bound_or_constraint=='abs_coeff_le_5')].iloc[0]
    full=d[(d.model=='projectors7') & (d.bound_or_constraint=='q_min_ge_0.5_and_abs_coeff_le_5')].iloc[0]
    assert one.max_error_pct > 1.0
    assert full.max_error_pct < 1e-5

def test_local_core_success_and_chirality():
    d=pd.read_csv(RES/'local_core_invariant_summary.csv')
    raw=d[d.core=='raw_local_gradient'].iloc[0]
    assert raw.max_error_pct < 1e-5
    assert raw.kappa_min > -5 and raw.kappa_max < 5
    s=json.loads((RES/'local_core_partner_spectrum_summary.json').read_text())
    assert s['desired_profiles_with_exactly_one_near_zero']==9
    assert s['opposite_profiles_with_zero_near_zero']==9
    assert s['minimum_opposite_eigenvalue'] > 0.5
