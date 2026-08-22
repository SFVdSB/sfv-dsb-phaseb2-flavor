import json
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'seed_sector_beta'
S=json.loads((OUT/'seed_sector_beta_summary.json').read_text())

def test_equal_running_protects_mu_over_m2():
    for b in S['blocks']:
        assert abs(b['mu_over_M2_ratio']-1.0) < 2e-9

def test_derivative_kernel_enhances_and_fails():
    row=next(x for x in S['scenario_results'] if x['scenario']=='composite derivative susceptibility')
    assert row['N_eff'] > 21.0
    assert row['max_flavor_error_pct'] > 1.0

def test_protected_scenarios_keep_exact_21():
    for name in ['O22 boundary imposed at matching scale','source-to-mass-squared protected ratio','dimension-six normalized kernel']:
        row=next(x for x in S['scenario_results'] if x['scenario']==name)
        assert abs(row['N_eff']-21.0) < 1e-8

def test_scale_window_is_short_for_derivative_route():
    assert S['maximum_scale_ratio_for_sub1pct']['derivative_enhancing_mu2_over_M2'] < 1.25

def test_outputs_exist():
    assert (OUT/'one_loop_seed_running_by_block.csv').exists()
    assert len(pd.read_csv(OUT/'effective_kernel_scenarios.csv')) >= 6
