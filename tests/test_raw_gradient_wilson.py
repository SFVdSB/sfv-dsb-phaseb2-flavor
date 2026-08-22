import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
RES=ROOT/'results/raw_gradient_wilson'
def test_raw_wilson_closure():
 d=json.loads((RES/'raw_gradient_wilson_summary.json').read_text())
 assert d['zero_fit_max_error_pct'] < 3.0
 assert d['five_fixed_two_fitted']['max_error_pct'] < 0.71
 assert d['four_fixed_three_fitted']['max_error_pct'] < 0.54
 assert d['six_fixed_one_fitted']['max_error_pct'] > 1.0
 assert d['jacobian_condition_number'] > 2000
 assert max(abs(x) for x in d['formula_relative_errors_vs_exact_raw_pct'].values()) < 1.6
def test_raw_formula_spectrum():
 d=json.loads((RES/'raw_formula_partner_spectrum_summary.json').read_text())
 for model in ['zero_fit_raw_rational','five_fixed_two_fit']:
  assert d[model]['desired_profiles_with_exactly_one_near_zero'] == 9
  assert d[model]['opposite_profiles_with_zero_near_zero'] == 9
  assert d[model]['minimum_opposite_eigenvalue'] > 0.8
def test_crosswall_local_stability():
 d=json.loads((RES/'raw_formula_crosswall_stability_summary.json').read_text())['local33']
 assert max(v['mean_abs_relative_error_pct'] for v in d.values()) < 1.5
 assert max(v['max_abs_relative_error_pct'] for v in d.values()) < 3.0
