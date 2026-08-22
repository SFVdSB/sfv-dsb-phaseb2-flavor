import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def test_coefficient_checkpoint():
 d=json.loads((ROOT/'results/coefficient_origin/frozen_architecture_coefficient_summary.json').read_text())
 assert d['four_fixed_three_fit_max_error_pct'] < 0.35
 assert d['zero_fit_candidate_max_error_pct'] < 2.2
 assert len(d['jacobian_singular_values']) == 7
 assert min(d['jacobian_singular_values']) > 0
 assert d['jacobian_condition_number'] > 1000
