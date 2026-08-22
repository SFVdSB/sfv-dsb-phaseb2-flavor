import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUMMARY = json.loads((ROOT / 'results/absolute_amplitude_closure/absolute_amplitude_closure_summary.json').read_text())


def test_zero_fit_flavor_below_one_percent():
    assert SUMMARY['zero_fit_flavor']['max_error_pct'] < 1.0


def test_amplitudes_match_four_parameter_solution():
    assert max(abs(v) for v in SUMMARY['baseline_amplitude_relative_errors_pct'].values()) < 0.1


def test_common_attenuation_clusters():
    assert SUMMARY['shared_attenuation']['relative_span_pct'] < 0.05


def test_multiplicity_attenuation_matches_cluster():
    pred = SUMMARY['shared_attenuation']['predicted_exp_minus_2eps_over_15']
    mean = SUMMARY['shared_attenuation']['mean']
    assert abs(pred / mean - 1.0) * 100 < 0.05


def test_chiral_spectrum_passes():
    spec = SUMMARY['chiral_spectrum']
    assert spec['desired_profiles_with_exactly_one_near_zero'] == 9
    assert spec['opposite_profiles_with_zero_near_zero'] == 9
    assert spec['minimum_opposite_eigenvalue'] > 0.8


def test_representation_ratios_are_protected():
    z = SUMMARY['zero_fit_flavor']['controls']
    assert abs(z[0] / z[5] - 22/21) < 1e-12
    assert abs(z[1] / z[5] - 23/21) < 1e-12
    assert abs(z[2] / z[6] - 1/21) < 1e-12


def test_claim_boundary_is_present():
    assert 'not yet' in SUMMARY['claim_boundary'].lower()
