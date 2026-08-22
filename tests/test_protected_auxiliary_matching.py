import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from protected_auxiliary_composite_matching import (
    kernel_from_weights, sector_contractions, canonical_basis_invariance,
    protected_running_kernel, finite_momentum_neff
)


def test_exact_gram_and_determinant():
    K=kernel_from_weights(1,1,1,1)
    assert np.allclose(K,[[1,1],[1,22]],rtol=0,atol=1e-13)
    assert abs(np.linalg.det(K)-21)<1e-12


def test_exact_sector_ratios():
    d=sector_contractions(kernel_from_weights(1,1,1,1))
    r=d['ratios_to_down']
    assert abs(r['Q_L']-22/21)<1e-12
    assert abs(r['u_R']-23/21)<1e-12
    assert abs(r['d_R']-1)<1e-12
    assert abs(d['odd_up_over_down']-1/21)<1e-12


def test_basis_invariance():
    t=canonical_basis_invariance(kernel_from_weights(1,1,1,1),ntrial=50)
    assert t['max_abs_cross_term_error']<1e-11


def test_one_loop_protected_N():
    _,d=protected_running_kernel()
    assert abs(d['N_eff']-21)<1e-10
    assert all(abs(v-1)<1e-10 for v in d['block_k_over_k0'].values())


def test_finite_momentum_universal_derivative_ratio():
    for q in (0,0.5,1,3,10):
        assert abs(finite_momentum_neff(q,1,1,1)-21)<1e-12


def test_summary_exists_after_reproduction():
    p=ROOT/'results/protected_auxiliary_matching/protected_auxiliary_matching_summary.json'
    if p.exists():
        s=json.loads(p.read_text())
        assert s['flavor_at_protected_Neff']['max_error_pct']<1.0
