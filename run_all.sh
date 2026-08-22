#!/usr/bin/env bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
python src/operator_origin_initial_audit.py >/dev/null
python src/exact_operator_and_two_spurion.py >/dev/null
python src/local_core_mass_extension.py >/dev/null
python src/raw_gradient_wilson_closure.py >/dev/null
python src/vectorlike_mediator_closure.py >/dev/null
python src/mediator_partner_spectrum.py >/dev/null
python src/integer_mediator_ps_audit.py >/dev/null
python src/explicit_ps_seed_lagrangian.py >/dev/null
python src/seed_sector_beta_function.py >/dev/null
python src/protected_auxiliary_composite_matching.py >/dev/null
python src/absolute_amplitude_closure.py >/dev/null
python src/residual_spurion_normalization.py >/dev/null
python src/radial_mode_seed_matching.py >/dev/null
python src/explicit_radial_uv_operator_matching.py >/dev/null
python src/geometric_modulus_embedding_audit.py >/dev/null
python src/constrained_internal_response_metric.py >/dev/null
pytest -q
