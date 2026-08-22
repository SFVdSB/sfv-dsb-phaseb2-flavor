# SFV/dSB Microphysical Closure — Phase B v0.1.0

This repository tests whether the seven successful Route-I chiral-localization controls can be replaced by functions of the Phase-A wall microphysics and a smaller set of fermion/flavor inputs.

## Scientific result

The result is **partial closure**, not a first-principles flavor prediction.

- The scalar wall fixes the odd transition mode `O(x)`, the even gradient-core mode `E(x)`, and several geometric landmarks.
- A robust relation is found for the down-sector charge spacing:

  `a_d1 ≈ alpha * (R_max_mixing - R_total_gradient_peak)`.

  Across 51 nearby wall solutions, this relation has correlation `0.96784`, mean absolute relative error `0.633%`, and maximum error `2.91%`.
- Fixing only this relation and fitting the remaining six controls reproduces all seven flavor observables within `0.0691%`.
- An exploratory seven-formula map built only from Phase-A quantities reproduces the canonical benchmark with **zero continuous flavor fit** and maximum observable error `0.6005%`.
- The zero-fit formula model preserves nine desired chiral zero modes, excludes all nine opposite-chirality zero modes, and has minimum opposite-sector eigenvalue `0.65055`.
- A model with four controls fixed by Phase-A formulas and only three controls fitted reproduces all seven observables within `0.2745%`.
- Simple universal one-slope or shared-charge alternatives fail, with best tested maximum error `10.94%`.

## Claim boundary

The seven-formula map was constructed after inspecting the successful seven-control benchmark. It is therefore a **frozen candidate ansatz**, not a prediction. The 51-wall test supports the locking-interval formula for `a_d1`, but it does not independently validate all seven formulas. Counterfactual wall changes, especially in the stiffness coordinate `Y`, produce large flavor shifts.

The scalar wall cannot by itself select generation-dependent fermion charges or Wilson coefficients. A flavor symmetry, charge quantization rule, or explicit ultraviolet fermion sector remains necessary for full first-principles closure.

## Exploratory Phase-A formula map

For the canonical wall:

```text
h_Q  = 3 / (2 I_Phi)
h_u0 = sqrt(pi) / (2 rho)
h_u1 = -1 / (3 r_tachyonic_end)
a_d0 = xi_true_Phi / 3
a_d1 = alpha (R_max_mixing - R_total_gradient_peak)
h_d0 = (5/3) epsilon_center
h_d1 = -(1/2) f_gradient_phi
```

where `I_Phi` is the integrated bulk-gradient weight, `epsilon_center` is the fractional central energy excess above the true vacuum, and `f_gradient_phi` is the brane share of the integrated gradient energy.

## Reproduce the benchmark analysis

```bash
python -m pip install -r requirements.txt
python src/phaseB_closure_analysis.py
python src/make_plots.py
pytest -q
```

The included all-wall tables are frozen extended-audit products. Regenerating them requires the full Part-I profile archive and the Phase-A dictionary; see `src/refit_controls_all_walls.py` and `src/evaluate_formula_all_walls.py`.

## Important outputs

- `docs/phaseB_microphysical_closure_results.tex`
- `results/phaseB_summary.json`
- `results/baseline_microphysical_control_map.csv`
- `results/constrained_model_results.csv`
- `results/all_wall_control_relation_summary.csv`
- `results/refit_controls_all_walls.csv`
- `results/formula_closure_all_walls.csv`
- `results/formula4_blind_tests.csv`
- `results/posthoc_formula_partner_spectrum.csv`
- `results/compact_charge_model_results.csv`

## Recommended next step

Freeze the candidate formula map before any new observable is inspected. Then introduce the smallest explicit flavor-spurion or charge-symmetry model capable of deriving the remaining generation dependence. Test that frozen structure against an external observable sector, rather than adding further benchmark-specific formulas.
