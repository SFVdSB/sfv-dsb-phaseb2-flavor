# SFV/dSB Phase B2 v1.7.0

## Pati-Salam flavor operator origin and constrained response-metric completion

This repository is the reproducible Phase B2 v1.7.0 release of the SFV/dSB quark-flavor construction. It starts from the frozen corrected O(4) wall background and develops an explicit operator-origin chain for the wall-localized flavor coefficients.

### Main result

The release combines:

- an exact local reduction of the chiral wall operator;
- a complete Pati-Salam adjoint seed with multiplicity `15+3+3=21`;
- a protected zero-momentum auxiliary response with Gram matrix `[[1,1],[1,22]]`;
- the rational flavor relations `22/21`, `23/21`, and `1/21`;
- a bounce-derived radial driver
  `Sigma_rad = 2(1-Phi_c/rho)^2(1+Phi_c/rho) = 0.0670851132`;
- a constrained internal-response metric with generator
  `T = diag(I21/21,-I15/15)`;
- exact response factors `Z0=exp(Sigma/21)` and `ZF=exp(-2 Sigma/15)`;
- a frozen seven-observable quark-flavor realization with **0.6111% maximum relative error** and **0.3580% RMS relative error**, with no continuous flavor amplitude adjusted at the v1.7.0 stage.

The full regression suite contains **74 tests**, all passing in the publication audit.

### Claim boundary

This release is deliberately **not** presented as proof that the original two-scalar SFV/dSB action uniquely forces the complete flavor sector. The determinant-preserving response compensator and its locking invariant are one additional structural principle in the augmented effective action. The historical development also used the flavor benchmark while discovering the architecture. Accordingly, the defensible claim is a reproducible zero-continuous-fit realization / minimal completion, not a blind prediction from the unextended action.

### Paper

- `paper/phaseB2_v1.7.0_paper.tex`
- `paper/phaseB2_v1.7.0_paper.pdf`

### Reproduce the v1.7 checkpoint

```bash
python -m pip install -r requirements.txt
python src/constrained_internal_response_metric.py
pytest -q
```

Expected headline outputs:

```text
Sigma_rad = 0.06708511319891537
max flavor error = 0.6110760482679778 %
74 tests passed
```

### Repository structure

- `src/` — analysis and operator-matching code
- `tests/` — regression tests
- `data/` — frozen wall/background and microphysical inputs
- `results/` — frozen numerical outputs and audits
- `docs/` — checkpoint reports documenting the derivation chain
- `paper/` — LaTeX manuscript and compiled preprint
- `historical/` — earlier partial-closure package retained for provenance only

### Historical Phase B v0.1.0

The older microphysical-closure package is included under `historical/phaseB-microphysical-closure-v0.1.0/`. It is **not required** to reproduce v1.7.0 and should not be interpreted as a second current model. It is retained because it records the earlier post-hoc/partial-closure stage and makes the progression to Phase B2 auditable.

### License

See `LICENSE`. The retained project license allows inspection and noncommercial scholarly reproduction with attribution and reserves redistribution/commercial use unless permission is granted.
