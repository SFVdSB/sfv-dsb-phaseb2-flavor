# SFV/dSB Phase B2 v2.0.0

## O(4) bounce-encoded daughter flavor: dimensional correction and conditional closure

This repository contains the reproducible SFV/dSB Phase B2 quark-flavor construction and its v2.0.0 physical-embedding correction.

**The published Phase B2 v1.7 numerical benchmark is unchanged.** The v2.0.0 update corrects the interpretation of the historical radial flavor calculation and freezes a dimensionally consistent intrinsic-3+1 daughter-flavor theorem.

### Unchanged numerical result

The v1.7 construction remains frozen:

- complete Pati-Salam adjoint multiplicity `15+3+3=21`;
- protected Gram matrix `[[1,1],[1,22]]`;
- rational relations `22/21`, `23/21`, and `1/21`;
- bounce-derived radial driver
  `Sigma_B = 2(1-Phi_c/rho)^2(1+Phi_c/rho) = 0.0670851132`;
- response factors `Z0=exp(Sigma_B/21)` and `ZF=exp(-2 Sigma_B/15)`;
- seven-observable quark-flavor benchmark with **0.6111% maximum relative discrepancy** and **0.3580% RMS relative discrepancy**.

No observable target, wall input, response coefficient, or final benchmark value is retuned in v2.0.0.

## What v2.0.0 corrects

The historical one-dimensional flavor coordinate should **not** be interpreted as an independent fifth spacetime dimension. The corrected ontology is

```text
4D Euclidean O(4) bounce
        |
        +--> dS3 codimension-one wall (early cosmological state/support)
        |
        +--> intrinsic 3+1 daughter host (physical quark spacetime)
```

Quarks are intrinsic 3+1 daughter fields. The historical radial overlap is reinterpreted as a **four-dimensional one-bounce flavor form factor**, not as permanent extra-dimensional or 2+1 wall localization.

The frozen flavor statement is therefore:

```text
O(4) bounce-encoded radial flavor response
        -> branch-conditioned daughter Yukawa kernel
        -> intrinsic 3+1 quark flavor
```

## v2.0.0 matrix-origin theorem

Let

```text
P1 = (1/3) * 1 1^T
```

be the projector onto the symmetric family direction `1=(1,1,1)`.

For a three-family permutation carrier, S3 invariance gives the most general primitive tensor `P=a I+b J`. If the symmetric-stage Higgs sector contains only one primitive collective family channel (rank one), then `a=0`, giving uniquely, up to overall normalization,

```text
P = P1 = J/3.
```

The frozen Phase-B2 family-breaking direction `n=(-1,0,+1)` is orthogonal to the singlet and supplies the nontrivial bounce-conditioned radial family dressing.

The corrected flavor kernels are

```text
Ku = Nu ∫ dr H(r) DQ(r) P1 Du(r)
Kd = Nd ∫ dr H(r) DQ(r) P1 Dd(r)
```

with diagonal radial response matrices `DA(r)=diag(fA1,fA2,fA3)`.

With `P=I`, both kernels remain simultaneously diagonal and CKM is trivial. With `P=P1` and distinct up/down radial dressing, the integrated kernels can be full-rank and noncommuting, generating nontrivial CKM mixing.

## Daughter-EFT inheritance

The response variables of the augmented Phase-B2 completion are interpreted as **global reduced matching coordinates** of the one-bounce functional, not local fields that track the late-time scalar value.

Their reduced saddle is

```text
A* = Sigma_B
B* = -Sigma_B
```

and these values enter the branch-conditioned 3+1 daughter EFT as fixed coupling data after matching.

This is the precise conditional meaning of **bounce-encoded, daughter-inherited flavor**.

## Claim boundary

Phase B2 v2.0.0 is a **conditionally closed effective flavor theorem**, not a blind first-principles prediction of the unextended two-scalar action.

The remaining quarantined first-principles targets are:

1. microscopic origin of exactly three visible families / the permutation carrier;
2. derivation of the singlet-only primitive Higgs channel rather than assuming the one-channel condition;
3. derivation of the reduced global response sector from the microscopic SFV/dSB action;
4. complete Euclidean-to-Lorentzian 1PI matching and counterterm/RG treatment;
5. absolute top and bottom Yukawa normalizations;
6. CKM CP phase and later lepton-sector extension.

Future FP work is not permitted to retune the frozen v1.7 numerical architecture merely to make those derivations easier.

## Superseded interpretations

The following interpretations are retired:

- an independent fifth-dimensional flavor-localization coordinate;
- permanent quark localization on a 2+1 dS3 wall;
- a 2+1 quark-flavor theory holographically lifted to 3+1;
- treating the Phase-B2 matrix as an ordinary primordial `10^13-10^14 GeV` Yukawa boundary condition;
- treating the late-time local scalar value as a persistent flavor-memory source.

## Paper

- `paper/phaseB2_v1.7.0_paper.tex` / `.pdf` - historical numerical/operator-origin publication
- `paper/phaseB2_v2.0.0_correction.tex` - v2 dimensional correction and conditional daughter-flavor theorem

## Reproduce the unchanged v1.7 numerical benchmark

```bash
python -m pip install -r requirements.txt
python src/constrained_internal_response_metric.py
pytest -q
```

Expected headline output remains:

```text
Sigma_rad = 0.06708511319891537
max flavor error = 0.6110760482679778 %
74 tests passed
```

## Publication identifiers

- GitHub: https://github.com/SFVdSB/sfv-dsb-phaseb2-flavor
- Phase B2 v1.7 Zenodo DOI: https://doi.org/10.5281/zenodo.22059294
- Phase B2 v2.0 Zenodo DOI: pending new-version deposit

## License

See `LICENSE`.
