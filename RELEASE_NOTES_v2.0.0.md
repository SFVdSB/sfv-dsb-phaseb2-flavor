# SFV/dSB Phase B2 v2.0.0 Release Notes

## Classification

**Major interpretive correction / conditional theorem closure.**

This release does not revise the published Phase-B2 numerical flavor benchmark. It corrects the spacetime interpretation and freezes a dimensionally consistent intrinsic-3+1 flavor theorem.

## Unchanged from v1.7.0

- corrected O(4) bounce inputs;
- Pati-Salam `15+3+3=21` seed multiplicity;
- protected `[[1,1],[1,22]]` Gram matrix;
- relations `22/21`, `23/21`, `1/21`;
- `Sigma_rad = 0.0670851132`;
- `Z0` and `ZF` response factors;
- seven-observable result: **0.6111% max**, **0.3580% RMS**;
- no continuous flavor amplitude adjusted at the frozen v1.7 stage.

## Corrected physical interpretation

The historical radial coordinate is an O(4) Euclidean relative coordinate, not an independent fifth spacetime direction. Quarks are intrinsic fields of the 3+1 daughter universe. The radial overlap is reclassified as a four-dimensional one-bounce transition/1PI flavor form factor.

## New conditional matrix-origin result

For a three-family S3 permutation carrier, the most general invariant primitive family tensor is

`P = a I + b J`.

Requiring one primitive collective family channel (rank one) forces `a=0`, yielding the normalized family-singlet projector

`P1 = J/3`.

The Phase-B2 direction `(-1,0,+1)` lies orthogonal to the singlet and supplies bounce-conditioned family breaking. The corrected kernels are

`Ku = Nu ∫ dr H DQ P1 Du`

`Kd = Nd ∫ dr H DQ P1 Dd`.

The counterfactual `P=I` gives simultaneously diagonal kernels and trivial CKM; `P=P1` plus distinct up/down radial dressings permits full-rank noncommuting kernels.

## New conditional daughter-inheritance result

The Phase-B2 response coordinates are interpreted as global finite-dimensional variables of the reduced one-bounce matching functional, with saddle

`A* = Sigma_B`, `B* = -Sigma_B`.

After matching they are constants of the daughter EFT, not local memory fields required to track the scalar after the wall passes.

## Retired interpretations

- independent fifth-dimensional flavor localization;
- permanent 2+1 wall quarks;
- 2+1-to-3+1 holographic flavor lift;
- ordinary primordial high-scale Yukawa-boundary interpretation;
- local late-time scalar memory as the source of present-day flavor.

## Claim boundary

v2.0.0 is a **conditional effective closure**. It does not derive the three-family carrier or reduced response sector uniquely from the unextended two-scalar action, and the historical architecture was developed with knowledge of flavor data.

## Relationship to v1.7.0

The v1.7.0 publication remains valid as the historical numerical/operator-origin record. v2.0.0 should be cited for the corrected spacetime embedding and daughter-flavor theorem.
