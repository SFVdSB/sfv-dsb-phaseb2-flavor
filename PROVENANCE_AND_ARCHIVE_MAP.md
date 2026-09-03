# Flavor provenance and archive map

**Status:** Gate-A1 provenance map  
**Last audited:** 2026-09-03

This file distinguishes the numerical/operator-origin publication from the later dimensional/host correction. The two layers are related, but they do not make identical claims.

| Layer | Identifier | Classification | Scientific role |
|---|---|---|---|
| Earlier Phase B microphysical closure | `historical/phaseB-microphysical-closure-v0.1.0/` | Historical / superseded as current Flavor authority | Earlier partial/post-hoc closure preserved for provenance |
| Internal Phase-B2 operator-origin work | local `sfv-dsb-fermion-wall-operator-origin-phaseB2-v1.7.0` archive | Historical/current derivation parent for v1.7; not itself the exact Zenodo publication package | Detailed operator-origin/numerical work; 216 SHA-256 entries verified in Gate A1 |
| Public Phase B2 v1.7.0 | GitHub commit `675b64905183efbcc7b187fc9360d86d8846d389`; Zenodo `10.5281/zenodo.22059294` | **Historical numerical/operator-origin authority** | Frozen `15+3+3=21`, protected Gram structure, bounce response, and seven-observable benchmark |
| Internal FLAVOR-HOST-0Z freeze | `SFV-dSB-FLAVOR-HOST-0Z-FROZEN-v1.0.0` | **Current internal dimensional/host theorem parent** | Post-publication audit that retires the extra-dimensional interpretation and freezes the intrinsic-3+1 daughter theorem; no numerical retuning |
| Public Phase B2 v2.0.0 | GitHub `main`; Zenodo `10.5281/zenodo.22113999` | **Current public corrected-interpretation layer** | Dimensional correction, `P1=J/3` conditional matrix-origin theorem, and reduced one-bounce daughter-EFT inheritance |

## Frozen numerical continuity

Phase B2 v2.0.0 does **not** replace the v1.7 numerical benchmark. It retains without retuning:

- Pati-Salam adjoint multiplicity `15+3+3=21`;
- protected Gram matrix `[[1,1],[1,22]]`;
- ratios `22/21`, `23/21`, `1/21`;
- `Sigma_B = 0.0670851132`;
- `Z0=exp(Sigma_B/21)` and `ZF=exp(-2 Sigma_B/15)`;
- seven-observable result: `0.6111%` maximum relative discrepancy and `0.3580%` RMS.

Thus the authority split is:

> **v1.7.0 = numerical/operator-origin record**  
> **v2.0.0 = corrected physical embedding and conditional daughter-flavor theorem**

## HOST-0Z relationship

FLAVOR-HOST-0Z was frozen against public Phase-B2 commit

`675b64905183efbcc7b187fc9360d86d8846d389`

and DOI

`10.5281/zenodo.22059294`.

Its central correction is that the historical radial flavor coordinate is the relative radial coordinate of the four-dimensional Euclidean O(4) saddle, not an independent fifth spacetime dimension. Quarks are intrinsic `3+1` daughter fields. The historical radial overlap is reinterpreted as a one-bounce form factor / reduced matching kernel.

The internal HOST freeze also makes explicit two conditional ingredients:

1. the three-family permutation carrier plus one primitive collective Higgs channel leading to `P1=J/3`;
2. the response coordinates as global reduced one-bounce matching variables whose saddle values become daughter-EFT coupling data.

The public v2.0.0 GitHub commits after the v1.7 baseline publish this correction/theorem rather than altering the frozen numerical solver.

## Retired interpretations

The following are historical and must not be presented as the current physical interpretation:

- independent fifth-dimensional flavor localization;
- permanent quark localization on the `2+1 dS3` wall;
- a `2+1 -> 3+1` holographic flavor lift;
- ordinary primordial high-scale Yukawa boundary-condition interpretation;
- local late-time scalar-memory interpretation of the response coordinates.

## Zenodo status

- `10.5281/zenodo.22059294` remains citable as the historical v1.7 numerical/operator-origin publication.
- `10.5281/zenodo.22113999` is the current v2.0 correction record supplied for this audit and should be cited for the corrected embedding/theorem.
- The v2 release-preparation metadata states that v2 was intended as a new version of v1.7. Gate A1 still requires a file-by-file verification of the live `22113999` record and confirmation of the final Zenodo version relation/metadata.

## Current GitHub publication notes

The v2 correction manuscript source was committed before the Zenodo deposit and therefore retains pre-publication wording such as “should be deposited.” Do not silently rewrite the scientific paper merely to remove that historical wording. The repository README and citation metadata carry the current DOI.

No GitHub Release object was observed during Gate A1. Formal tags/releases can be addressed later in the publication-infrastructure pass after the repository/content audit is complete.

## Authority rules for Volkas/publication-readiness work

1. Use the frozen v1.7 solver/results for the numerical seven-observable benchmark and operator-origin history.
2. Use FLAVOR-HOST-0Z for the complete internal post-publication dimensional/host audit.
3. Use public v2.0 / DOI `10.5281/zenodo.22113999` for the corrected physical interpretation and daughter-flavor theorem.
4. When discussing `21=15+3+3`, preserve the distinction between the exact Pati-Salam adjoint multiplicity and the additional equal-per-component response/threshold condition needed to obtain the exact protected 21-channel kernel.
5. Do not describe v2.0 as a blind first-principles prediction of the unextended two-scalar action.
