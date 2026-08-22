# Publication guide — Phase B2 v1.7.0

## Recommended public identity

### GitHub repository name

`sfv-dsb-phaseb2-flavor`

### GitHub repository title / release title

**SFV/dSB Phase B2 v1.7.0 — Pati–Salam Flavor Operator Origin and Constrained Response-Metric Completion**

### GitHub “About” description

Reproducible Phase B2 v1.7.0 flavor construction for SFV/dSB: O(4)-wall chiral localization, Pati–Salam adjoint-seed origin of 21, protected auxiliary matching, a bounce-derived radial driver, and a constrained response metric yielding a zero-continuous-fit 0.6111% maximum discrepancy across seven quark-flavor observables.

### GitHub release description

Phase B2 v1.7.0 is the frozen reproducible release of the SFV/dSB quark-flavor operator-origin program. A complete Pati–Salam adjoint seed supplies the multiplicity 15+3+3=21. Protected auxiliary matching produces the exact Gram matrix [[1,1],[1,22]] and the rational relations 22/21, 23/21, and 1/21. The corrected O(4) bounce independently fixes Sigma_rad=0.0670851132. The v1.7.0 constrained internal-response action then yields Z0=exp(Sigma/21) and ZF=exp(-2 Sigma/15), preserving a frozen seven-observable result with 0.6111% maximum and 0.3580% RMS relative discrepancy. The preferred auxiliary realization adds no spacetime dimension and no propagating scalar. The repository includes source code, frozen inputs, numerical outputs, derivation reports, 74 passing regression tests, the LaTeX paper, and the earlier Phase B v0.1.0 package as historical provenance. The compensator locking invariant remains one additional structural principle beyond the original two-scalar SFV/dSB action, so this release is described as a zero-continuous-fit minimal completion rather than a blind prediction of the unextended action.

## Zenodo metadata

### Recommended record type

**Publication → Preprint** if the main object you want cited is the scientific paper and you want one DOI covering the paper plus its reproducibility archive.

This is preferable to the automatic GitHub-Zenodo software route if your goal is one publication DOI. The GitHub integration creates a software-oriented archive record. You can use that instead if you want the DOI to identify the software release.

### Zenodo title

**Pati–Salam Seed Multiplicity and a Constrained Internal-Response Metric for Quark Flavor in SFV/dSB: Phase B2 v1.7.0**

### Creator

Steven Hoffmann — Independent Researcher

Add your ORCID in Zenodo if desired.

### Publication date

2026-08-22

### Version

1.7.0

### Zenodo description / abstract

We present Phase B2 v1.7.0 of the SFV/dSB quark-flavor program, a reproducible construction built on a frozen corrected O(4) scalar-wall background. A complete adjoint seed of SU(4)_C × SU(2)_L × SU(2)_R supplies the exact multiplicity 15+3+3=21. In a protected zero-momentum auxiliary matching, one singlet seed plus the 21-component adjoint generates the Gram matrix [[1,1],[1,22]], whose inverse produces the rational flavor-sector relations 22/21, 23/21, and 1/21. The corrected bounce independently fixes the radial driver Sigma_rad=2(1-Phi_c/rho)^2(1+Phi_c/rho)=0.0670851132. A minimal determinant-preserving internal-response action with 21- and 15-dimensional isotropic blocks then uniquely produces the generator T=diag(I21/21,-I15/15) and the response factors Z0=exp(Sigma/21) and ZF=exp(-2 Sigma/15).

With the architecture and frozen inputs fixed, the resulting seven-observable quark-flavor benchmark has a maximum relative discrepancy of 0.6111% and an RMS discrepancy of 0.3580%, with no continuous flavor amplitude adjusted at the v1.7.0 stage. The preferred response variables are auxiliary and add neither a spacetime dimension nor a propagating scalar. The full repository contains frozen inputs, source code, numerical outputs, derivation reports, and 74 passing regression tests.

The claim boundary is explicit. The determinant-preserving response compensator and its locking invariant are one additional structural principle beyond the original two-scalar SFV/dSB action, and the historical model-development process inspected the flavor benchmark. The result is therefore presented as a zero-continuous-fit realization and minimal effective completion, not an independently blind first-principles prediction of the unextended action. The earlier Phase B v0.1.0 microphysical-closure package is included only as historical provenance.

### Keywords

SFV/dSB; quark flavor; Pati–Salam; chiral localization; O(4) bounce; CKM matrix; flavor hierarchy; effective field theory; domain wall; reproducible research

## Final publication identifiers

- GitHub: https://github.com/SFVdSB/sfv-dsb-phaseb2-flavor
- Zenodo DOI: https://doi.org/10.5281/zenodo.22059294

## What to upload

### GitHub

Upload the **unpacked publication folder**, not only the ZIP. GitHub should expose the README, paper, code, data, tests, and results as browsable files.

The publication package contains more than 100 files, so the GitHub browser uploader is not the best route. Use GitHub Desktop or Git instead.

The historical Phase B v0.1.0 folder is included under `historical/` because it is useful for provenance, but it is not required to reproduce v1.7.0.

### Zenodo — recommended single-record route

Upload only two files:

1. `phaseB2_v1.7.0_paper.pdf`
2. `sfv-dsb-phaseB2-v1.7.0-publication.zip`

Zenodo itself recommends ZIP packaging when an upload contains many files. The PDF is the primary preprint; the ZIP is the complete reproducibility archive.

## Minimal GitHub workflow

### Easiest: GitHub Desktop

1. Extract `sfv-dsb-phaseB2-v1.7.0-publication.zip` on your computer.
2. Open GitHub Desktop. Create or initialize a repository in that extracted folder, commit all files once, and click **Publish repository**. Make it public when ready.
3. On GitHub, create a release using tag `v1.7.0` and paste the release description above. GitHub automatically provides source ZIP/tar archives for the tagged repository.

### Command-line alternative

Create an empty repository on GitHub named `sfv-dsb-phaseb2-flavor`, then open a terminal inside the extracted publication folder and run:

```bash
git init
git add .
git commit -m "Publish SFV/dSB Phase B2 v1.7.0"
git branch -M main
git remote add origin https://github.com/SFVdSB/sfv-dsb-phaseb2-flavor.git
git push -u origin main
```

Then create the GitHub release/tag `v1.7.0` in the web interface.

## Minimal Zenodo workflow for one paper DOI

1. Create a **New upload** in Zenodo and choose **Publication → Preprint**.
2. Before final publication, click **Get a DOI now!** to reserve the DOI.
3. The final identifiers are already inserted in `paper/phaseB2_v1.7.0_paper.tex`:
   - GitHub: `https://github.com/SFVdSB/sfv-dsb-phaseb2-flavor`
   - Zenodo DOI: `10.5281/zenodo.22059294`
   Commit the final PDF/source before tagging `v1.7.0`.
4. Upload the final PDF and the final publication ZIP to Zenodo, paste the metadata above, set the PDF as the default preview, and publish.

This produces one Zenodo DOI for the paper/reproducibility package while GitHub remains the browsable living repository.

## Alternative: automatic GitHub-Zenodo software DOI

If you prefer the DOI to identify the software release instead of the preprint, connect your GitHub account to Zenodo, enable the repository, then create the GitHub `v1.7.0` release. Zenodo will archive the release as software. The included `CITATION.cff` supplies software citation metadata. Do not use both this route and a manual Zenodo upload for the exact same object unless you intentionally want separate records with different resource identities.

## Before publication

- GitHub and Zenodo identifiers are already inserted in the LaTeX source.
- Decide whether to retain the current custom restrictive license or replace it with standard licenses. The package preserves the license already present in the earlier project; it has not been changed automatically.
- Recompile the paper after inserting the final identifiers.
- Run `pytest -q` one final time.
- Tag the exact final GitHub commit as `v1.7.0`.
