# v1.4.0

- Added absolute-amplitude closure for the protected 21 mediator kernel.
- Identified one shared attenuation across `c_d1`, `a_d0`, and `a_d1` with only 0.042% relative span.
- Added the frozen multiplicity-normalized formulas:
  - `c_d0 = 0.5*(g_L^-2+g_R^-2)*(1+epsilon_c/21)`
  - `Z_F = exp(-2 epsilon_c/15)`
  - `c_d1 = -Z_F/(4 g_L g_R)`
  - `a_d0 = Z_F*m_true_Phi/4`
  - `a_d1 = Z_F*L_lock/2`
- Reproduced all four benchmark amplitudes within 0.05%.
- Obtained a zero-continuous-fit seven-observable result with 0.6329% maximum error.
- Preserved all 9 desired chiral zero modes and excluded all 9 opposite-chirality zero modes.
- Added 51-wall amplitude tracking and a central response-matrix audit.
- Added seven automated tests.

## v1.4.0
- Derived a local residual-relaxation spurion form for `Z0` and `ZF`.
- Reduced the four absolute amplitudes to one common order-one coefficient `lambda_res`.
- Unit normalization gives a zero-fit maximum flavor error of 0.6321%.
- Best common coefficient is 1.07573; best two-coefficient diagnostic is 0.5943%.
- Proved the claim boundary: Pati-Salam/O(22) multiplicities fix normalized projector dimensions, not the common Wilson coefficient.
- Reclassified the 51-wall amplitude scan as a compensation-fit diagnostic rather than an independent UV-coefficient measurement.

## v1.4.0
- Extracted the incomplete-settling mode directly from the corrected two-field bounce.
- Verified that the center residual is a pure bulk soft/radial mode and that
  `epsilon_c = (1-(Phi_c/rho)^2)^2` across all 51 walls.
- Derived the finite-amplitude radial conversion `lambda_radial = 1/(1-u/2) = 1.071855028` without flavor data.
- Obtained a zero-continuous-fit seven-observable result with 0.6111% maximum error.
- Preserved the claim boundary: identifying this canonical radial factor with the seed Wilson normalization remains a UV matching principle, not a symmetry theorem.

## v1.5.0

- Added explicit local radial dependence for seed masses, seed sources, and fermion sources.
- Proved that canonical radial normalization does not fix the independent Wilson derivatives entering `Z0` and `ZF`.
- Identified the leading phase-invariant renormalizable portal variable `chi = 1-|Phi|^2/rho^2 = sqrt(epsilon)`.
- Showed that a unit-strength `chi` portal overcorrects the frozen flavor model (7.32% worst error).
- Reclassified the 0.611% radial result as a minimal compensator matching hypothesis rather than a uniqueness theorem.

## v1.6.0

- Solved the coupled O(4) fluctuation operator and verified one negative breathing mode plus the near-zero `ell=1` translation mode.
- Measured a substantial but non-unit overlap between the incomplete-settling deformation and the geometric breathing mode.
- Derived the induced-measure/canonical-normalization Weyl weights for auxiliary and propagating mediator realizations in three and four worldvolume dimensions.
- Proved that ordinary embedding geometry does not generate the required opposite `+1/21` and `-2/15` compensator charges.
- Audited the numerical `S^3` shell-volume near-match and showed that it is not a universal 51-wall law.
- Introduced the minimal internal unimodular `O(21) x O(15)` block-volume symmetry, which reproduces `Z0` and `ZF` exactly.
- Isolated the remaining structural condition: the internal response metric must be induced by the canonical radial bounce mode.
- Added five regression tests for the geometric-modulus checkpoint.

## v1.7.0

- Replaced the informal block-volume ansatz with an explicit constrained auxiliary response-metric action.
- Proved the uniqueness of `T=diag(I21/21,-I15/15)` under block isotropy, determinant one, and the convention `Sigma=ln det C21`.
- Derived the exact action minimum `A=Sigma_rad`, `B=-Sigma_rad` independently of the positive constraint stiffnesses.
- Expressed the bounce driver in closed form as `Sigma_rad=2(1-Phi_c/rho)^2(1+Phi_c/rho)`.
- Preserved the zero-continuous-fit `0.6111%` seven-observable result.
- Demonstrated that the preferred auxiliary realization adds no spacetime dimension and no propagating scalar.
- Quantified a propagating-modulus alternative: sub-1% flavor requires `m_Sigma/p >= 1.55`.
- Verified one-loop gauge stability of the determinant constraint and normalized radial charges.
- Clarified the claim boundary: the compensator locking is a consistent minimal extension, but remains one added structural postulate beyond the original two-scalar action.
