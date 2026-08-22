# Phase B2 v0.4.0: Frozen Architecture and Coefficient-Origin Checkpoint

## Question addressed

Does the successful local flavor operator leave only one missing order-one coefficient, or is a larger matching law still required?

## Correction to the one-coefficient interpretation

The successful operator architecture is now fixed at the functional level:

\[
B_{Ai}(x)=q_{Ai}O(x)+\kappa_{Ai}\,\mathcal I_G(x)+C_{\rm geo}(x),
\]

where `O` is the kink-like wall transition, `I_G` is the actual normalized local gradient density, and `C_geo` is the independently calculated Hessian-channel correction.  This architecture has a local Hermitian action origin and preserves the required chiral spectrum.

The remaining unknown is **not one scalar coefficient**.  In the present Route-I coordinates there are seven effective combinations.  They can be represented by conventional generation and up/down projectors with order-one coefficients, but the law selecting those coefficients has not been derived.

This is better described as a missing UV matching rule or flavor sector than as an inaccurate bounce equation.

## Identifiability

The logarithmic observable Jacobian at the exact local-core solution has seven nonzero singular values:

\[
72.10,\ 32.40,\ 12.63,\ 3.395,\ 1.299,\ 0.324,\ 0.0341.
\]

Thus all seven local directions are identifiable in principle, although the condition number is approximately `2111`.  The weakest direction is mainly a correlated change of `h_u0`, `h_d0`, and `h_Q`.  The data therefore do not reduce the problem to one undetermined number, but they do reveal one weakly constrained combination.

## Exploratory four-fixed/three-fitted hypothesis

A physically motivated, but post-hoc, benchmark hypothesis was frozen:

\[
\epsilon_c = \text{center residual-energy fraction},\qquad
k_{\rm geo}=\text{two-channel Hessian correction},
\]

\[
L_{\rm lock}=\alpha(R_{\rm mix}-R_{\rm grad}),
\]

and

\[
h_Q=3-\frac{2}{5}\epsilon_c,
\]

\[
h_{u0}=h_Q+k_{\rm geo}+\frac14\epsilon_c,
\]

\[
h_{d0}=h_Q-\frac25 k_{\rm geo},
\]

\[
a_{d1}=L_{\rm lock}\left(\frac12-\frac13\epsilon_c\right).
\]

Only `h_u1`, `a_d0`, and `h_d1` were then adjusted.  The maximum error among the seven flavor observables was

\[
0.3234\%.
\]

This meets the numerical compression objective of four fixed controls plus three fitted controls.

A tentative seven-formula completion, with no continuous flavor fit, reaches a maximum error of `2.09%`.  It is close but does not meet the one-percent strong-realization standard.

## Cross-wall result

The four-fixed hypothesis was also evaluated on the 51-wall family while fitting only the remaining three controls.  It works reasonably near the benchmark but is not a universal law across large stiffness displacements:

- 14 of 51 walls remain below 1%;
- 31 of 51 remain below 2%;
- 44 of 51 remain below 3%;
- the large failures occur primarily on the extreme `Y` corridor.

This result prevents the benchmark formulas from being promoted to first-principles laws.  It also does not falsify constant UV flavor coefficients: a physically different wall is expected to predict different flavor observables rather than allowing its coefficients to be re-fitted to our universe.

## Interpretation

What is frozen:

1. the corrected bounce and physical wall;
2. the local gradient/portal core operator;
3. the universal Hessian-channel correction;
4. the chiral sign and normalizability structure;
5. the ordinary up/down representation projectors.

What remains open:

1. the UV matching law for the order-one flavor matrices;
2. whether a discrete/non-Abelian symmetry relates their eigenvalues;
3. whether integrating out a small heavy-fermion mediator sector generates the observed rational-looking relations;
4. a blind prediction after the matching law is frozen.

## Next step

The next calculation should introduce the smallest explicit heavy-fermion/flavon mediator sector and integrate it out.  The target is not another seven-parameter fit.  It is to derive the effective diagonal matrices from:

- one common generation representation;
- up/down representation projectors;
- one or two mediator masses/couplings;
- the already fixed local wall invariants.

Any proposed structure must reduce the independent coefficient count and be frozen before evaluating a new observable.
