# Phase B2 v1.3.0 — Residual-spurion derivation of \(Z_0\) and \(Z_F\)

## Objective

Phase B2 v1.2.0 found the benchmark normalization factors

\[
Z_0\simeq 1+\frac{\epsilon_c}{21},
\qquad
Z_F\simeq \exp\!\left(-\frac{2\epsilon_c}{15}\right),
\]

where \(\epsilon_c=0.0625878658\) is the fractional excess energy of the bounce center above the stationary true vacuum.  Those formulas reproduced the four remaining amplitudes within 0.05%, but were identified after inspecting the benchmark.

This checkpoint asks whether the explicit protected seed action actually derives those factors and whether multiplicity alone fixes their numerical coefficients.

## 1. Residual-relaxation spurion

Introduce a real dimensionless scalar spurion \(\mathcal R\) representing the incomplete post-tunneling relaxation of the bounce center, normalized at the benchmark by

\[
\langle\mathcal R\rangle=\epsilon_c.
\]

The protected seed sector contains the complete Pati–Salam adjoint multiplicity

\[
d_{\rm PS}=15+3+3=21,
\]

while the quark-carrying color block has

\[
d_4=\dim\operatorname{Adj}(SU(4)_C)=15.
\]

The lowest-order linked response compatible with the previously derived mediator structure can dress:

1. the Pati–Salam-singlet collective source, normalized by the complete adjoint trace;
2. the two fermion-source vertices, normalized by the \(SU(4)_C\) adjoint trace.

Writing the corresponding dimensionless Wilson coefficients as \(\lambda_0\) and \(\lambda_F\), repeated local insertions exponentiate to

\[
\boxed{Z_0=\exp\!\left(\lambda_0\frac{\epsilon_c}{21}\right)},
\]

\[
\boxed{Z_F=\exp\!\left(-2\lambda_F\frac{\epsilon_c}{15}\right)}.
\]

The factor of two in \(Z_F\) counts the two source vertices in the bilinear wall response.  In the common \(O(22)\)-normalized limit,

\[
\lambda_0=\lambda_F\equiv\lambda_{\rm res}.
\]

This yields one residual-response coefficient for all four absolute amplitudes.

## 2. Resulting amplitude map

The protected Pati–Salam relations remain

\[
\frac{h_Q}{h_{d0}}=\frac{22}{21},\qquad
\frac{h_{u0}}{h_{d0}}=\frac{23}{21},\qquad
\frac{h_{u1}}{h_{d1}}=\frac1{21}.
\]

The absolute amplitudes become

\[
c_{d0}=\frac12\left(g_L^{-2}+g_R^{-2}\right)Z_0,
\]

\[
c_{d1}=-\frac{Z_F}{4g_Lg_R},
\]

\[
a_{d0}=Z_F\frac{m_{\Phi,T}}4,
\qquad
a_{d1}=Z_F\frac{L_{\rm lock}}2.
\]

Thus the previous four continuous amplitudes reduce to one common residual-response coefficient.

## 3. Unit-normalized prediction

Taking the canonical unit choice

\[
\lambda_{\rm res}=1
\]

gives

\[
Z_0=1.0029848203,
\qquad
Z_F=0.9916896745.
\]

No continuous flavor quantity is then optimized.  The seven observable errors are:

| Observable | Relative error |
|---|---:|
| \(m_c/m_t\) | -0.53891% |
| \(m_u/m_t\) | +0.07940% |
| \(m_s/m_b\) | +0.36025% |
| \(m_d/m_b\) | -0.25121% |
| \(|V_{us}|\) | +0.27339% |
| \(|V_{cb}|\) | -0.10871% |
| \(|V_{ub}|\) | +0.63207% |

Therefore

\[
\boxed{\max|\Delta O/O|=0.63207\%}
\]

with zero continuous flavor fitting after the action and unit normalization are specified.

## 4. One- and two-coefficient diagnostics

Allowing one common coefficient to adjust gives

\[
\lambda_{\rm res}=1.0757296,
\]

with

\[
\max|\Delta O/O|=0.61497\%.
\]

Allowing the even and fermionic responses to have separate coefficients gives

\[
\lambda_0=1.1424280,
\qquad
\lambda_F=1.0366577,
\]

and

\[
\max|\Delta O/O|=0.59434\%.
\]

The independently inferred benchmark values from the four-amplitude reference are

\[
\lambda_0\simeq1.0694,
\qquad
\lambda_F\simeq1.0292.
\]

The model is therefore not hiding large or unnatural constants.  The remaining normalization is one order-one coefficient close to unity.

## 5. Why multiplicity does not finish the derivation

The dimensions 21 and 15 determine natural normalized projectors, for example

\[
\frac1{21}\operatorname{Tr}_{\rm Adj(PS)}X^2,
\qquad
\frac1{15}\operatorname{Tr}_{\rm Adj(SU4)}X_4^2.
\]

But the action may contain

\[
\lambda\,\mathcal R\frac1d\operatorname{Tr}X^2
\]

with an arbitrary dimensionless Wilson coefficient \(\lambda\).  Rescaling the definition of the normalized invariant shifts the coefficient in the opposite direction.  Consequently, group multiplicity fixes the **projector and channel count**, but not the absolute residual-to-seed response.

The allowed residual-spurion operator ledger contains independent coefficients for

\[
\mathcal R X_0^2,
\quad
\mathcal R\operatorname{Tr}X_4^2,
\quad
\mathcal R\operatorname{Tr}X_L^2,
\quad
\mathcal R\operatorname{Tr}X_R^2,
\]

and for wall- and fermion-source dressing.  An exact \(O(22)\) boundary can equate some of these coefficients, reducing them to one, but does not set their common value to one.

Therefore, declaring \(\lambda_{\rm res}=1\) is a natural canonical matching condition, not yet a theorem derived from the current action.

## 6. Resummation ambiguity is numerically small

Several physically reasonable resumptions were compared at unit coefficient:

| Scheme | \(Z_0\) | \(Z_F\) | Worst flavor error |
|---|---:|---:|---:|
| linear \(Z_0\), exponential \(Z_F\) | 1.0029804 | 0.9916897 | 0.63290% |
| exponential both | 1.0029848 | 0.9916897 | 0.63207% |
| propagator resummation | 1.0029893 | 0.9917069 | 0.63291% |
| determinant/geometric mean | 1.0028950 | 0.9919384 | 0.67308% |

All remain below 1%.  The success is not dependent on a finely chosen resummation convention.

## 7. Correct interpretation of the 51-wall audit

The 51-wall amplitude tables were generated by changing the wall and then refitting the amplitudes so that every altered wall continued to reproduce the same observed flavor targets.  Those are **compensation fits**, not independent measurements of UV Wilson coefficients.

A fixed UV action is not required to follow those fitted amplitudes: a universe with different \(X,Y,Z\) should generally predict different flavor observables.  Therefore the failure of one \(\epsilon_c\)-only law to track all compensation fits does not falsify the residual-spurion action.

The cross-wall exercise remains diagnostically useful.  It shows that the compensating effective coordinates are highly basis- and wall-dependent and cannot be used to claim independent validation of \(Z_0\) and \(Z_F\).

## 8. Wall-dependent gauge normalization

The inherited Pati–Salam values \(g_L,g_R\) provide the baseline gauge normalization.  The earlier gauge-normalization module used calibrated group weights and a chosen local kernel, so it does not yet provide an independent \(X,Y,Z\)-dependent prediction for \(g_L,g_R\).

An exploratory comparison found strong correlations between some global bubble-radius quantities and the compensation-fit \(c_{d0}\), but these were inspected post hoc and are not expected to determine a local four-dimensional gauge kinetic coefficient.  They were therefore not promoted to physics laws.

A genuine wall-dependent gauge prediction requires a frozen local gauge kernel and its zero-mode normalization to be recomputed on each wall.

## 9. Scientific status

### Demonstrated

- One local residual-relaxation spurion can generate both \(Z_0\) and \(Z_F\).
- The four amplitude problem reduces to one common natural coefficient.
- The unit-normalized action gives a zero-continuous-fit result within 0.633%.
- The best common coefficient is only 7.6% above unity.
- The conclusion is stable under several resummation choices.

### Not yet demonstrated

- The current symmetries do not force \(\lambda_{\rm res}=1\).
- The mapping from the bounce energy fraction \(\epsilon_c\) to the seed response has not been computed from a specified microscopic interaction.
- Wall-dependent \(g_L,g_R\) are not independently predicted by the existing gauge-normalization module.

## Verdict

\[
\boxed{
\text{The four amplitudes collapse to one natural residual-response coefficient,}
}
\]

but

\[
\boxed{
\text{a fully parameter-free first-principles flavor prediction has not yet been proved.}
}
\]

The strongest current statement is:

> With the protected Pati–Salam seed kernel and the canonically unit-normalized residual coupling, the SFV/dSB bounce predicts the seven real quark-flavor observables at the sub-percent level without continuous flavor fitting.  The remaining theoretical task is to derive, rather than normalize, the residual-seed Wilson coefficient.

## Next calculation

The next controlled calculation should choose one explicit microscopic origin for \(\mathcal R\)—most naturally the radial incomplete-settling mode of the two-field bounce—and calculate its coupling to the protected seed susceptibility.  That computation will determine \(\lambda_{\rm res}\) and decide whether the unit-normalized closure is a prediction or an additional boundary condition.
