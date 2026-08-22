# Phase B2 v1.5.0 — Explicit Radial UV Operator Matching

## Purpose

The v1.4.0 checkpoint showed that the corrected bounce supplies a pure bulk radial incomplete-settling mode and a bounce-only finite-amplitude number

\[
\lambda_{\rm radial}=1.071855028,
\]

which, when inserted into the frozen flavor construction, gives a maximum quark-observable error of \(0.6111\%\).  The remaining question was whether an explicit local seed action **forces** this number or merely permits it.

This checkpoint writes the most general local radial dependence of the protected seed and fermion-source sectors and integrates it out exactly.

## 1. General local UV action

Let \(\mathcal R\) denote any canonically defined dimensionless radial invariant.  The gauge- and \(O(22)\)-compatible seed sector may contain

\[
\mathcal L_{\rm seed}\supset
-\frac12 M_X^2(\mathcal R)X^IX^I
-\mu_X(\mathcal R)S\,X^IX^I,
\]

while the fermion source can contain

\[
\mathcal L_{\rm source}\supset
-y_F(\mathcal R)S\,\bar\Psi\Psi.
\]

The protected auxiliary response gives

\[
Z_0(\mathcal R)
=
\frac{\mu_X(\mathcal R)/M_X^2(\mathcal R)}
     {\mu_X(0)/M_X^2(0)},
\]

and the two fermion-source vertices give

\[
Z_F(\mathcal R)
=
\left[\frac{y_F(\mathcal R)}{y_F(0)}\right]^2.
\]

Expanding around the true vacuum,

\[
\frac{d\ln Z_0}{d\mathcal R}
=
\frac{d\ln\mu_X}{d\mathcal R}
-
\frac{d\ln M_X^2}{d\mathcal R},
\]

\[
\frac12\frac{d\ln Z_F}{d\mathcal R}
=
\frac{d\ln y_F}{d\mathcal R}.
\]

These derivatives are independent Wilson coefficients in the general local action.

### Central conclusion

Canonical normalization fixes the normalization of the radial field or radial invariant.  It does **not** fix how an unrelated seed mass, seed source, or fermion source depends on that radial variable.

Thus neither canonical normalization, \(O(22)\), nor Pati–Salam alone forces

\[
\lambda_{\rm res}=\lambda_{\rm radial}.
\]

\(O(22)\) still enforces equality among the 22 seed channels, and Pati–Salam still fixes the multiplicity factors 21 and 15.  They do not fix the absolute radial derivative of the seed/source functions.

## 2. The lowest-dimension radial invariant

Because the bulk SFV field is complex, the leading local phase-invariant radial departure is

\[
\chi_\Phi
=
1-\frac{|\Phi|^2}{\rho^2}.
\]

At the bubble center,

\[
\chi_c=0.250175660.
\]

The normalized residual energy is

\[
\epsilon_c
=
\left(1-\frac{|\Phi_c|^2}{\rho^2}\right)^2
=\chi_c^2
=0.062587861.
\]

Therefore a lowest-dimension portal such as \(|\Phi|^2X^2\) responds linearly to \(\chi\sim\sqrt\epsilon\), not to \(\epsilon\).

This distinction matters numerically.  Feeding the different bounce-derived radial drivers into the frozen flavor construction gives:

| Radial driver \(D\) in the exponent | Worst observable error |
|---|---:|
| exact residual energy \(D=\epsilon_c\) | 0.6321% |
| finite-amplitude radial prescription \(D=\lambda_{\rm radial}\epsilon_c\) | **0.6111%** |
| harmonic energy fraction \(D=4u^2\) | 0.6884% |
| canonical displacement \(D=u\) | 2.6078% |
| unit lowest-dimension portal \(D=\chi_c\) | 7.3177% |

A unit-strength renormalizable radial portal is therefore incompatible with the successful mild residual dressing.  A smaller portal coefficient or a compensator relation is needed.

## 3. What the previous radial matching means in the lowest-dimension basis

Let

\[
u=\frac{\rho-\Phi_c}{\rho}=0.134076019.
\]

The exact identities are

\[
\chi_c=u(2-u),
\qquad
\epsilon_c=\chi_c^2,
\qquad
\lambda_{\rm radial}=\frac{1}{1-u/2}.
\]

Therefore

\[
\frac{\lambda_{\rm radial}\epsilon_c}{\chi_c}
=2u
=0.268152039.
\]

So the v1.4.0 radial prescription is equivalent to a lowest-dimension \(\chi\)-portal with a specific common weight

\[
w_\chi=2u.
\]

That weight is bounce-derived, but the explicit local action does not force the seed and fermion-source functions to carry it.  Doing so is a **radial compensator matching condition**.

The common coefficient preferred diagnostically by the seven flavor observables corresponds to

\[
w_\chi\simeq0.2691,
\]

only about \(0.4\%\) above the bounce value \(2u\).  This explains the excellent numerical result but does not remove the UV matching assumption.

## 4. Renormalizable rational-portal diagnostic

For illustration, take simple rational responses

\[
Z_0=\frac{1}{1+\eta_0\chi},
\qquad
Z_F=\frac{1}{(1+\eta_F\chi)^2}.
\]

To reproduce the bounce radial matching at the benchmark requires approximately

\[
\eta_0=-0.01275,
\qquad
\eta_F=+0.01792.
\]

These are small, technically plausible coefficients, but they are independent and have opposite signs.  A single universal portal coefficient does not generate both factors.

This is a useful no-go statement:

> The minimal action requires either distinct radial charges for the seed ratio and fermion source, or one compensator symmetry that enforces the required opposite weights.

## 5. What remains valid

The negative uniqueness result does not undo the earlier successes:

1. The local wall-core operator and chiral spectrum remain valid.
2. The rank-two flavor structure remains valid.
3. The protected Pati–Salam matrix and the integer 21 remain valid.
4. The rational relations \(22/21\), \(23/21\), and \(1/21\) remain valid within the protected construction.
5. The bounce-derived canonical-radial matching still gives a zero-continuous-fit result at \(0.6111\%\).

What changes is the claim boundary: that result is not forced by the already stated symmetries.

## 6. Scientific classification

The explicit calculation establishes

\[
\boxed{
\text{one genuine radial portal/compensator coefficient remains.}
}
\]

The model has therefore reached a highly compressed, physically allowed flavor construction, but not a unique parameter-free UV prediction.

The correct next step is to search for a symmetry that makes the radial mode a compensator and fixes the relative weights

\[
q_0=+\frac1{21},
\qquad
q_F=-\frac1{15},
\]

with common strength \(w_\chi=2u\), or to derive those weights from a geometric modulus/dilaton-like origin.  Without such an added principle, the final coefficient is an allowed Wilson input.

## Reproducibility

```bash
python src/explicit_radial_uv_operator_matching.py
pytest -q tests/test_explicit_radial_uv_matching.py
```
