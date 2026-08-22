# Phase B2 v1.7.0 — Constrained Internal-Response-Metric Action

## Purpose

The preceding checkpoint found that an internal response metric with two isotropic blocks,

\[
C_{21}=e^{\Sigma/21}\mathbf 1_{21},\qquad
C_{15}=e^{-\Sigma/15}\mathbf 1_{15},
\]

reproduces the compensator factors required by the frozen flavor construction.  This checkpoint asks whether that structure can be written as a precise action, whether it adds physical dimensions or particles, and whether its determinant-one relation survives quantum running.

The answer is:

1. a minimal algebraic action exists;
2. its block generator is unique once the block dimensions, isotropy, and determinant constraint are declared;
3. it introduces no new spacetime dimensions and need not introduce a new propagating field;
4. its determinant relation and normalized radial charges survive one-loop gauge running;
5. the locking of the response metric to the bounce radial driver is an explicit additional compensator principle, not a theorem of the original two-scalar SFV/dSB action.

## 1. The smallest symmetry is one rescaling, not an `SL(36)` world

Let the logarithmic eigenvalue of the 21-component seed-response block be \(a\), and that of the 15-component fermion-source block be \(b\):

\[
C_{21}=e^a\mathbf 1_{21},\qquad C_{15}=e^b\mathbf 1_{15}.
\]

The determinant-one condition is

\[
21a+15b=0.
\]

Define the common logarithmic volume coordinate by

\[
\Sigma\equiv\ln\det C_{21}=21a.
\]

Then the solution is unique:

\[
\boxed{a=\frac{\Sigma}{21},\qquad b=-\frac{\Sigma}{15}.}
\]

Equivalently, the generator in the 36-component bookkeeping representation is

\[
T=\operatorname{diag}\left(\frac1{21}\mathbf 1_{21},-\frac1{15}\mathbf 1_{15}\right),
\]

with

\[
\operatorname{Tr}T=0,\qquad
\operatorname{Tr}T^2=\frac1{21}+\frac1{15}=\frac4{35}.
\]

No physical symmetry is required to rotate the 21 block into the 15 block.  Such a full rotation would conflict with their different gauge roles.  The actual symmetry is only a one-parameter determinant-preserving internal rescaling that acts uniformly within each block and commutes with Pati–Salam.

Thus the number 36 is only the total number of response components used to display the traceless generator.  It is not a spacetime dimension and not a compact internal coordinate space.

## 2. Explicit auxiliary action

Write

\[
A\equiv\ln\det C_{21}=21a,
\qquad
B\equiv\ln\det C_{15}=15b.
\]

The minimal invariant auxiliary potential is

\[
\boxed{
V_{\rm resp}
=
\frac{\Lambda_U^4}{2}(A+B)^2
+
\frac{\Lambda_L^4}{2}
\left[\frac{A-B}{2}-\Sigma_{\rm rad}\right]^2
+V_{\rm shape}.
}
\]

Here:

- the first term enforces the determinant-one condition;
- the second locks the relative block response to the existing bounce radial driver;
- \(V_{\rm shape}\) is a positive, block-isotropic penalty that removes traceless distortions within each block.

For any positive \(\Lambda_U^4\) and \(\Lambda_L^4\), minimization gives

\[
A+B=0,
\qquad
\frac{A-B}{2}=\Sigma_{\rm rad},
\]

and therefore

\[
\boxed{A=\Sigma_{\rm rad},\qquad B=-\Sigma_{\rm rad}.}
\]

The mass scales control how stiff the constraints are, but they do not change the minimum.

When \(C_{21}\) and \(C_{15}\) are algebraic response metrics, they carry no kinetic terms.  They can be eliminated exactly and add no propagating particle.

## 3. The bounce radial source

The earlier bounce calculation gives

\[
x_c\equiv\frac{\Phi_c}{\rho}=0.8659239806,
\qquad
u\equiv1-x_c=0.1340760194.
\]

The exact residual-energy identity is

\[
\epsilon_c=(1-x_c^2)^2=0.0625878607,
\]

and the finite-amplitude radial conversion is

\[
\lambda_{\rm radial}=\frac1{1-u/2}=1.0718550277.
\]

Their product has the closed form

\[
\boxed{
\Sigma_{\rm rad}
=
\lambda_{\rm radial}\epsilon_c
=
2(1-x_c)^2(1+x_c)
=
0.0670851132.
}
\]

This identity is bounce-derived and uses no flavor observable.

The constrained action then gives

\[
Z_0=e^{\Sigma_{\rm rad}/21}=1.0031996371,
\]

\[
Z_F=e^{-2\Sigma_{\rm rad}/15}=0.9910952029.
\]

## 4. Frozen flavor consequence

Using the protected mediator ratios and the previously derived amplitude map gives the seven observable errors

| Observable | Error |
|---|---:|
| \(m_c/m_t\) | \(-0.6111\%\) |
| \(m_u/m_t\) | \(-0.0305\%\) |
| \(m_s/m_b\) | \(+0.1720\%\) |
| \(m_d/m_b\) | \(-0.2491\%\) |
| \(|V_{us}|\) | \(+0.3545\%\) |
| \(|V_{cb}|\) | \(-0.1418\%\) |
| \(|V_{ub}|\) | \(+0.5341\%\) |

The worst error remains

\[
\boxed{0.6111\%}.
\]

No continuous flavor amplitude is adjusted in this result.

## 5. Is the locking fragile?

The flavor calculation is more tolerant than the near-exact benchmark agreement might suggest.  Scaling the entire radial response by a factor \(r\),

\[
\Sigma\rightarrow r\Sigma_{\rm rad},
\]

keeps every tested observable below 1% over

\[
\boxed{0.706\lesssim r\lesssim1.361.}
\]

Changing only the seed response allows approximately

\[
0.398\lesssim r_0\lesssim1.495,
\]

while changing only the fermion response allows

\[
0.463\lesssim r_F\lesssim1.438.
\]

Thus ordinary loop-sized corrections do not immediately endanger the sub-percent result.

## 6. A propagating modulus is optional and constrained

If the relative block mode has a kinetic term and mass \(m_\Sigma\), its finite-momentum response is approximately

\[
r(p)=\frac{m_\Sigma^2}{m_\Sigma^2+p^2}.
\]

The 1% flavor requirement gives

\[
\boxed{m_\Sigma/p\gtrsim1.55.}
\]

Using the two inherited wall-width scales gives:

| Wall scale | Minimum dimensionless \(m_\Sigma\) for sub-1% flavor |
|---|---:|
| action FWHM | 0.885 |
| gradient FWHM | 0.596 |

A sufficiently heavy propagating mode is viable.  The simpler auxiliary realization avoids this condition completely and creates no fifth-force or light-scalar concern.

## 7. One-loop stability

### Determinant constraint

The determinant condition is algebraic:

\[
\ln\det C_{21}+\ln\det C_{15}=0.
\]

It is exact in the constrained path integral and cannot drift through ordinary multiplicative wave-function renormalization.

### Seed block

For every Pati–Salam adjoint block, the previous beta-function calculation found

\[
\frac{d}{d\ln Q}\ln\frac{\mu_a}{M_a^2}=0
\]

at one-loop gauge order.  Therefore all 21 directions retain the same protected zero-momentum response despite their different gauge running.

### Fermion-source block

Gauge contributions to a source vertex are multiplicative and independent of \(\Sigma\).  Consequently the normalized ratio

\[
\frac{y_F(\Sigma,Q)}{y_F(0,Q)}
\]

retains its exponential radial charge under gauge running.

A generic nonlinear source term, modeled by

\[
\frac{dy}{d\ln Q}=\frac{a y^3}{16\pi^2},
\]

can renormalize the effective radial exponent.  The numerical audit shows that the direct fermion-charge result remains under 1% until the accumulated nonlinear parameter

\[
\zeta=\frac{2a y^2\ln(\Lambda/Q)}{16\pi^2}
\]

reaches approximately 1.165.  This is already a strong, not ordinary small-loop, correction.

Therefore the determinant relation is exact and the radial weights are comfortably stable in the perturbative regime.

## 8. Precise claim boundary

This checkpoint establishes a fully explicit, internally consistent extension:

\[
\boxed{
\text{SFV/dSB bounce}
+
\text{protected Pati–Salam seed sector}
+
\text{one determinant-preserving response compensator}.
}
\]

Within that augmented action:

- no extra spacetime dimension is added;
- no new propagating field is required;
- the generator is unique;
- the block determinants are exactly locked;
- the weights \(+1/21\) and \(-1/15\) are forced;
- the two fermion legs force \(-2/15\);
- the radial source is fixed by the bounce;
- the zero-fit flavor prediction remains at 0.611%.

The remaining qualification is conceptual rather than numerical:

> The locking invariant is one new structural principle.  It is not derivable from the original two-scalar SFV/dSB Lagrangian without adding this response-compensator relation.

Thus this is a consistent minimal completion, not proof that the unextended original action uniquely demanded it.

## Reproduction

```bash
python src/constrained_internal_response_metric.py
pytest -q tests/test_constrained_internal_response_metric.py
```
