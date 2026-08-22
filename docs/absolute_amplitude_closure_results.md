# Phase B2 v1.2.0: Absolute-Amplitude Closure

## Objective

The protected Pati–Salam auxiliary construction fixed the representation ratios

\[
\frac{h_Q}{h_{d0}}=\frac{22}{21},\qquad
\frac{h_{u0}}{h_{d0}}=\frac{23}{21},\qquad
\frac{h_{u1}}{h_{d1}}=\frac1{21},
\]

but left four continuous benchmark amplitudes:

\[
h_{d0},\qquad h_{d1},\qquad a_{d0},\qquad a_{d1}.
\]

In the raw local-gradient basis, the first two are written

\[
h_{d0}=c_{d0}G_{\max},\qquad
h_{d1}=c_{d1}G_{\max}.
\]

This checkpoint asks whether those four amplitudes can be organized by independently available gauge and wall quantities rather than fitted separately.

## Reference four-parameter solution

The protected four-parameter refit gives

| Amplitude | Value |
|---|---:|
| \(c_{d0}\) | 3.54179145 |
| \(c_{d1}\) | -0.87435738 |
| \(a_{d0}\) | 0.26133958 |
| \(a_{d1}\) | 0.24186493 |

It reproduces the seven flavor targets with maximum error

\[
0.59571\%.
\]

## Shared-normalization discovery

Three amplitudes independently define nearly the same attenuation factor:

\[
Z_F^{(c)}=\frac{c_{d1}}{-1/(4g_Lg_R)}=0.99124427,
\]

\[
Z_F^{(a0)}=\frac{a_{d0}}{m_{\Phi,T}/4}=0.99166339,
\]

\[
Z_F^{(a1)}=\frac{a_{d1}}{L_{\rm lock}/2}=0.99143537,
\]

where

\[
L_{\rm lock}=\alpha(R_{\rm mix}-R_{\rm grad})=0.48790861.
\]

Their full relative span is only

\[
0.04227\%.
\]

This is evidence that \(c_{d1},a_{d0},a_{d1}\) do not represent three unrelated normalizations. They appear to share one common seed/wall renormalization.

## Frozen multiplicity-normalized hypothesis

The following leading-order map was frozen:

\[
\overline{g_W^{-2}}
=\frac12\left(g_L^{-2}+g_R^{-2}\right),
\]

\[
Z_0=1+\frac{\epsilon_c}{21},
\]

\[
Z_F=\exp\left(-\frac{2\epsilon_c}{15}\right),
\]

\[
\boxed{c_{d0}=\overline{g_W^{-2}}Z_0},
\]

\[
\boxed{c_{d1}=-\frac{Z_F}{4g_Lg_R}},
\]

\[
\boxed{a_{d0}=Z_F\frac{m_{\Phi,T}}4},
\qquad
a_{d1}=Z_F\frac{L_{\rm lock}}2}.
\]

The interpretation of the integer factors is:

- 21: the complete Pati–Salam adjoint multiplicity already derived in the protected seed kernel;
- 15: the number of \(SU(4)_C\) adjoint directions participating in the quark-carrying seed block;
- 2: the two fermion legs of the bilinear response;
- 4: the \(SU(4)_C\) fundamental multiplicity;
- 2 in \(L_{\rm lock}/2\): the two-channel locking interval.

The numerical inputs are

\[
g_L=0.54185306,\qquad g_R=0.52305841,
\]

\[
\epsilon_c=0.06258787,\qquad
m_{\Phi,T}=1.05414633.
\]

These give

| Amplitude | Formula | Prediction | Difference from four-parameter refit |
|---|---|---:|---:|
| \(c_{d0}\) | \(\frac12(g_L^{-2}+g_R^{-2})(1+\epsilon_c/21)\) | 3.54104355 | -0.02112% |
| \(c_{d1}\) | \(-Z_F/(4g_Lg_R)\) | -0.87475026 | +0.04493% |
| \(a_{d0}\) | \(Z_Fm_{\Phi,T}/4\) | 0.26134651 | +0.00265% |
| \(a_{d1}\) | \(Z_FL_{\rm lock}/2\) | 0.24192697 | +0.02565% |

All four amplitudes are therefore reproduced within 0.05% at the benchmark.

## Zero-continuous-fit flavor result

After applying the protected 21 relations, the resulting seven controls are

\[
(2.81149921,\ 2.93929463,\ -0.03156953,\ 0.26134651,\ 0.24192697,\ 2.68370379,\ -0.66296010).
\]

No continuous flavor amplitude is optimized after inserting the formulas.

| Observable | Relative error |
|---|---:|
| \(m_c/m_t\) | -0.53864% |
| \(m_u/m_t\) | +0.07920% |
| \(m_s/m_b\) | +0.36216% |
| \(m_d/m_b\) | -0.25271% |
| \(|V_{us}|\) | +0.27159% |
| \(|V_{cb}|\) | -0.10969% |
| \(|V_{ub}|\) | +0.63290% |

Thus

\[
\boxed{\max |\Delta O/O|=0.63290\%}
\]

with RMS error

\[
0.37374\%.
\]

Allowing only the common \(Z_F\) to be adjusted improves the maximum error slightly to

\[
0.60070\%,
\]

at

\[
Z_F=0.99135924.
\]

The multiplicity formula predicts \(Z_F=0.99168967\), so the zero-fit point is already effectively at the one-parameter optimum.

## Chiral-spectrum audit

The zero-fit amplitude map retains:

- 9/9 desired profiles with exactly one near-zero chiral state;
- 9/9 opposite-chirality profiles with no near-zero state;
- minimum opposite-sector eigenvalue 0.876875.

The result is not obtained by weakening localization or introducing additional light states.

## Cross-wall audit

The formulas were applied to all 51 available Phase-A walls and compared with independently refitted four-amplitude values.

For the 33-point local design, mean absolute amplitude errors are

| Amplitude | Mean absolute error | Worst error |
|---|---:|---:|
| \(c_{d0}\) | 0.475% | 1.357% |
| \(c_{d1}\) | 0.745% | 2.027% |
| \(a_{d0}\) | 1.103% | 2.768% |
| \(a_{d1}\) | 1.313% | 2.843% |

Over the entire 51-wall set, including the extreme \(\pm8\%\) stiffness corridor, the worst errors reach approximately 3%, 5%, 11%, and 12%, respectively.

The central response audit also shows that the candidate formulas do not reproduce all fitted \(X,Y,Z\) derivatives. This prevents them from being promoted to universal wall laws.

The likely missing ingredients are wall-dependent running of the inherited Pati–Salam gauge normalizations and higher-order seed/wall response terms.

## What has been achieved

1. The protected 21 kernel remains exact.
2. The four benchmark amplitudes are reproduced within 0.05% by a zero-continuous-fit multiplicity map.
3. Three amplitudes share one independently visible normalization at the 0.042% level.
4. The seven flavor observables remain within 0.633%.
5. The complete desired chiral spectrum survives.

## Claim boundary

The representation ratios \(22/21,23/21,1/21\) have an explicit protected auxiliary origin.

The new absolute-amplitude formulas do **not yet** have the same status. They were identified after examining the benchmark and use inherited high-scale gauge couplings. The factors \(1+\epsilon_c/21\) and \(\exp(-2\epsilon_c/15)\) are physically organized leading-order matching hypotheses, not yet coefficients obtained from a complete seed-sector loop or bound-state calculation.

The appropriate classification is:

> **Zero-fit benchmark closure and a strong candidate amplitude law, not yet an independent first-principles prediction.**

## Next calculation

The next calculation should derive the two normalization factors from the explicit seed action:

\[
Z_0=1+\epsilon_c/21,
\qquad
Z_F=e^{-2\epsilon_c/15},
\]

and include the wall dependence of \(g_L,g_R\). The target is to reproduce both the benchmark values and the local \(X,Y,Z\) response matrix without using flavor data.
