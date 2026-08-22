# Phase B2 v1.4.0 — Canonical Radial Mode and Residual Seed Matching

## Purpose

The previous checkpoint reduced the four absolute flavor amplitudes to one common residual-response coefficient,

\[
Z_0=\exp\!\left(\lambda_{\rm res}\epsilon_c/21\right),\qquad
Z_F=\exp\!\left(-2\lambda_{\rm res}\epsilon_c/15\right).
\]

This checkpoint extracts the incomplete-settling mode directly from the corrected two-field bounce and asks whether a bounce-only canonical normalization determines \(\lambda_{\rm res}\).

No quark observable is used to calculate the radial matching coefficient.

## 1. The center residual is a single bulk radial mode

The stationary true vacuum and numerical bounce center are

\[
(\Phi_T,\phi_T)=(2.357142857,0),
\]

\[
(\Phi_c,\phi_c)=(2.041106526,-2.73623\times10^{-7}).
\]

The canonical field-space displacement is

\[
\sigma_c=\sqrt{(\Phi_c-\rho)^2+\phi_c^2}=0.316036331.
\]

Projection onto the true-vacuum Hessian eigenbasis gives a soft-mode norm fraction

\[
P_{\rm soft}=0.999999999999250,
\]

and projection onto the local center Hessian gives

\[
P_{\rm soft,center}=0.999999999999655.
\]

Thus the incomplete settling is not a significant mixture of the bulk and brane modes. It is, to numerical accuracy, the true-bulk radial amplitude mode.

The center residual-energy decomposition is

| Component | Fraction of residual |
|---|---:|
| bulk quartic | 0.999999999985 |
| portal | \(1.49\times10^{-11}\) |
| brane terms | negligible |

## 2. Exact residual identity

Along the center trajectory \(\phi_c\simeq0\), the bulk potential is

\[
V_\Phi(\Phi)=\frac{\lambda_\Phi}{4}(\Phi^2-\rho^2)^2.
\]

Using the false–true gap \(\Delta V=\lambda_\Phi\rho^4/4\) at zero bias gives the exact normalized radial-energy operator

\[
\mathcal R_\Phi(\Phi)
=\frac{V_\Phi(\Phi)-V_\Phi(\rho)}{\Delta V}
=\left(1-\frac{\Phi^2}{\rho^2}\right)^2.
\]

At the bounce center,

\[
\epsilon_c=\mathcal R_\Phi(\Phi_c)=0.0625878626.
\]

This identity was checked on all 51 wall solutions. The maximum absolute discrepancy is

\[
2.51\times10^{-9}.
\]

Therefore the residual spurion used in the amplitude formulas is not an independently invented scalar: it is the normalized bulk radial potential already contained in the bounce action.

## 3. Finite-amplitude canonical normalization

Define the fractional canonical displacement

\[
u=\frac{\rho-\Phi_c}{\rho}=0.134076019.
\]

The true-vacuum harmonic approximation gives

\[
V_{\rm harm}=\frac12m_{\Phi,T}^2(\rho-\Phi_c)^2,
\]

while the exact quartic potential gives

\[
V_{\rm exact}
=V_{\rm harm}\left(1-\frac{u}{2}\right)^2.
\]

Numerically,

\[
V_{\rm harm}=0.0554939747,
\qquad
V_{\rm exact}=0.0483029587,
\]

so

\[
\frac{V_{\rm exact}}{V_{\rm harm}}=0.870418075.
\]

The conversion from the harmonic energy-normalized radial amplitude to the exact finite-displacement canonical amplitude is therefore

\[
\boxed{
\lambda_{\rm radial}
=\sqrt{\frac{V_{\rm harm}}{V_{\rm exact}}}
=\frac{1}{1-u/2}
=1.071855028.
}
\]

This number is calculated entirely from the bounce center and the original quartic potential.

## 4. Frozen seed matching

Identifying the common residual seed normalization with this finite-amplitude conversion gives

\[
\lambda_{\rm res}=\lambda_{\rm radial}=1.071855028.
\]

The protected amplitude formulas then give

\[
Z_0=1.003199637,
\qquad
Z_F=0.991095202,
\]

\[
c_{d0}=3.541817665,
\qquad
c_{d1}=-0.874225883,
\]

\[
a_{d0}=0.261189843,
\qquad
a_{d1}=0.241781942.
\]

Together with the protected Pati–Salam relations

\[
h_Q=\frac{22}{21}h_{d0},\qquad
h_{u0}=\frac{23}{21}h_{d0},\qquad
h_{u1}=\frac1{21}h_{d1},
\]

this gives the following zero-continuous-fit flavor errors:

| Observable | Relative error |
|---|---:|
| \(m_c/m_t\) | \(-0.6111\%\) |
| \(m_u/m_t\) | \(-0.0305\%\) |
| \(m_s/m_b\) | \(+0.1720\%\) |
| \(m_d/m_b\) | \(-0.2491\%\) |
| \(|V_{us}|\) | \(+0.3545\%\) |
| \(|V_{cb}|\) | \(-0.1418\%\) |
| \(|V_{ub}|\) | \(+0.5341\%\) |

The maximum error is

\[
\boxed{0.6111\%}.
\]

The RMS error is

\[
0.3580\%.
\]

The independently fitted one-coefficient benchmark preferred \(\lambda_{\rm res}=1.07573\). The bounce-derived radial value differs by only

\[
-0.360\%.
\]

## 5. Comparison with other bounce-only normalizations

| Normalization | \(\lambda\) | Worst flavor error |
|---|---:|---:|
| unit energy operator | 1.000000 | 0.6321% |
| canonical radial anharmonic conversion | 1.071855 | **0.6111%** |
| squared energy conversion | 1.148873 | 0.6884% |
| true/center soft-mass ratio | 1.265178 | 0.8050% |
| square root of path tortuosity | 1.098744 | 0.6381% |
| path tortuosity | 1.207238 | 0.7469% |

The canonical radial conversion is the best of the predeclared direct bounce normalizations in RMS error and gives a sub-percent result without using flavor data.

## 6. Claim boundary

### Exact results

1. The incomplete-settling displacement is a pure bulk soft/radial mode to better than one part in \(10^{12}\) in norm fraction.
2. The center residual is the normalized bulk quartic energy,
   \[
   \epsilon_c=[1-(\Phi_c/\rho)^2]^2,
   \]
   across all 51 walls.
3. The finite-amplitude conversion
   \[
   \lambda_{\rm radial}=1/(1-u/2)
   \]
   follows exactly from the quartic radial potential.
4. Using this number in the previously frozen flavor theory gives a 0.6111% zero-continuous-fit realization.

### Remaining conditional statement

The seed action must identify its residual-response normalization with the canonical radial displacement normalization. The present symmetries still permit an additional order-one Wilson coefficient multiplying the normalized radial operator.

Therefore this checkpoint establishes a **bounce-derived canonical matching** but not an absolute uniqueness theorem. A different UV seed functional could multiply the same normalized radial operator by another order-one constant.

## 7. Scientific interpretation

The situation has improved from

\[
\text{one unidentified }O(1)\text{ number}
\]

to

\[
\text{one explicit bounce-only number generated by the radial anharmonicity,}
\]

with a clearly stated matching assumption.

The next decisive test is to derive the seed mass/source dependence on the canonical radial coordinate from an explicit local UV operator. If the seed threshold is linear in the canonically normalized radial response, \(\lambda_{\rm res}=\lambda_{\rm radial}\) is fixed. If the UV action allows an unrelated coefficient, that coefficient remains the final flavor-sector input.
