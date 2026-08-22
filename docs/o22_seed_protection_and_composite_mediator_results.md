# Phase B2 v0.9.0 — O(22) seed protection, gauge splitting, and composite mediator audit

## 1. What the numbers 21 and 22 physically mean

The Pati–Salam gauge group is a direct product with gauge-algebra dimensions

\[
\dim SU(4)_C+\dim SU(2)_L+\dim SU(2)_R=15+3+3=21.
\]

A field in the complete adjoint therefore contains exactly 21 real internal components, one for every independent generator. This is the physical reason for **21 once a complete adjoint seed is chosen**. Pati–Salam does not force the UV theory to contain that seed; its presence remains a model-building hypothesis.

Adding the singlet seed gives

\[
1\oplus\operatorname{Adj}(G_{PS}),\qquad 1+21=22.
\]

This is a 22-dimensional **internal field space**, not 22 spacetime dimensions. It has no direct relation to the 10- or 11-dimensional spacetime used in string/M-theory constructions.

## 2. O(22) boundary condition

Collect the seed fields into

\[
X^I=(X_0,X_4^a,X_L^i,X_R^j),\qquad I=1,\ldots,22.
\]

At a seed-unification scale \(\Lambda_{seed}\), impose

\[
\mathcal L_{O(22)}\supset
\frac12(\partial X^I)^2
-\frac12 M_X^2 X^I X^I
-\mu_2 S_2 X^I X^I.
\]

These terms treat all 22 components identically. The controlled rank-one spurion

\[
-\mu_1 S_1 X_0^2
\]

selects the singlet and creates the incidence vector \((1,1)\). Pati–Salam gauges the remaining 21 directions and explicitly breaks the global seed-space O(22).

## 3. Gauge-matching inputs

Using the inherited SFV matching scale and high-scale gauge values,

- \(M_{match}=1.6280e+13\,\mathrm{GeV}\),
- \(\Lambda_{seed}=2.4100e+14\,\mathrm{GeV}\),
- \(\ln(\Lambda_{seed}/M_{match})=2.694860\),
- \(g_4=0.580379\), \(g_L=0.541853\), \(g_R=0.523058\).

The Pati–Salam matching used

\[
g_Y=\sqrt{3/5}g_1,\qquad
\frac1{g_Y^2}=\frac1{g_R^2}+\frac2{3g_4^2}.
\]

## 4. Leading-log viability audit

A full seed-sector RGE requires all scalar quartics, trilinears, and masses. The present checkpoint therefore freezes a transparent leading-log marker,

\[
\frac{k_a}{k_0}=
\exp\left[
 s c_g\frac{C_2(\operatorname{Adj}_a)g_a^2}{16\pi^2}
 \ln\frac{\Lambda_{seed}}{M_{match}}
\right],
\]

where \(s=\pm1\) records the unresolved sign and \(c_g\) is the model-specific coefficient.

For \(c_g=1\):

- suppressive sign: \(N_{eff}=20.6012\), mismatch \(-1.899\%\), interpolated worst flavor error \(0.893\%\);
- enhancing sign: \(N_{eff}=21.4073\), mismatch \(+1.939\%\), outside the sub-1% flavor interval.

The allowed coefficient ranges at the existing scale interval are

\[
c_g\lesssim 1.283\quad(s=-1),\qquad
c_g\lesssim 0.395\quad(s=+1)
\]

for the 1% flavor criterion.

For \(c_g=1\), the maximum allowed scale ratios are approximately

\[
\left(\frac{\Lambda}{M}\right)_{max}
=31.72\quad(s=-1),\qquad
2.90\quad(s=+1).
\]

The current ratio is \(14.80\). Thus an O(1) suppressive correction is compatible with the 1% tolerance, while an enhancing correction would require a much shorter running interval or additional protection.

This is a viability result, not a completed RGE derivation. General gauge-theory RGEs for dimensionful scalar parameters are model dependent and require the complete interaction content.

## 5. Composite/induced mediator alternative

Define two gauge-invariant seed bilinears

\[
\mathcal O_1=X_0^2,\qquad
\mathcal O_2=X_0^2+\sum_{\mathcal A=1}^{21}X_\mathcal A^2.
\]

Their leading two-point susceptibility is proportional to

\[
\Pi_{ab}\propto
\begin{pmatrix}1&1\\1&22\end{pmatrix}
\]

when all seed propagators share one normalization. Auxiliary/composite mediator fields \(S_a\) introduced for these bilinears therefore inherit the required matrix directly.

The advantage is important: a compositeness boundary condition can set the unrelated bare \(S\)-mass matrix to zero, so the mediator matrix is induced by the seed correlator rather than freely chosen. Gauge running still perturbs the equality of seed kernels, but it cannot generate an arbitrary new flavor matrix without additional spurions.

## 6. Verdict

1. **There is a physical reason for 21:** it is the exact number of internal directions in one complete Pati–Salam adjoint.
2. **There is no extra-spacetime claim:** 22 is the dimension of an internal seed multiplet, not the dimension of the universe.
3. **O(22) can protect the common boundary condition:** common mass and \(S_2\) coupling are symmetry enforced at \(\Lambda_{seed}\).
4. **Pati–Salam gauging breaks it:** the size is loop suppressed and, under an O(1) suppressive marker over the existing SFV interval, remains inside the 1% flavor tolerance.
5. **The induced/composite route is preferred:** it removes the arbitrary bare mediator matrix and makes the \(21\) Gram structure the seed susceptibility.
6. **Remaining work:** derive the exact sign and coefficient from the full scalar-sector RGEs, then derive the four absolute amplitudes.
