# Phase B2 v1.6.0 — Geometric Modulus and Brane-Breathing Audit

## Purpose

The v1.5.0 checkpoint isolated one remaining requirement in the frozen quark-flavor construction.  The residual radial driver

\[
\Sigma_{\rm radial}=\lambda_{\rm radial}\epsilon_c
=0.0670851132
\]

must dress the protected seed response and the fermion source as

\[
Z_0=e^{+\Sigma_{\rm radial}/21},
\qquad
Z_F=e^{-2\Sigma_{\rm radial}/15}.
\]

This checkpoint asks whether ordinary codimension-one brane embedding geometry—normal displacement, wall breathing, induced measure, and canonical field normalization—forces these weights without another symmetry.

It also tests a narrower alternative: an internal block-volume modulus of the already present Pati–Salam response space.  This is an internal field-space geometry, not a compact spacetime dimension.

## 1. The corrected bounce has a genuine breathing mode

For the coupled two-field background, the O(4) radial fluctuation operator is

\[
\mathcal L_\ell
=-\partial_r^2-\frac3r\partial_r
+\frac{\ell(\ell+2)}{r^2}
+H(r),
\]

where \(H(r)\) is the local two-field Hessian.  After the standard transformation \(q=r^{-3/2}u\), the numerical operator is

\[
-\partial_r^2+H(r)
+\frac{\ell(\ell+2)+3/4}{r^2}.
\]

The lowest eigenvalues are

\[
\lambda_{\ell=0,0}=-0.0978735,
\qquad
\lambda_{\ell=0,1}=+0.0214274,
\]

and

\[
\lambda_{\ell=1,0}=1.993\times10^{-4}.
\]

Thus the spectrum contains exactly the expected single negative O(4) breathing mode.  The nearly zero \(\ell=1\) state has correlation

\[
|\langle q_{\ell=1,0},\Phi'\rangle|=0.99652,
\]

confirming that it is the numerical translation zero mode.  The negative mode has correlation \(0.93990\) with the wall-translation tangent and \(0.67656\) with the simple dilation tangent \(r\Phi'\).

A family of reasonable interior-settling deformation windows has correlation \(0.638\)–\(0.844\) with the negative mode.  Therefore the incomplete-settling deformation has a substantial breathing component, but it is not identical to the geometric wall-breathing eigenmode.

## 2. Uniform brane breathing gives universal Weyl weights

A uniform normal displacement of a maximally symmetric wall acts locally as

\[
\gamma_{\mu\nu}\rightarrow e^{2\omega}\gamma_{\mu\nu}.
\]

For a canonical scalar seed \(X\), canonical normalization gives

\[
M_X^2\rightarrow e^{2\omega}M_X^2.
\]

If the protected mediator \(S\) is auxiliary, the source coefficient \(\mu SX^2\) receives the same factor.  Consequently,

\[
\boxed{
\frac{\mu}{M_X^2}\rightarrow\frac{\mu}{M_X^2}
}
\]

under ordinary uniform breathing.  This is the same protected ratio that survived the gauge beta-function audit.  Pure embedding geometry therefore predicts

\[
Z_0=1,
\]

not \(e^{\Sigma/21}\).

For the fermion source, canonical normalization produces ordinary engineering-dimension Weyl weights.  The complete possibilities tested were:

| worldvolume dimension | mediator treatment | \(d\ln Z_0/d\omega\) | \(d\ln Z_F/d\omega\) | relative mismatch with required pair |
|---:|---|---:|---:|---:|
| 3 | auxiliary | 0 | 2 | 33.6% |
| 3 | canonical scalar | −1/2 | 1 | 12.0% |
| 4 | auxiliary | 0 | 2 | 33.6% |
| 4 | canonical scalar | −1 | 0 | 94.2% |

No ordinary canonical-normalization case generates both required factors.

### Why normalized traces do not solve this

For

\[
\mathcal O_N=\frac1N\sum_{A=1}^N X_A^2,
\]

if every component transforms with the same geometric weight,

\[
X_A\rightarrow e^{-w\omega}X_A,
\]

then

\[
\mathcal O_N\rightarrow e^{-2w\omega}\mathcal O_N.
\]

The exponent is independent of \(N\).  The normalized trace removes the overall multiplicity from the operator normalization; it does not divide the geometric Weyl charge by 21 or 15.

This gives a clean no-go statement:

> Ordinary brane breathing supplies a genuine modulus, but it does not by itself generate the representation-normalized charges \(+1/21\) and \(-2/15\).

Extrinsic-curvature portals such as \(K\mathcal O_i\) can generate the desired response only after assigning separate curvature coefficients to the two sectors.  Those coefficients are the Wilson data the calculation was intended to eliminate.

## 3. A numerically interesting O(4)-shell coincidence

The O(4) solver has a three-dimensional Euclidean shell measure.  If one makes the additional, unproved identification

\[
\delta R=u,
\]

between the field-space fractional displacement

\[
u=0.134076019
\]

and a coordinate displacement, the shell-volume strain at the Hessian-mixing radius is

\[
\Sigma_{S^3}
=3\ln\left(1+\frac{u}{R_{\rm mix}}\right)
=0.0672580865.
\]

This differs from the required radial driver by only

\[
+0.2578\%.
\]

Inserted into the frozen flavor calculation it gives

\[
\max|\Delta O/O|=0.61385\%.
\]

This is striking numerically, but it is not a derivation.  The bulk field displacement and transverse coordinate displacement live in different spaces, and no existing term in the action identifies them.

The relation is also not a stable 51-wall law.  Using the same formula at \(R_{\rm mix}\):

- mean absolute driver discrepancy across all 51 walls: \(4.00\%\);
- maximum discrepancy: \(26.09\%\);
- mean absolute discrepancy in the local non-corridor 33-wall design: \(2.47\%\);
- local maximum: \(6.03\%\).

For the physical \((3+1)\)-dimensional worldvolume, the analogous four-dimensional volume strain is

\[
4\ln\left(1+\frac{u}{R_{\rm mix}}\right)
=0.0896774,
\]

which is \(33.68\%\) above the required driver.  It still happens to remain just below the broad 1% flavor threshold, but it does not reproduce the radial normalization.

Therefore the near equality in the O(4) shell is recorded as a numerical clue, not promoted to a law.

## 4. Minimal internal block-volume symmetry

Although external embedding geometry does not produce the charges, a very small internal geometric structure does.

Let the common radial driver \(\Sigma\) control two isotropic response blocks:

\[
C_{21}=e^{\Sigma/21}\mathbf 1_{21},
\]

for the complete 21-component Pati–Salam adjoint seed, and

\[
C_{15}=e^{-\Sigma/15}\mathbf 1_{15},
\]

for each fermion-source leg in the 15-dimensional \(SU(4)_C\) adjoint block.

Then

\[
\det C_{21}=e^{\Sigma},
\qquad
\det C_{15}=e^{-\Sigma},
\]

so

\[
\boxed{
\det C_{21}\det C_{15}=1.
}
\]

This is a unimodular block-volume constraint.  O(21) and O(15) isotropy distribute one total logarithmic volume deformation equally among the corresponding directions.

It gives immediately

\[
Z_0=e^{\Sigma/21},
\]

and because the fermion bilinear contains two source legs,

\[
Z_F=\left(e^{-\Sigma/15}\right)^2
=e^{-2\Sigma/15}.
\]

At the bounce value \(\Sigma=0.0670851132\), the exact outputs are

\[
Z_0=1.0031996371,
\qquad
Z_F=0.9910952029,
\]

and the determinant product equals one to numerical precision.

### Interpretation

This does not introduce 21, 15, or 36 compact spacetime dimensions.  It is a geometry of the internal kinetic/response metric of already present field channels.

The structure required is

\[
O(21)\times O(15)
\]

isotropy plus one shared unimodular volume modulus.  It is the most economical symmetry found so far that explains:

- why the seed response is divided by 21;
- why one fermion source leg is divided by 15;
- why two legs give \(-2/15\);
- why the signs are opposite.

## 5. What remains

The external brane calculation does not force the internal modulus to equal the existing bounce driver.  The final condition is

\[
\Sigma_{\rm internal}
=\Sigma_{\rm radial}
=\lambda_{\rm radial}\epsilon_c.
\]

This can be imposed without a continuous fitting coefficient if the internal response metric is *induced* by the canonical radial mode rather than represented by an independent scalar.  But that induced-metric identification is an additional geometric-compensator principle.

The present result is therefore:

\[
\boxed{
\begin{aligned}
&\text{ordinary brane embedding alone: insufficient,}\\
&\text{genuine breathing mode: verified,}\\
&\text{internal unimodular block-volume symmetry: exact charge derivation,}\\
&\text{radial-to-internal modulus identification: still to be established.}
\end{aligned}
}
\]

This is narrower than a complete parameter-free UV theorem, but it improves the target considerably.  The missing ingredient is no longer an arbitrary pair of portal coefficients.  It is one explicit structural question: whether the Pati–Salam seed/source response metric is induced by the same canonical radial mode that controls the bounce residual.

## Reproducibility

```bash
python src/geometric_modulus_embedding_audit.py
pytest -q tests/test_geometric_modulus_embedding.py
```
