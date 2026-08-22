# Exact Fermion-Wall Operator Reduction and Minimal Two-Spurion Test

## Scope

This checkpoint follows the initial Phase-B2 operator audit. It derives the chiral zero-mode equations from a Hermitian local fermion action before attempting to identify the effective Route-I controls. It then tests the predeclared smallest flavor algebra consisting of one generation spurion and one Standard-Model representation spurion.

The calculation does **not** modify the bounce, Higgs interface, target observables, or the accepted chiral-overlap solver.

## 1. Hermitian local action

For commuting positive kinetic weights, take

\[
S_\Psi=\int d^4x\,dy\left[
\frac{i}{2}Z_\parallel(y)\,\bar\Psi\gamma^\mu\!\stackrel{\leftrightarrow}{D}_\mu\Psi
+\frac{i}{2}Z_\perp(y)\,\bar\Psi\Gamma^y\!\stackrel{\leftrightarrow}{D}_y\Psi
-\bar\Psi W(y)\Psi
\right].
\]

Canonical normalization of the four-dimensional kinetic term uses

\[
\Psi=Z_\parallel^{-1/2}\chi,
\qquad
A(y)=\frac{Z_\perp(y)}{Z_\parallel(y)},
\qquad
M(y)=Z_\parallel^{-1/2}WZ_\parallel^{-1/2}.
\]

With the convention \(i\Gamma^y=\gamma^5\), the canonically normalized equation is

\[
\left[
i\gamma^\mu D_\mu
+\gamma^5\left(A D_y+\frac12A'\right)
-M
\right]\chi=0.
\]

If a local internal basis \(U(y)\) is used, then

\[
D_y=\partial_y+\mathcal A_y,
\qquad
\mathcal A_y=U^\dagger U'.
\]

For a massless four-dimensional mode,

\[
f_L'=-\mathcal A_y f_L-\frac12(\ln A)'f_L-A^{-1}Mf_L,
\]

\[
f_R'=-\mathcal A_y f_R-\frac12(\ln A)'f_R+A^{-1}Mf_R.
\]

The sign of the mass term reverses between the two chiralities. Reversing the transverse coordinate or the sign convention for the right-handed wall mass changes the displayed signs but not this invariant statement.

## 2. Correction to the initial kinetic-kernel interpretation

A common position-dependent factor

\[
Z_\parallel=Z_\perp=K(y)
\]

gives \(A=1\). Its derivative cancels under proper Hermitian canonical normalization. Therefore, an isotropic scalar kinetic weight cannot by itself generate the canonical even wall-core term.

There are three conventional ways to obtain such a term:

1. transverse/brane kinetic anisotropy, \(Z_\perp/Z_\parallel\ne1\);
2. an explicit local wall-core contribution in the mass kernel \(W\);
3. coupled-channel mixing followed by controlled elimination of a heavier channel.

The strong correlation found in v0.1.0 remains a useful shape observation, but it cannot be interpreted as arising from a single isotropic kinetic factor.

## 3. Exact embedding of the Route-I equation

Let

\[
F_E'(x)=E(x),
\]

and for each sector and generation define

\[
\kappa_{Ai}=q_{Ai}h_{Ai},
\qquad
A_{Ai}(x)=\exp\!\left[2\kappa_{Ai}F_E(x)\right],
\]

\[
W_{Ai}(x)=A_{Ai}(x)q_{Ai}O(x).
\]

Then

\[
A_{Ai}^{-1}W_{Ai}+\frac12(\ln A_{Ai})'
=q_{Ai}O+\kappa_{Ai}E
=q_{Ai}(O+h_{Ai}E).
\]

Thus the accepted scalar Route-I equation is exactly embeddable in a local anisotropic action. The numerical reconstruction error is

\[
8.88\times10^{-16}.
\]

This establishes **operator existence**, not naturalness. If the entire even term is assigned to kinetic anisotropy, the outer \(Q_L\) generations require

\[
\frac{A_{\max}}{A_{\min}}\simeq3.67\times10^4.
\]

That is too large to accept without a microscopic mechanism. A more plausible final operator may split the even-core contribution between an explicit localized mass invariant and a smaller kinetic anisotropy.

## 4. Hessian/Berry connection

For the real two-field Hessian rotation

\[
U(\theta)=
\begin{pmatrix}
\cos\theta&-\sin\theta\\
\sin\theta&\cos\theta
\end{pmatrix},
\]

one obtains

\[
\mathcal A_y=U^TU'=
\begin{pmatrix}
0&-\theta'\\
\theta'&0
\end{pmatrix}.
\]

The connection is off diagonal and antisymmetric. It has no diagonal scalar component in this real basis.

A numerical basis-covariance test gives:

- including \(\mathcal A_y\): relative maximum discrepancy \(7.10\times10^{-9}\);
- omitting \(\mathcal A_y\): relative maximum discrepancy \(0.387\).

Therefore, the connection is real and necessary, but it is not an independent scalar potential that may simply be replaced by \(|\theta'|\).

The raw shape \(|\theta'|\) correlates with the canonical core profile at \(0.92484\). A leading adiabatic-elimination candidate,

\[
\frac{(\theta')^2}{\Delta_H},
\]

where \(\Delta_H\) is the Hessian eigenvalue gap, correlates at \(0.77191\). This weaker but still nontrivial relation is the physically safer quantity because a real off-diagonal connection affects a retained channel only after coupling to and eliminating another channel.

## 5. Minimal two-spurion model

The predeclared flavor algebra uses

\[
F=\operatorname{diag}(-1,0,+1)
\]

for generation and Standard-Model hypercharge

\[
Y_Q=\frac16,\qquad Y_u=\frac23,\qquad Y_d=-\frac13
\]

for representation dependence. The six-coefficient bilinear is

\[
q_{Ai}=\exp\left(a_SY_A+a_{FS}Y_An_i\right),
\]

\[
\kappa_{Ai}=b_0+b_Fn_i+b_SY_A+b_{FS}Y_An_i,
\]

with

\[
B_{Ai}=q_{Ai}O+\kappa_{Ai}E.
\]

This is one parameter fewer than the seven-control Route-I realization and has a clear operator/spurion interpretation.

### Results

| Test | Maximum observable error | Interpretation |
|---|---:|---|
| \(|c_j|\le3\) | 156.55% | fails decisively with order-one coefficients |
| \(|c_j|\le5\) | 18.47% | still fails; coefficients hit bounds and profile ratios become extreme |
| seeded unbounded diagnostic | 5.54% | not acceptable; \(h_{\max}\simeq5980\) and \(q\) spans \(0.0016\) to \(52\) |

The unbounded result is not a controlled global optimum; its purpose is to show the direction in which the model tries to improve. It does so by generating pathological effective charges rather than a natural flavor explanation.

## 6. Interpretation

The checkpoint separates two questions that had previously been mixed together:

1. **Can a conventional local action produce the successful chiral equations?**  
   Yes. An exact local operator reduction exists.

2. **Does the smallest conventional flavor algebra determine the needed generation pattern?**  
   No. One generation generator plus hypercharge is insufficient.

The failure is therefore in the compact flavor-charge assignment, not in chiral localization or in the wall background.

## 7. Next controlled extension

The next model should add only one new structural ingredient at a time:

1. an actual two-channel bulk/brane fermion mass matrix, allowing the Hessian connection to contribute through adiabatic elimination;
2. one additional independent representation spurion, or a small non-Abelian/discrete flavor generator, to distinguish the \(Q_L\), \(u_R\), and \(d_R\) slope/intercept patterns;
3. a local mass-core invariant so that the large \(Q_L\) even term need not be generated entirely by kinetic anisotropy.

No CP phases should be introduced until this real operator structure is frozen.
