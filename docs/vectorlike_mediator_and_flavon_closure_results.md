# Phase B2 v0.6.0 — Vectorlike mediator and flavon closure

## Scope and claim boundary

The wall-dependent fermion operator was frozen before this test:

\[
\mathcal M_{Ai}(y)=q_{Ai}O(y)+q_{Ai}h_{Ai}\,\mathcal I_G(y)+C_{\rm geo}(y),
\]

where \(O\) is the kink-like wall mode, \(\mathcal I_G=G/G_{\max}\) with
\(G=(\partial_y\Phi)^2+(\partial_y\phi)^2\), and \(C_{\rm geo}\) is the independently calculated universal Hessian-channel envelope.

This checkpoint asks whether the remaining representation/generation coefficients can arise by integrating out a small, ordinary heavy-fermion sector.

The tree-level matching and matrix-rank statements below are conventional EFT results.  The small rational values

\[
\frac{22}{21},\qquad \frac{23}{21},\qquad \frac{1}{21},\qquad \frac{15}{14}
\]

were discovered after examining the successful coefficients.  They are now frozen *targets for a UV representation or mediator mass matrix*, not claimed first-principles group factors.

## 1. Explicit heavy-vectorlike matching

For each Standard Model representation \(A\in\{Q,u,d\}\), introduce a heavy vectorlike partner \(F_A\):

\[
\mathcal L_{\rm UV}\supset
\bar F_A(i\!\not\!D-M_A)F_A
-\lambda_{L,A}\sqrt{G(y)}\,\bar\psi_{Ai,L}F_{A,R}
-\sqrt{G(y)}\,\bar F_{A,L}
\left(\mu_{A0}+n_i\mu_{A1}\right)\psi_{Ai,R}
+\text{h.c.},
\]

with

\[
n_i=(-1,0,+1).
\]

At energies below \(M_A\), the leading equation of motion is

\[
F_A\simeq M_A^{-1}\,\sqrt{G}\,(\cdots)\psi_A.
\]

Substitution gives the local wall operator

\[
\delta\mathcal L_{\rm eff}
\supset
-c_{Ai}G(y)\bar\psi_{Ai,L}\psi_{Ai,R}+\text{h.c.},
\]

with

\[
\boxed{
 c_{Ai}=b_A+s_A n_i,
 \qquad
 b_A=\frac{\lambda_{L,A}\mu_{A0}}{M_A},
 \qquad
 s_A=\frac{\lambda_{L,A}\mu_{A1}}{M_A}.
}
\]

A symmetry can set \(b_Q=0\), producing the pure generation-odd \(Q_L\) pattern found numerically.  The down-sector kink factor

\[
q_{d,i}=\exp(a_{d0}+a_{d1}n_i)
\]

can arise from repeated heavy-mediator/flavon insertions or exponential canonical normalization.  This establishes a standard local UV route to the complete *form* of the effective operator.

The construction is four-dimensionally Lorentz invariant at tree level.  Its variation is transverse to the brane, so it remains compatible with the earlier charged-sector Lorentz-invariance assumption.

## 2. How many flavor directions are actually required?

The fitted core-coupling matrix before the down-sector \(q\) prefactor is

\[
H_{Ai}=
\begin{pmatrix}
-h_Q&0&h_Q\\
h_{u0}-h_{u1}&h_{u0}&h_{u0}+h_{u1}\\
h_{d0}-h_{d1}&h_{d0}&h_{d0}+h_{d1}
\end{pmatrix}.
\]

Its singular values are

\[
6.88422,\qquad 3.97553,\qquad 7.25\times10^{-17}.
\]

It therefore has **exact rank two**.  The two required generation directions are simply

\[
\mathbf 1=(1,1,1),
\qquad
\mathbf n=(-1,0,+1).
\]

This is a major simplification: no third independent generation function is needed for the real flavor hierarchy.  The difficult information is in the representation-dependent matching coefficients \(b_A,s_A\), not in a complicated generation-space function.

The down-kink matrix also has rank two because the \(Q_L\) and \(u_R\) rows are universal and only the \(d_R\) row carries the exponential flavon weight.

## 3. Strict one-mediator test

A genuinely rank-one hypothesis forces every representation to use the same generation vector:

\[
h_{Ai}=s_A(1+r n_i).
\]

It has six continuous parameters after adding the two down-kink parameters.  Its best benchmark result is

\[
\boxed{25.55\%\text{ maximum observable error}.}
\]

The optimization also drives two sector amplitudes to the imposed limit \(8\).  One rank-one mediator cannot simultaneously generate:

- a purely odd \(Q_L\) core;
- an almost generation-independent \(u_R\) core;
- a more strongly sloped \(d_R\) core.

Thus **one shared generation channel is ruled out**.  A flavor-singlet channel plus one traceless flavon channel is the minimum.

## 4. One frozen mediator relation

The exact benchmark gives

\[
\frac{h_{u1}}{h_{d1}}=0.0479527.
\]

The simple frozen ratio

\[
\boxed{h_{u1}=\frac{1}{21}h_{d1}}
\]

differs by only \(0.696\%\) at the coefficient level.

With that one UV relation fixed and the other six quantities fitted, the seven observables are reproduced with

\[
\boxed{0.0926\%\text{ maximum error}.}
\]

More importantly, the same six-parameter relation was refitted independently on all 51 available wall solutions:

- 51/51 points remain below 1%;
- mean worst error: \(0.160\%\);
- median worst error: \(0.0959\%\);
- worst case: \(0.918\%\), at the extreme \(-8\%\) stiffness corridor point.

This is the strongest new result.  The ratio \(1/21\) was found post hoc, but its successful 51-wall survival makes it a serious target for an explicit UV Clebsch factor or mediator-mass relation.

A purely geometric propagator estimate,

\[
\frac{h_{u1}}{h_{d1}}
\simeq
\frac{\xi_{\Phi,F}}{R_{\rm peak}}
=0.0467514,
\]

also works well:

- benchmark maximum error: \(0.355\%\);
- 50/51 walls below 1%;
- worst case: \(1.15\%\).

This is suggestive of a finite-size/heavy-propagator suppression, but its cross-wall correlation is weaker than the fixed \(1/21\) relation.

## 5. Five- and four-parameter mediator models

Restoring the raw gradient invariant revealed additional near-rational coefficient ratios:

\[
\frac{h_Q}{h_{d0}}\simeq\frac{22}{21},
\qquad
\frac{h_{u0}}{h_{d0}}\simeq\frac{23}{21},
\qquad
\frac{h_{u1}}{h_{d1}}\simeq\frac{1}{21}.
\]

At the benchmark their relative mismatches are only

\[
-0.060\%,\qquad +0.196\%,\qquad -0.696\%.
\]

### Five parameters

Impose

\[
h_{u0}=\frac{23}{21}h_{d0},
\qquad
h_{u1}=\frac{1}{21}h_{d1},
\]

while retaining independent \(h_Q,h_{d0},h_{d1},a_{d0},a_{d1}\).  The result is

\[
\boxed{0.1068\%\text{ maximum error}.}
\]

Across the 51 walls, 47 remain below 1%.  The failures occur only in the far \(\pm5\%\) and \(\pm8\%\) stiffness corridor, where the entire wall is physically different.

### Four parameters

Impose the additional core relation

\[
h_Q=\frac{22}{21}h_{d0}.
\]

Only

\[
h_{d0},\quad h_{d1},\quad a_{d0},\quad a_{d1}
\]

remain adjustable.  The benchmark result is

\[
\boxed{0.5957\%\text{ maximum error}.}
\]

This reaches the preregistered 1% strong-realization threshold with four parameters for seven observables.  Across all 51 altered walls, 33 remain below 1%; the relation is useful locally but not universal over the full \(\pm8\%\) stiffness corridor.

### Three parameters

Finally imposing

\[
\frac{a_{d0}}{a_{d1}}=\frac{15}{14}
\]

gives

\[
1.277\%\text{ maximum error}.
\]

This is a near miss.  The present compression frontier is therefore **four continuous parameters**, not three.

## 6. Chiral spectrum

The six-, five-, and four-parameter mediator models all retain:

- 9/9 desired operators with exactly one near-zero chiral state;
- 9/9 opposite-chirality operators with no near-zero state;
- minimum opposite-sector eigenvalue \(0.87688\).

The compression does not work by weakening chirality or introducing extra light partners.

## 7. What has and has not been derived

### Established

1. The frozen wall operator can be generated by a local Hermitian vectorlike-fermion action.
2. Exactly two generation directions, \(\mathbf1\) and \(\mathbf n\), span the required real flavor structure.
3. One rank-one mediator is insufficient.
4. A singlet plus one traceless flavon direction is sufficient.
5. One frozen relation, \(h_{u1}/h_{d1}=1/21\), reduces the model to six parameters and survives all 51 wall tests below 1%.
6. A four-parameter mediator/Clebsch model reproduces the benchmark within 0.60% while retaining the exact chiral spectrum.

### Not established

1. No known symmetry or group representation has yet derived the numbers \(22/21,23/21,1/21\).
2. The rational relations were discovered post hoc.
3. The four-parameter result is not yet a blind prediction.
4. The heavy mediator masses, gauge representations, and flavon potential have not yet been uniquely fixed.

## 8. Interpretation and next target

The outcome favors a conventional explanation rather than a missing bounce equation:

\[
\boxed{
\text{wall geometry}
+\text{two heavy mediator channels}
+\text{one traceless flavon direction}
\longrightarrow
\text{observed real flavor pattern}.
}
\]

The next sharp task is no longer to search arbitrary flavor formulas.  It is to find the smallest anomaly-safe gauge/flavor representation or heavy mass matrix that produces the frozen Clebsch targets

\[
\boxed{
C_Q:C_u:C_d=22:23:21,
\qquad
C_{u,\rm odd}:C_{d,\rm odd}=1:21.
}
\]

Those ratios provide a concrete algebraic target for a discrete or unified mediator sector.  They must be derived without changing them after examining a new observable.
