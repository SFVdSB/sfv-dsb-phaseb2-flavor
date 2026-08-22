# Phase B2 v0.3.0: Three Independent Operator/Flavor Extensions

## Scope and test discipline

This checkpoint tests the three proposed extensions **one at a time**.  No result from a later extension is inserted into an earlier test.

1. A genuine two-channel bulk/brane fermion mass matrix with the scalar-Hessian rotation connection.
2. One additional representation distinction in the generation/charge algebra.
3. A genuinely local wall-core mass invariant, independent of kinetic anisotropy.

The bounce, Higgs profile, seven MZ targets, and accepted overlap solver are unchanged.

## Extension 1: exact two-channel Hessian-rotation model

We take the local scalar-Hessian eigenbasis `U_H(x)` and define the fermion mass matrix

\[
\mathcal M_f(x)=U_H(x)
\begin{pmatrix}
 b(x)&0\\
 0&b(x)+\beta\Delta_H(x)
\end{pmatrix}
U_H^T(x),
\]

where

\[
\Delta_H(x)=\frac{\sqrt{\lambda_{\rm hard}(x)-\lambda_{\rm soft}(x)}}{\alpha},
\]

and `beta=1` is the natural benchmark.  In the rotating basis,

\[
\ell'=-b\ell+\theta_H' h,
\qquad
h'=-\theta_H'\ell-(b+\Delta_H)h.
\]

Writing `r=h/l`, the induced scalar correction in the light equation is

\[
C_{\rm geo}(x)=-\theta_H'(x)r(x).
\]

The common diagonal function `b(x)` cancels from the Riccati equation for `r`:

\[
r'=-\theta_H'-\Delta_H r-\theta_H' r^2.
\]

Therefore, in this minimal aligned construction the Hessian-rotation correction is universal: it cannot distinguish QL, uR, dR or the three generations.

### Baseline result

For `beta=1`:

- effective projection onto the canonical core mode: `k_geo = 0.199166`;
- correlation with `E(x)`: `0.812493`;
- heavy-channel norm fraction: `0.052155`;
- overlap between the exact total profile and its best scalar-core approximation: `0.998683`;
- leading local adiabatic approximation differs from the exact correction by about `30.6%` in shape norm.

The effect is real, modest, and well behaved.  It supplies a positive universal core envelope but not generation structure.

### 51-wall test

Across the 33-point local design and 18 corridor walls:

- mean `k_geo = 0.199189`;
- standard deviation `0.002112`;
- range `0.191629` to `0.207096`;
- mean heavy-channel fraction `5.215%`;
- mean exact/scalar profile overlap `0.998683`.

This is a robust geometric property of the wall family, not a benchmark-only accident.

### Flavor consequence

Refitting the original seven Route-I controls with the exact two-channel envelope still reproduces the targets, but the controls merely shift by a few percent.  Repeating the failed generation-plus-hypercharge two-spurion test does **not** improve it: the best bounded versions remain at roughly `46%` worst error.

**Extension-1 verdict:** physically successful coupled-channel effect; no solution to the flavor algebra.

## Extension 2: additional representation structure

The first attempt added only one extra down-sector or right-handed representation term to the hypercharge-based algebra, while retaining seven or fewer coefficients.

The best one-extra-spurion case was a down-sector slope term.  It reduced the worst error from the previous tens-to-hundreds of percent to:

\[
5.43\%\quad (|c_j|\le 5),
\]

but it hit coefficient bounds and remained outside the 1% target.

Thus one added representation distinction gets appreciably closer but does not close the model.

### Two independent sector projectors

A conventional representation-projector algebra using independent up- and down-singlet projectors was then tested:

\[
P_U=(0,1,0),\qquad P_D=(0,0,1)
\]

in the ordered representation space `(QL,uR,dR)`.  A seven-coefficient model can reproduce all seven observables to numerical precision with:

- all fitted coefficients below about `2.22` in magnitude;
- minimum kink strength `q_min = 0.871` on the selected well-localized branch;
- no need for near-delocalized or exponentially tiny charges.

This demonstrates that a standard local flavor algebra **can** encode the required sector pattern naturally.  It is not yet explanatory compression: seven coefficients still determine seven observables, and the projectors are assumed rather than derived.

**Extension-2 verdict:** one extra spurion is insufficient; two independent singlet-sector projectors are sufficient and physically allowed, but not yet predictive.

## Extension 3: local wall-core mass invariant

The large even-core term need not come from extreme kinetic anisotropy.  We tested local mass operators of the form

\[
\delta \mathcal L_M=-\bar\Psi\Psi\;c_A\,\mathcal I_{\rm core}(x),
\]

with several wall-derived invariants.

The preferred analytic candidate is the raw local gradient density

\[
\mathcal I_G(x)=
\frac{(\partial_y\Phi)^2+(\partial_y\phi)^2}
{\max[(\partial_y\Phi)^2+(\partial_y\phi)^2]}.
\]

Unlike the canonical reporting mode `E(x)`, this function is not parity-symmetrized.  It is the actual local wall invariant.

Refitting the unchanged seven-control sector structure with `I_G` gives:

- maximum observable error `5.2e-10%` (numerical zero);
- effective local core coefficients `kappa` between `-2.975` and `+3.727`;
- no active fit bounds;
- exact preservation of the original Yukawa condition numbers.

The partner-spectrum audit gives:

- `9/9` desired operators with exactly one near-zero chiral mode;
- `9/9` opposite operators with no near-zero mode;
- minimum unwanted eigenvalue `0.81128`.

Two other local candidates also work:

- normalized field-space gradient speed: exact fit with coefficients within about `1.71`, but it is a nonanalytic square-root EFT proxy;
- normalized portal cross-curvature `|4gPhi phi|`: exact fit with coefficients within about `4.26`.

Pure Hessian-mixing support alone misses at `1.53%`, so the wall-core operator must track gradient/portal strength more directly than just the mixing angle.

**Extension-3 verdict:** success.  A conventional local wall mass operator removes the need for enormous kinetic anisotropy while retaining exact flavor and chirality.

## Combined interpretation

The independent tests identify a much clearer architecture:

\[
\boxed{
\text{local gradient/portal mass core}
+\text{small universal Hessian-channel correction}
+\text{representation projectors}
}
\]

The wall supplies the shapes, scales, and a robust universal mixing correction.  The remaining generation information resides in the representation/flavor matrices.

What is now demonstrated:

1. The successful chiral equations have a local Hermitian operator origin.
2. The even core can be generated by an ordinary local wall invariant with order-few coefficients.
3. Hessian rotation produces a real 5% heavy-channel admixture and a stable universal core correction.
4. A conventional two-projector representation algebra can reproduce the required sector pattern with order-one coefficients.

What is not yet demonstrated:

1. Why nature chooses the specific up/down projector coefficients.
2. A reduction below seven continuously determined flavor coefficients.
3. A blind flavor prediction.

The next project should combine only the successful pieces, then ask whether the Phase-B wall formulas or a discrete/non-Abelian flavor symmetry can fix at least four of the seven projector coefficients before any new observable is examined.
