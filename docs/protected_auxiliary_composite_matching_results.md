# Phase B2 v1.1.0 — protected auxiliary/composite matching

## Scope

The previous checkpoint established the one-loop gauge relation

\[
\frac{d}{d\ln Q}\ln\!\left(\frac{\mu_a}{M_a^2}\right)=0
\]

for each charged Pati–Salam seed block.  It did not yet prove that the physical mediator kernel actually uses the protected combination \(\mu/M^2\), rather than the excluded loop-bubble combination \(\mu^2/M^2\).

This checkpoint writes the zero-momentum auxiliary/Legendre matching explicitly, eliminates the collective channels, audits canonical normalization, and quantifies the additional condition required if the collective field is promoted to a propagating composite.

## 1. Gauge-invariant collective channels

The 22 seed components are represented by two gauge-invariant collective responses:

\[
Y_0\quad\leftrightarrow\quad X_0^2,
\]

and

\[
Y_{\rm A}\quad\leftrightarrow\quad
\frac1{21}\operatorname{Tr}_{\rm Adj}(X_{\cal A}X_{\cal A}).
\]

The second channel is one normalized adjoint **norm**, not 21 gauge-noninvariant scalar fields.  Its multiplicity is retained explicitly because the trace contains 21 equally normalized adjoint directions.

Define the source combinations

\[
u_0=S_1+S_2,
\qquad
u_{\rm A}=S_2.
\]

At zero momentum the local collective response is taken to satisfy

\[
M_0^2Y_0=\mu_0 u_0,
\qquad
M_{\rm A}^2Y_{\rm A}=\mu_{\rm A}u_{\rm A}.
\]

Thus

\[
Y_i=k_i u_i,
\qquad
k_i\equiv\frac{\mu_i}{M_i^2}.
\]

Equivalently, the zero-momentum Legendre action has inverse response coefficients

\[
\Gamma_{\rm coll}^{(0)}=
\frac12\frac{M_0^2}{\mu_0}Y_0^2-Y_0u_0
+21\left[
\frac12\frac{M_{\rm A}^2}{\mu_{\rm A}}Y_{\rm A}^2-Y_{\rm A}u_{\rm A}
\right].
\]

The overall Hubbard–Stratonovich sign or contour is conventional; the physical statement is the real zero-momentum response \(Y_i=k_i u_i\).

## 2. Eliminating the collective channels

The connected quadratic response in the mediator coordinates is

\[
K_S=
 k_0v_0v_0^T+21k_{\rm A}v_{\rm A}v_{\rm A}^T,
\]

where

\[
v_0=(1,1),\qquad v_{\rm A}=(0,1).
\]

Therefore

\[
K_S=
\begin{pmatrix}
k_0&k_0\\
k_0&k_0+21k_{\rm A}
\end{pmatrix}.
\]

At the protected equality point

\[
k_{\rm A}=k_0\equiv k,
\]

we obtain

\[
K_S=k
\boxed{
\begin{pmatrix}
1&1\\
1&22
\end{pmatrix}}
\]

and

\[
K_S^{-1}=\frac1{21k}
\begin{pmatrix}
22&-1\\
-1&1
\end{pmatrix}.
\]

The common response \(k\) changes the overall strength only.  It cancels from all representation ratios.

## 3. Flavor relations remain exact

Using the Pati–Salam source vector

\[
r_A=(1,-2T_{3R}^{(A)}),
\]

and the wall source \(e_1=(1,0)\), the cross contractions are

\[
e_1^TK_S^{-1}r_{d_R}=\frac1k,
\]

\[
e_1^TK_S^{-1}r_{Q_L}=\frac1k\frac{22}{21},
\]

\[
e_1^TK_S^{-1}r_{u_R}=\frac1k\frac{23}{21}.
\]

The family-odd up path gives

\[
e_2^TK_S^{-1}e_2=\frac1k\frac1{21}.
\]

Hence the frozen relations remain

\[
\boxed{
\frac{h_Q}{h_{d0}}=\frac{22}{21},\qquad
\frac{h_{u0}}{h_{d0}}=\frac{23}{21},\qquad
\frac{h_{u1}}{h_{d1}}=\frac1{21}.}
\]

The four-parameter flavor model remains at

\[
\max|\Delta O/O|=0.5957126\%.
\]

## 4. Canonical normalization does not spoil the result

For a general invertible change of mediator coordinates

\[
S=RS',
\]

the kernel and sources transform as

\[
K'=R^TKR,
\qquad
J'=R^TJ.
\]

The physical cross term is invariant:

\[
J_W^TK^{-1}J_F
=J_W'^TK'^{-1}J_F'.
\]

A 250-trial random numerical audit gave a maximum absolute discrepancy of

\[
9.06\times10^{-14}.
\]

Thus diagonalizing or canonically normalizing the two mediator coordinates cannot secretly turn the protected \(\mu/M^2\) response into \(\mu^2/M^2\).

## 5. One-loop protection survives the explicit matching

The inherited one-loop results are

\[
\frac{\mu_4/M_4^2}{\mu_0/M_0^2}=1,
\qquad
\frac{\mu_L/M_L^2}{\mu_0/M_0^2}=1,
\qquad
\frac{\mu_R/M_R^2}{\mu_0/M_0^2}=1.
\]

Consequently

\[
N_{\rm eff}
=15+3+3
=21
\]

throughout the audited running interval at one-loop gauge order.

This is the central positive result:

> The coefficient used by the explicit zero-momentum auxiliary matching is exactly the one-loop protected coefficient.

## 6. Why this is not the excluded derivative bubble

Two physically different calculations must not be conflated.

### Auxiliary/Legendre response

\[
k_i^{(0)}=\frac{\mu_i}{M_i^2}.
\]

This is the static source-to-response coefficient.  It enters once and is protected.

### Propagating scalar-loop derivative susceptibility

\[
k_i^{(2)}\propto\frac{\mu_i^2}{M_i^2}.
\]

This contains two trilinear vertices and was shown in v1.0.0 to enhance the adjoint contribution to \(N_{\rm eff}\simeq23.73\), producing a 6.5% flavor error.

The auxiliary construction is not a rewriting of that loop.  It is a distinct zero-momentum EFT matching.

## 7. Propagating-composite caveat

If the collective coordinate becomes dynamical, write its response as

\[
k_i(p)=\frac{\mu_i}{M_i^2+Z_ip^2}.
\]

The exact 21 remains valid at all momenta provided

\[
\frac{Z_i}{M_i^2}
\]

is common across the singlet and adjoint channels.

If it is not common, finite transverse momentum changes \(N_{\rm eff}\).  With a common relative mismatch

\[
\rho=\frac{(Z/M^2)_{\rm adj}}{(Z/M^2)_0}=1+\delta,
\]

we obtain

\[
N_{\rm eff}(q)=21\frac{1+q^2}{1+(1+\delta)q^2},
\qquad q=p/M_{\rm comp}.
\]

The inherited sub-1% interval gives the following allowed mismatch ranges:

| \(q\) | Allowed \(\delta\) |
|---:|---:|
| 0.5 | \(-3.78\%,+12.45\%\) |
| 1.0 | \(-1.51\%,+4.98\%\) |
| 2.0 | \(-0.95\%,+3.11\%\) |
| 3.0 | \(-0.84\%,+2.77\%\) |
| 3.6 | \(-0.81\%,+2.68\%\) |

Using the Phase A inverse wall widths and identifying the composite mass with the inherited matching scale would give approximately

\[
q_{\rm grad}\simeq2.41,
\qquad
q_{\rm action}\simeq3.59.
\]

Therefore a propagating composite at that relatively low mass requires percent-level universality of \(Z/M^2\).  There are two clean ways around this:

1. retain the strict nondynamical auxiliary/Legendre realization, for which \(Z=0\);
2. place the composite mass well above the inverse wall width, or derive an O(22)-universal derivative term.

## 8. What has and has not been proved

### Demonstrated

- A gauge-invariant zero-momentum collective matching exists.
- Its physical coefficient is \(\mu/M^2\), not \(\mu^2/M^2\).
- The one-loop gauge beta functions protect that coefficient.
- The exact Gram matrix, determinant 21, and rational flavor ratios survive.
- Canonical normalization and mediator-basis changes do not alter the result.
- The four-parameter flavor result remains at 0.596%.

### Still conditional

- The auxiliary/Legendre field has not yet been derived as a specific propagating non-supersymmetric bound state.
- A propagating composite needs controlled derivative normalization.
- Non-gauge scalar interactions that break O(22) could split \(\mu/M^2\) and require a separate beta-function audit.
- The four absolute flavor amplitudes remain to be derived.

## Verdict

\[
\boxed{
\text{The protected coefficient survives explicit auxiliary matching.}}
\]

The model does not fail at this checkpoint.  The strict auxiliary realization is internally complete at the EFT level.  A fully dynamical composite completion remains possible but is more constrained: its derivative term must either respect the same O(22) structure or lie above the wall-resolution scale.

The next physics task is to derive the four remaining absolute amplitudes from the common response scale, wall normalization, mediator masses, and flavor-breaking vacuum expectation values.
