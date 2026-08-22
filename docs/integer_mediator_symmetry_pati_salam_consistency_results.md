# Phase B2 v0.7.0 — Integer mediator, flavor symmetry, Pati–Salam, and consistency audit

## Scope and claim boundary

The wall operator and the three rational targets were frozen before this checkpoint:

\[
\frac{h_Q}{h_{d0}}=\frac{22}{21},\qquad
\frac{h_{u0}}{h_{d0}}=\frac{23}{21},\qquad
\frac{h_{u1}}{h_{d1}}=\frac1{21}.
\]

This checkpoint tests five questions, in the declared order:

1. Can an ordinary heavy-mediator mass matrix generate the ratios?
2. Can a flavor symmetry enforce the required generation direction?
3. Can a known gauge-representation contraction explain the sector pattern?
4. Is the candidate gauge invariant and anomaly safe?
5. Does it preserve the existing Lorentz-violation and strong-CP architecture?

The rational targets were discovered post hoc in v0.6.0. Nothing below converts that historical fact into a blind prediction. The new result is that they now admit one compact, explicit and testable UV candidate.

## 1. Integer two-channel mediator mass matrix

Consider two heavy vectorlike channels with dimensionless mass matrix

\[
\widehat M_N=
\begin{pmatrix}
1&1\\
1&N+1
\end{pmatrix},
\qquad \det\widehat M_N=N.
\]

It is positive definite for every positive integer \(N\), with inverse

\[
\widehat M_N^{-1}=\frac1N
\begin{pmatrix}
N+1&-1\\
-1&1
\end{pmatrix}.
\]

All light-heavy coupling vectors below contain only \(0,+1,-1\). Let \(e_1=(1,0)\), \(e_2=(0,1)\), and choose

\[
r_d=(1,1),\qquad r_Q=(1,0),\qquad r_u=(1,-1).
\]

Then

\[
e_1^T\widehat M_N^{-1}r_d=1,
\]

\[
e_1^T\widehat M_N^{-1}r_Q=\frac{N+1}{N},
\]

\[
e_1^T\widehat M_N^{-1}r_u=\frac{N+2}{N}.
\]

A second contraction gives

\[
e_2^T\widehat M_N^{-1}e_2=\frac1N.
\]

Therefore one matrix produces

\[
\boxed{
\frac{h_Q}{h_{d0}}=\frac{N+1}{N},\qquad
\frac{h_{u0}}{h_{d0}}=\frac{N+2}{N},\qquad
\frac{h_{u1}}{h_{d1}}=\frac1N.
}
\]

For \(N=21\), these are exactly the frozen \(22/21\), \(23/21\), and \(1/21\) relations.

The numerical matrix is

\[
\widehat M_{21}=
\begin{pmatrix}1&1\\1&22\end{pmatrix},
\]

with eigenvalues

\[
0.9524884,\qquad22.0475116,
\]

and condition number \(23.15\). It is hierarchical but neither singular nor pathologically tuned.

### Finite minimality search

An exact finite search was performed over symmetric positive-definite integer \(2\times2\) matrices with nonnegative entries and light-heavy coupling vectors in \(\{-1,0,1\}^2\). The three \(21{:}22{:}23\) channels were required to share one left coupling vector, and a separate \(1{:}21\) contraction was required.

No matrix with maximum entry \(\leq21\) succeeds in this declared class. The first solutions appear at maximum entry \(22\). This is a computational minimality statement for the stated search class, not a theorem covering arbitrary UV theories.

## 2. Does the flavor symmetry work?

The required generation direction remains

\[
F=\operatorname{diag}(-1,0,+1).
\]

A gauged family symmetry \(U(1)_F\) can enforce it by assigning the three complete fermion families charges

\[
q_F=(-1,0,+1).
\]

Using one complete Standard Model family plus a right-handed neutrino in left-handed Weyl notation, the full three-family anomaly ledger gives

\[
[SU(3)_c]^2U(1)_F=0,
\]

\[
[SU(2)_L]^2U(1)_F=0,
\]

\[
[U(1)_Y]^2U(1)_F=0,
\]

\[
U(1)_Y[U(1)_F]^2=0,
\]

\[
[U(1)_F]^3=0,
\qquad
[\mathrm{gravity}]^2U(1)_F=0.
\]

The cancellation occurs because the coefficients linear in family charge are proportional to

\[
-1+0+1=0,
\]

while the cubic anomaly is proportional to

\[
(-1)^3+0^3+(+1)^3=0.
\]

There are 12 \(SU(2)_L\) doublets and, in a Pati–Salam completion, 12 \(SU(2)_R\) doublets, so the global Witten anomalies are also absent.

A discrete \(Z_k\) remnant obtained by Higgsing this anomaly-free gauged \(U(1)_F\) has a consistent UV parent. The symmetry explains the traceless family direction, but **does not by itself explain why \(N=21\)**.

## 3. Pati–Salam representation contraction

The Pati–Salam group is

\[
G_{PS}=SU(4)_C\times SU(2)_L\times SU(2)_R,
\]

with fermion multiplets

\[
(4,2,1)\oplus(\bar4,1,2).
\]

The right-handed up and down channels form an \(SU(2)_R\) doublet:

\[
T_{3R}(u_R)=+\frac12,\qquad
T_{3R}(d_R)=-\frac12,
\]

while \(Q_L\) is an \(SU(2)_R\) singlet with \(T_{3R}=0\).

The representation factor

\[
1+2T_{3R}
\]

therefore gives

\[
Q_L:1,\qquad u_R:2,\qquad d_R:0.
\]

Adding a universal integer \(N\) produces

\[
Q_L:N+1,\qquad u_R:N+2,\qquad d_R:N.
\]

For \(N=21\), this is exactly

\[
\boxed{Q_L:u_R:d_R=22:23:21.}
\]

### Why 21 is especially notable here

Pati–Salam contains two independent group-theoretic appearances of 21:

1. Its total gauge-algebra dimension is
   \[
   15+3+3=21.
   \]

2. The sum of quadratic Casimirs for either fermion multiplet is
   \[
   C_2(4)+C_2(2)=\frac{15}{8}+\frac34=\frac{21}{8}.
   \]
   Thus
   \[
   8[C_2(4)+C_2(2)]=21.
   \]

This makes Pati–Salam a far more relevant candidate than an arbitrary group that merely has a 21-dimensional representation. It simultaneously supplies the integer and the exact \(T_{3R}\) sector offsets.

It is still a **conditional derivation**: a concrete mediator calculation must show that the \(N\) entry in \(\widehat M_N\) is proportional to the Pati–Salam invariant or to a full adjoint/multiplicity contribution. Group-number coincidence alone is not enough.

## 4. Gauge invariance and anomaly safety

The candidate passes the structural audit under ordinary assignments:

- The wall gradient invariant
  \[
  G(y)=(\partial_y\Phi)^2+(\partial_y\phi)^2
  \]
  is a Standard Model and Pati–Salam singlet.
- Each heavy mediator can be introduced as a vectorlike pair in the same gauge representation as the light channel.
- Vectorlike pairs make no net contribution to chiral gauge anomalies.
- The singlet and \(T_{3R}\)-adjoint contractions are gauge invariant when the breaking spurion is assigned the corresponding Pati–Salam representation.
- The family \(U(1)_F\) anomaly ledger vanishes exactly.
- The Pati–Salam fermion content has even numbers of both left and right \(SU(2)\) doublets.

Thus no anomaly obstruction has been found.

## 5. Lorentz and strong-CP consistency

### Lorentz sector

The mediator terms are conventional Lorentz-scalar mass and mixing operators with dependence only on the transverse wall coordinate. They do not add higher spatial derivatives along the observable brane. Integrating out a heavy vectorlike field also produces Lorentz-covariant kinetic corrections suppressed by \(M_F^{-2}\).

The construction therefore preserves the earlier assumption that charged fermions have a standard Lorentz-invariant tree-level brane kinetic term, while photon-sector Lorentz violation arises from the gauge tower.

### Strong CP

At the real-amplitude stage, the integer mediator matrix, the family flavon and the Pati–Salam spurions can all be chosen real and PQ neutral. They then introduce no explicit Peccei–Quinn breaking.

A later complex spurion needed for the CKM phase may generate an order-one

\[
\arg\det(Y_uY_d).
\]

That does not invalidate the existing Route-II strong-CP architecture: its bulk axion is specifically intended to relax this effective phase. The future complex extension must still verify that the new mediator sector does not create an unsuppressed explicit PQ-breaking operator.

## Integer scan: does flavor select 21?

The frozen matrix family was tested for every integer

\[
2\leq N\leq80.
\]

For each \(N\), only the four already-allowed continuous amplitudes

\[
h_{d0},\quad h_{d1},\quad a_{d0},\quad a_{d1}
\]

were reoptimized.

| Integer | Maximum observable error |
|---:|---:|
| 20 | 1.946% |
| **21** | **0.596%** |
| 22 | 2.910% |
| 23 | 5.048% |

**N=21 is the only integer from 2 through 80 that remains below the preregistered 1% threshold.**

This is a sharp selection result inside the frozen model. It is not a blind discovery because the ratios were originally noticed from the successful coefficients before the scan.

## Overall verdict

The five ordered tests give:

1. **Mediator mass matrix:** strong success. One positive integer \(2\times2\) matrix generates all three ratios exactly.
2. **Flavor symmetry:** partial success. An anomaly-free \(U(1)_F\) or discrete remnant enforces \((-1,0,+1)\), but does not alone derive 21.
3. **Representation contraction:** strong conditional success. Pati–Salam gives the exact \((1,2,0)\) sector offsets and two physically relevant appearances of 21.
4. **Gauge/anomaly audit:** pass under the explicit vectorlike and complete-family assignments.
5. **Lorentz/strong-CP audit:** pass conditionally; no conflict appears if the real mediator sector is PQ neutral.

The result is best classified as

> **a concrete conditional UV derivation of the rational flavor relations, with one remaining dynamical matching problem: proving that the SFV/Pati–Salam mediator mass invariant is exactly \(N=21\).**

The four absolute amplitudes remain independent and must still be derived from mediator masses, wall normalization, or additional UV matching.

## Primary references

- J. C. Pati and A. Salam, “Lepton Number as the Fourth Color,” *Phys. Rev. D* **10**, 275 (1974), DOI 10.1103/PhysRevD.10.275.
- R. N. Mohapatra and J. C. Pati, “A Natural Left-Right Symmetry,” *Phys. Rev. D* **11**, 2558 (1975), DOI 10.1103/PhysRevD.11.2558.
- L. E. Ibáñez, “More About Discrete Gauge Anomalies,” arXiv:hep-ph/9210211.
