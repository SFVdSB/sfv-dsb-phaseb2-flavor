# Phase B2 v0.5.0: Raw-Gradient Wilson Closure

## Question

After freezing the successful local operator architecture, are the remaining order-one flavor coefficients genuinely independent, or did the numerical normalization of the wall-core invariant obscure a simpler matching law?

## Frozen operator architecture

The benchmark calculation uses

\[
B_{Ai}(y)=q_{Ai}O(y)+c_{Ai}G(y)+B_{\rm geo}(y),
\]

where

\[
G(y)=(\partial_y\Phi)^2+(\partial_y\phi)^2
\]

is the **raw local gradient density**, `O` is the chirality-producing kink, and `B_geo` is the separately calculated universal two-channel Hessian-rotation correction.

Earlier numerical work instead used

\[
I_G(y)=G(y)/G_{\max}.
\]

Consequently, the reported normalized coefficients were

\[
\kappa_{Ai}=c_{Ai}G_{\max}.
\]

As the wall stiffness changes, `Gmax` changes strongly. Treating `kappa` as the fundamental coefficient therefore made a mostly constant Wilson coefficient appear wall dependent.

## Exact benchmark in the raw basis

At the zero-bias benchmark,

\[
G_{\max}=0.7578850.
\]

The exact local-plus-geometric refit gives the raw coefficient coordinates

| Quantity | Exact value |
|---|---:|
| \(c_Q\) | 3.684647 |
| \(c_{u0}\) | 3.842283 |
| \(c_{u1}\) | -0.0416952 |
| \(a_{d0}\) | 0.259641 |
| \(a_{d1}\) | 0.242283 |
| \(c_{d0}\) | 3.515057 |
| \(c_{d1}\) | -0.869506 |

The refit reproduces the seven target observables to numerical precision.

## Frozen rational/geometric hypothesis

A simple post-hoc map was frozen:

\[
c_Q=\frac{11}{3},\qquad
c_{u0}=\frac{23}{6},\qquad
c_{u1}=-\frac1{24},
\]

\[
a_{d0}=\frac{m_{\Phi,T}}4,
\qquad
a_{d1}=\frac{\alpha}{2}(R_{\rm mix}-R_{\rm grad}),
\]

\[
c_{d0}=\frac72,
\qquad c_{d1}=-\frac78.
\]

At the benchmark, each formula lies within 1.51% of the exact raw coefficient; five of the seven lie within 0.69%.

These fractions were recognized after examining the solution. They are therefore candidate matching relations, not derived group factors or predictions.

## Cross-wall stability

The full local coefficients were independently refitted on all 51 available walls. Over the 33-point local design (baseline, axes, interactions, and held-out walls), the frozen formulas track the refitted raw coefficients with the following mean absolute errors:

| Relation | Mean error | Worst error |
|---|---:|---:|
| \(c_Q=11/3\) | 0.654% | 1.705% |
| \(c_{u0}=23/6\) | 1.071% | 2.860% |
| \(c_{u1}=-1/24\) | 0.647% | 1.928% |
| \(a_{d0}=m_{\Phi,T}/4\) | 1.476% | 1.968% |
| \(a_{d1}=L_{\rm lock}/2\) | 1.244% | 2.908% |
| \(c_{d0}=7/2\) | 0.843% | 1.745% |
| \(c_{d1}=-7/8\) | 0.680% | 1.455% |

This is much more stable than the normalized coefficients. It supports the interpretation that the raw-gradient Wilson basis is the physically appropriate coordinate system.

The relations are not universal across the full \(\pm8\%\) corridor, especially for large stiffness changes. The extreme-`Y` walls show nonlinear departures. This prevents the rational values from being promoted to exact laws.

## Flavor compression

Using all seven frozen formulas with no continuous flavor fit gives a maximum observable error of

\[
2.9416\%.
\]

The complete subset audit gives:

- four frozen relations and three fitted controls: **0.5304%** maximum error;
- five frozen relations and two fitted controls: **0.7021%** maximum error;
- six frozen relations and one fitted control: **1.3803%** maximum error.

The best five-fixed model freezes

\[
\left\{c_Q,c_{u0},c_{u1},a_{d1},c_{d1}\right\}
\]

and fits only

\[
\left\{a_{d0},c_{d0}\right\}.
\]

Its fitted values are

\[
a_{d0}=0.257738,
\qquad c_{d0}=3.514368,
\]

and every target remains within 0.703%.

Thus, within this frozen hypothesis, the residual continuous coefficient space has been reduced from seven to **two** at the one-percent strong-realization level.

This is an exploratory result because the rational formulas and the best subset were selected after studying the benchmark.

## Chiral spectral audit

Both the zero-fit rational model and the five-fixed/two-fitted model retain:

- 9 of 9 desired profiles with exactly one near-zero chiral state;
- 9 of 9 opposite-chirality partners with no near-zero state;
- minimum unwanted eigenvalue approximately 0.877.

The coefficient compression therefore does not spoil localization or chirality.

## Internal holdout diagnostic

With five relations frozen, the two remaining controls were calibrated to pairs of observables and the other five were withheld. Among the 12 balanced mass-plus-CKM pairs:

- 4 predicted every withheld observable within 1%;
- 9 were within 2%;
- the median worst withheld error was 1.40%;
- the worst balanced case was 4.26%.

This is not a blind prediction because the formulas and branch prior were developed using the complete benchmark. It does show that the two-control model has nontrivial predictive structure rather than behaving as a seven-target interpolation.

## Identifiability

The seven singular values of the local residual Jacobian are

\[
70.70,\ 29.16,\ 11.41,\ 3.60,\ 1.36,\ 0.340,\ 0.0273.
\]

All seven directions remain locally identifiable, but the condition number is about 2591. The weakest direction is mainly a correlated shift of the up- and down-sector core intercepts. This explains why two fitted controls can absorb small imperfections in several frozen relations.

## Interpretation

What is now strongly supported:

1. The local gradient operator is the natural core operator.
2. Its raw Wilson coefficients are much more stable than coefficients multiplying a normalized profile.
3. Several required coefficients are close to simple rational values and remain close across the local 33-wall design.
4. Five frozen relations plus only two continuous controls reproduce the full seven-observable benchmark below one percent.
5. The correct chiral spectrum survives.

What is not established:

1. The rational numbers have not been derived from a symmetry, mediator representation, or UV matching calculation.
2. The best subset was selected post hoc.
3. The formulas do not remain exact on extreme stiffness walls.
4. No genuinely unseen flavor observable has yet been predicted.

## Next logical step

The remaining task is now narrower. A minimal UV mediator/flavon sector should be constructed that can generate the raw coefficients

\[
\left(\frac{11}{3},\frac{23}{6},-\frac1{24},\frac72,-\frac78\right)
\]

or a nearby correlated set after integrating out one or two heavy fermion channels. The two residual quantities \(a_{d0}\) and \(c_{d0}\) provide the immediate matching targets. The operator content and representations must be frozen before evaluating a new observable.
