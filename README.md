# psdness

> Certify whether a *rounded* correlation matrix can be positive semidefinite.

`psdness` is a research-integrity / forensic-metascience tool. Given a reported
correlation matrix whose off-diagonal entries have been rounded to a fixed
number of decimals, it decides whether **any** positive semidefinite (PSD)
matrix is consistent with the resulting *rounding box*. The payload is
certifying **impossibility**: when no PSD matrix fits the box, the reported
matrix cannot be a genuine correlation matrix, and the tool returns a
machine-checkable certificate.

## The problem

A reported correlation matrix `R` (`p × p`, symmetric, unit diagonal) has its
off-diagonals rounded to `d` decimals. Let `δ = 0.5·10⁻ᵈ` (or a user-supplied
absolute `δ`). This induces a **rounding box** of candidate true matrices `X`:

```
X_ii = 1                    (exactly)
X_ij ∈ [R_ij − δ, R_ij + δ] ∩ [−1, 1]     (i ≠ j)
```

**Question.** Does there exist a symmetric `X` in this box that is PSD?

| answer | verdict |
|---|---|
| no such `X` | **`impossible`** (a genuine correlation matrix cannot round to `R`) |
| some `X` exists | **`possible`** |
| can't tell at this precision | **`undecided`** |

## The method — verified witness-vector bound

For any test direction `v ≠ 0`, the most PSD-favourable in-box matrix along `v`
gives, *exactly*,

```
max_{X in box}  vᵀ X v  =  vᵀRv  +  δ·(‖v‖₁² − ‖v‖₂²).
```

(Clipping the box to `[−1, 1]` only shrinks it, so ignoring that clip keeps this
an **upper** bound — sound for impossibility.) If this maximum is negative, then
no in-box matrix has a nonnegative quadratic form along `v`, so **none is PSD**,
and `v` is the certificate:

```
vᵀRv + δ·(‖v‖₁² − ‖v‖₂²) < 0   ⟹   IMPOSSIBLE.
```

The algorithm takes `v` from the smallest eigenvector of `R` (plus cheap
coordinate-pair and 3×3-submatrix directions that connect to the exact minor
conditions), and reports the most impossibility-favourable one. `v` need not be
an exact eigenvector — it is only a *test direction*, so its inaccuracy never
invalidates the certificate.

### Tiers

1. **`precheck`** — if any reported `|R_ij| − δ > 1`, the entry is out of range
   even at its nearest box edge: immediately impossible. `O(p²)`, no
   eigendecomposition.
2. **`witness`** — the rigorous bound above (primary path, all `p`).
3. **`pocs`** — a self-contained escalation (alternating projections onto the
   box and the PSD cone) used only in the precision-limited ambiguous zone; see
   below.

## The rigor guarantee

> **If `verdict == "impossible"`, then no PSD matrix exists in the rounding box,
> subject to the stated floating-point error model.**

The verdict does **not** depend on any solver tolerance. Because the witness
bound `B = vᵀRv + δ·(‖v‖₁² − ‖v‖₂²)` is a short computation, its rounding error
is bounded *a priori* with standard backward/forward-error analysis (Higham,
*Accuracy and Stability of Numerical Algorithms*, 2002, Ch. 3), using the
`γₙ = n·u/(1 − n·u)` bounds with unit roundoff `u = .Machine$double.eps/2`. That
error bound is added as slack, so the reported `B_upper` rigorously exceeds the
true real-arithmetic maximum. Impossibility is declared only when
`B_upper < 0`. R offers no portable rounding-mode control, so this a-priori
slack approach — not directed rounding — is the intended design.

A corollary of soundness: because `B_upper` over-estimates the true maximum, a
rounding box that contains **any** PSD matrix can never be reported
`impossible`. In particular, a genuine PSD correlation matrix rounded to `d`
decimals always leaves the original inside the box, so it is never a false
positive. (The package's property test checks exactly this over thousands of
random PSD matrices.)

### What impossibility does — and does not — prove

The precise claim behind an `impossible` verdict is:

> **the reported values cannot be the complete-data, single-sample Pearson
> correlation matrix of any dataset, up to the stated rounding.**

The inference from there to a *reporting error* requires ruling out legitimate
generators of genuinely non-PSD reported matrices:

- **pairwise deletion** — each correlation computed on different available
  cases; the assembled matrix need not be PSD;
- **polychoric / tetrachoric estimation** — estimated cell-by-cell and famously
  non-positive-definite in ordinary use;
- **meta-analytically assembled matrices** — each cell pooled from a different
  set of studies (the standard MASEM situation);
- **upstream disattenuation** — some papers report reliability-corrected
  matrices, which may legitimately leave the PSD cone;
- **splicing** — cells transcribed from different samples or subgroups.

The `severity` / `excusable_delta` output is the discriminator: these mechanisms
typically produce *small* violations, so an impossibility that would require
corrections tens of times the rounding step survives all of these excuses,
while a near-boundary violation does not. The impossible-verdict print
enumerates these cautions, and `check_corr_psd()` asks the analyst to confirm
the matrix is Pearson, complete-case, single-sample, and uncorrected before
reading the flag forensically.

### One headline number: the excusable imprecision

`excusable_delta(R)` returns the smallest uniform reporting half-width at which
the box admits a PSD matrix — the verdict and the severity folded into one
statistic in precision units. The matrix is impossible at precision `δ` iff
`excusable_delta(R) > δ`, and a value like `0.4` reads as *"inconsistent with
any correlation matrix unless each entry were mis-stated by more than ±0.4 — no
conventional rounding could excuse it."* The same quantity appears in
`localize_psd_fault()`'s severity and print.

### Verdict certification

Every result carries a `certified` flag. `impossible` is always certified (a
witness vector). A `possible` with `certified = TRUE` exhibits a Rump-verified
in-box PSD matrix; `certified = FALSE` marks a *presumed* possible ("not shown
impossible" — the witness margin is comfortably positive but no certificate was
found). The print labels the two cases explicitly.

### Rounding rules, mixed precision, and missing cells

`check_corr_psd()` also supports:

- **asymmetric rounding rules** — `rounding = "truncate" / "floor" / "ceiling"`
  builds the correct one-sided box of width `10^(-d)`. Using the rule that
  actually produced the values is a soundness requirement: a mismatched
  symmetric box can exclude the true value and mislabel a valid matrix (this is
  demonstrated in validation Sim 2);
- **per-cell precision** — `decimals` (or `delta`) may be a `p × p` matrix for
  mixed-precision tables;
- **missing cells** — symmetric `NA` off-diagonals are freed to `[-1, 1]`, so a
  verdict holds for *every* value the missing entries could take; an
  `impossible` with missing cells is the stronger claim "impossible whatever
  the unreported values were".

## The possible side — verified in-box PSD matrix

Possibility is certified constructively: the tool exhibits an explicit in-box
matrix (the reported matrix itself; the box-edge matrix shrunk toward the
diagonal; the box-clipped PSD projection of `R`) and **verifies** it PSD with a
one-sided Cholesky test (Rump, *Verification of positive definiteness*, BIT
46:433–452, 2006): `A` is provably PSD if `chol(A − c·I)` completes, where `c`
rigorously bounds the Cholesky backward error. Any success is a rigorous proof
that a PSD matrix exists in the box.

When no cheap construction succeeds and the witness margin sits in the
precision-limited zone `[0, τ]` (default `τ = 10·δ`), the result escalates to a
self-contained **projections-onto-convex-sets (POCS)** search: alternating
Euclidean projections onto the box and onto the PSD cone (tightened by a small
margin `μ` for a strictly positive-definite result). For two closed convex sets
this converges to a point in their intersection whenever one exists; the
recovered in-box matrix is then independently Rump-verified before a `possible`
verdict is accepted. A failure to find a verifiable point is inconclusive
(`undecided`), never a rigorous impossibility. This tier needs no external
solver.

## Installation

```r
# install.packages("remotes")
remotes::install_github("ianhussey/psd-ness")
```

There are no hard third-party dependencies beyond `tibble` (for the batch
helper's output).

## Usage

```r
library(psdness)

# A strongly inconsistent 3x3: r12 = r13 = 0.9, r23 = -0.9.
R <- matrix(c(1,   0.9,  0.9,
              0.9, 1,   -0.9,
              0.9, -0.9, 1), 3, 3)

check_corr_psd(R, decimals = 2)
#> <corr_psd_check>
#>   verdict : IMPOSSIBLE  (no PSD matrix fits the rounding box)
#>   tier    : witness
#>   p       : 3      delta : 0.005
#>   margin  : B = -0.79  (B_upper = -0.79 < 0)
#>   witness : v = (0.5774, -0.5774, -0.5774)
#>   certificate: v'Rv + delta*(||v||_1^2 - ||v||_2^2) = -0.79 < 0

# Accessor for the certificate.
certificate(check_corr_psd(R, decimals = 2))

# A genuine correlation matrix rounded to 2dp is never impossible.
G <- round(cor(matrix(rnorm(400), ncol = 4)), 2)
check_corr_psd(G, decimals = 2)$verdict
#> [1] "possible"
```

### Screening a corpus

```r
mats <- list(good = diag(3), bad = R)
check_corr_psd_batch(mats, decimals = 2)
#> # A tibble: 2 x 8
#>   id    verdict    tier    margin b_upper delta     p message
#>   <chr> <chr>      <chr>    <dbl>   <dbl> <dbl> <int> <chr>
#> 1 good  possible   witness  1       1     0.005     3 the reported matrix ...
#> 2 bad   impossible witness -0.79   -0.79  0.005     3 NA
```

The batch helper logs how often the POCS escalation tier fired — it is expected
to be rare, since the cheap rigorous tiers settle most cases.

## Fault localization

When a matrix is impossible-given-rounding, `localize_psd_fault()` supports
semi-automated inference about **where** the non-PSDness comes from and **how
severe** it is. Non-PSD is a *global* property, so a unique culprit is generally
under-identified — the tool reports an honest localization class and never
manufactures a single culprit when the evidence only supports a set.

Everything is box-aware: every feasibility/interval computation respects the
rounding box, never the reported point values. Independent, box-aware components
feed a deterministic ruleset:

- **A. per-cell interval / sole culprit** — for each cell, the interval it could
  take given the others within rounding; a nonempty interval marks a sole-culprit
  candidate with a `required_edit` beyond rounding.
- **B. impossible triples** — triples whose 3×3 box-max determinant is provably
  negative, via a sound interval bound on `det = (1−b²)(1−c²) − (a−bc)²`. Pure
  arithmetic; exposed on its own as `impossible_triples()`.
- **C. leave-one-out** — variables whose removal restores feasibility.
- **D. sparse correction / severity** — the smallest set of cells whose joint
  correction restores feasibility, plus severity measures.
- **R² localizer** — per-variable squared multiple correlation
  (van Tilburg & van Tilburg 2023): variables that are over-determined by the
  others (`R² > 100%`). Implemented as the regression-residual **witness
  direction** `v = (M₋ᵢ⁻¹cᵢ ; −1)`, for which `vᵀMv = 1 − R²ᵢ`, so it reuses the
  box-witness bound exactly — no interval determinant, no verified inversion. A
  variable-level voter, not an arbiter; clean blame requires a PD complement.

The verdict is one of `cell`, `cell_tentative`, `triad`, `variable`, `joint`,
`diffuse`, or `none`.

### Causes: benign structural vs substantive violation

Non-PSDness is *often benign* — a composite plus its subscores, or a full set of
category dummies, is exactly (rank-)deficient by construction (van Tilburg &
van Tilburg 2023; Lorenzo-Seva & Ferrando 2021). For a forensic tool the report
must surface these before anything that reads as a data-integrity inference.
Every impossible result carries a `structural` diagnosis that separates a benign
near-boundary rank-deficiency from a substantive over-the-boundary violation
(using the severity relative to the rounding step) and names the likely
generator (ipsative / composite / localized) from the near-dependency direction.

```r
# A composite ~ its subscores, nudged just past the rounding boundary:
localize_psd_fault(Cn, decimals = 2)
#>   Verdict: DIFFUSE - Diffuse: no small explanation
#>   Structural: The matrix sits within rounding of a singular (boundary)
#>     correlation matrix. The near-dependency has one variable opposite a cluster
#>     of {1, 2, 3, 4}, consistent with a composite built from its subscores.
#>     Structural rank-deficiency of this kind is often a benign modelling
#>     artifact, not evidence of a data problem (van Tilburg & van Tilburg 2023;
#>     Lorenzo-Seva & Ferrando 2021).
#>   Severity: largest single edit needed 0.0117 (beyond rounding; delta = 0.005) ...
```

When a matrix is **possible**, the result reports a `plausibility` gradient — the
per-variable reported-point `R²` and its closest approach to the 100% ceiling —
so "passes but sits at 96%" (or "exceeds 100% at face value but survives
rounding") is visible without implying a verdict.

### Implied intervals and single-cell imputation

`implied_interval()` gives the interval a cell could occupy for the matrix to be
valid, holding the others at their reported points (a fast closed form,
van Tilburg & van Tilburg 2023) or letting them roam their boxes. Pointed at a
genuinely missing (`NA`) entry, it imputes the interval that entry must have
occupied.

**Same engine, no solver.** The localizer reuses the base package's primitives
throughout — POCS as a *search* (finding feasible points, bisecting for
intervals and severities) and the witness bound + Rump verification as the
*certificates*. Impossibility sub-claims (an empty per-cell interval; "removal
still impossible"; an impossible triple; an over-determined variable) are sound;
feasibility sub-claims exhibit a Rump-verified witnessing matrix. Pass
`verify = FALSE` to skip the certification and run search-based only. No CVXR or
other optimizer is required.

**Scope.** This is a *forensic, post-hoc, rounding-aware* tool: it interrogates
only the cells actually present in a given reported matrix. It deliberately does
**not** import van Tilburg & van Tilburg's design-time arguments — that a
correlation structure is impossible because of *unmeasured population* variables,
or the "add a variable and it becomes impossible" construction — which concern
populations and hypotheses rather than a specific reported matrix.

## Disattenuated correlations

Correlations corrected for measurement (un)reliability —
`D_ij = R_ij / sqrt(rho_i·rho_j)` — are *inflated*, so a valid observed matrix can
become an impossible construct matrix under the reliabilities the authors report.
`check_disattenuated_psd()` decides whether the disattenuated matrix is
**impossible** (no PSD matrix fits, or a corrected correlation exceeds 1),
**implausible** (a corrected correlation exceeds a `max_plausible_r` cutoff, e.g.
0.9, but stays ≤ 1), or **possible**, taking the rounding of both the correlations
and the reliabilities into account. When reliabilities are unreported it returns
the **critical reliability** `rho* = max(1 − lambda_min(R), max|R_ij|)` — how
reliable the measures must have been — with an optional plausibility floor for a
reachability statement. Alongside the point closed form, the tool reports a
**box-sound** critical reliability (`rho_impossible_box`, located by bisection
over the rigorous box witness): the certified statement "the disattenuation is
impossible for any common reliability below X *even allowing for rounding*",
which is the number quoted in the print and used for the headroom. The
impossibility side reuses the same witness/Rump engine, so it carries the same
floating-point guarantee. See the vignette for worked examples.

## API

| function | purpose |
|---|---|
| `check_corr_psd(R, decimals = 2, delta = NULL, tau = NULL, rounding = "nearest")` | main entry point; `decimals`/`delta` may be per-cell matrices, `NA` cells are freed, asymmetric rounding rules supported; returns a `corr_psd_check` with a `certified` flag |
| `certificate(x)` | extract the witness vector and margin |
| `excusable_delta(R)` | smallest reporting imprecision that could excuse the matrix (verdict + severity in precision units) |
| `check_corr_psd_batch(mats, ...)` | screen a list of matrices → tibble |
| `localize_psd_fault(x, verify = TRUE, sparse_k = 3, ...)` | localize the fault in an impossible matrix; returns a `psd_fault` |
| `fault_evidence(x)` | the raw A–D + R² evidence behind a `psd_fault` |
| `impossible_triples(R, decimals = 2, ...)` | standalone component B (rounding-robust impossible triples) → tibble |
| `implied_interval(R, cells, ...)` | interval a cell could occupy; imputes a missing (`NA`) entry → tibble |
| `localize_psd_fault_batch(mats, ...)` | localize over a list of matrices → tibble |
| `disattenuate(R, reliability)` | Spearman reliability disattenuation `D_ij = R_ij/sqrt(rho_i·rho_j)` |
| `check_disattenuated_psd(R, reliability, ..., max_plausible_r)` | is the disattenuated matrix impossible / implausible / possible, rounding-aware; or the critical reliability when unreported |

## Provenance of the methods

This package assembles established results from verified numerical computing; the
table records what was drawn from each source, and what is original here.

| Method (module) | Source | What was used |
|---|---|---|
| **Verified certificate of infeasibility.** The impossibility verdict is a rigorous certificate that the semidefinite feasibility problem "does a PSD matrix exist in the rounding box?" is *infeasible*, made independent of any solver tolerance by bounding all floating-point rounding *a priori* (no rounding-mode control). | Jansson, Chaykin & Keil (2008); Jansson (2009) — the verified-SDP paradigm | The **strategy**: post-process an approximate/heuristic direction into a rigorous infeasibility verdict via a-priori error bounds instead of trusting a solver's "infeasible" status. (That a dual direction certifies primal infeasibility is classical SDP duality; Jansson's contribution is the *verified* version.) |
| A-priori floating-point error bounds `γₙ = nu/(1−nu)` for the witness quadratic form (`R/rigor.R`, `R/witness.R`, `R/box.R`), and the Cholesky backward-error bound behind the PSD test | Higham (2002), *Accuracy and Stability of Numerical Algorithms* | The γₙ inner-product / quadratic-form error bounds and the Cholesky backward-error bound used to size the rigorous slack. |
| Verified positive semidefiniteness by a one-sided Cholesky of `A − cI` (`R/verify.R`) | Rump (2006) | The test itself and the constant `c` that bounds the Cholesky rounding error. |
| Possibility search / feasibility oracle by alternating projection onto the box and the PSD cone; nearest-correlation-matrix helper (`R/pocs.R`, `R/box.R`, `R/localize-sparse.R`) | Higham (2002, *IMA J. Numer. Anal.*); Cheney & Goldstein (1959) | Alternating projection onto the PSD cone and the unit-diagonal set (Higham's nearest-correlation-matrix iteration); convergence of projections onto two convex sets (Cheney–Goldstein). |
| Fault-localization layer: the R² / VIF variable localizer, the per-cell implied-interval closed form, single-cell imputation, the plausibility gradient, and the benign-vs-substantive causes taxonomy (`R/localize-rsquared.R`, `R/localize-causes.R`, `R/implied-interval.R`) | van Tilburg & van Tilburg (2023); Lorenzo-Seva & Ferrando (2021) | The squared-multiple-correlation impossibility criterion and its VIF / determinant identities (V1–V2); the closed-form per-cell interval (V3); single-unknown solving / imputation (V8); the possibility-vs-plausibility gradient (V4); and the taxonomy of benign structural causes to surface before any integrity inference (V5–V6). Realized on this package's box-aware witness/POCS engine rather than via matrix inversion or interval determinants (see below). |

**Original to this package** (derivations, not taken from the above): the
closed-form box quadratic-form maximum `v'Rv + δ(‖v‖₁²−‖v‖₂²)` and its use as a
box-aware witness; the witness-polishing fixed point (alternate the most
PSD-favourable in-box matrix along `v` with its bottom eigenvector — sound by
construction, since every iterate is re-evaluated through the rigorous bound);
the sound interval determinant bound
`det = (1−b²)(1−c²) − (a−bc)²` for the impossible-triple scan; the greedy
minimal-cardinality fault localization; the excusable-imprecision statistic;
the tiered checker and localization
ruleset; and the observation that van Tilburg's R²>1 criterion is the
regression-residual **witness vector** `v = (M₋ᵢ⁻¹cᵢ ; −1)` (since
`vᵀMv = 1 − R²ᵢ`), which is what lets the R² localizer reuse the existing
box-witness rigor model instead of an interval determinant or a verified
inverse. These residual directions also seed the main witness search.

## References

- W. Cheney and A. A. Goldstein (1959), "Proximity maps for convex sets,"
  *Proceedings of the American Mathematical Society*, 10(3):448–450.
- N. J. Higham (2002), "Computing the nearest correlation matrix—a problem from
  finance," *IMA Journal of Numerical Analysis*, 22(3):329–343.
- N. J. Higham (2002), *Accuracy and Stability of Numerical Algorithms*, 2nd ed.,
  SIAM. (Inner-product / quadratic-form and Cholesky error bounds.)
- C. Jansson, D. Chaykin, and C. Keil (2008), "Rigorous error bounds for the
  optimal value in semidefinite programming," *SIAM Journal on Numerical
  Analysis*, 46(1):180–200. [doi:10.1137/050622870](https://doi.org/10.1137/050622870)
- C. Jansson (2009), "On verified numerical computations in convex programming,"
  *Japan Journal of Industrial and Applied Mathematics*, 26(2–3):337–363.
- U. Lorenzo-Seva and P. J. Ferrando (2021), "Not positive definite correlation
  matrices in exploratory item factor analysis: causes, consequences and a
  proposed solution," *Structural Equation Modeling*, 28(1):138–147.
  [doi:10.1080/10705511.2020.1735393](https://doi.org/10.1080/10705511.2020.1735393)
- S. M. Rump (2006), "Verification of positive definiteness," *BIT Numerical
  Mathematics*, 46(2):433–452.
- W. A. P. van Tilburg and L. J. A. van Tilburg (2023), "Impossible hypotheses
  and effect-size limits," *Advances in Methods and Practices in Psychological
  Science*, 6(4).
  [doi:10.1177/25152459231197605](https://doi.org/10.1177/25152459231197605)

## License

MIT © Ian Hussey
