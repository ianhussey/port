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

## API

| function | purpose |
|---|---|
| `check_corr_psd(R, decimals = 2, delta = NULL, tau = NULL)` | main entry point; returns a `corr_psd_check` |
| `certificate(x)` | extract the witness vector and margin |
| `check_corr_psd_batch(mats, ...)` | screen a list of matrices → tibble |
| `print(x)` | human-readable summary, with the certificate arithmetic when impossible |

## References

- N. J. Higham (2002), *Accuracy and Stability of Numerical Algorithms*, 2nd ed.,
  SIAM. (Inner-product and Cholesky error bounds.)
- S. M. Rump (2006), *Verification of positive definiteness*, BIT Numerical
  Mathematics 46:433–452.

## License

MIT © Ian Hussey
