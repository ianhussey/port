# Sparse witnesses for rounded correlation tables

## Summary

`port` proves a rounded correlation table inconsistent by exhibiting one test vector whose box bound is negative. Near the boundary of the PSD cone, the directions it used to try (eigenvectors of the full matrix, coordinate pairs, triples, regression residuals) often find nothing, and the verdict stays `undecided`. This note reports a simulation of such cases and the search that came out of it.

1. **Sparse directions certify most of these cases.** The rounding penalty in the bound grows with the number of variables a direction involves. A greedy backward elimination over principal submatrices, with each round's best direction polished on its own support, certifies 198 of the 202 cases the previous search left undecided, at polynomial cost. It is implemented as `.greedy_subset_witness()` (Sections 1 to 3).
2. **A single vector was always enough in the simulation.** Enumeration shows that all 202 cases have a single-vector certificate. A second direction, taken from alternating projections, certifies the four cases that greedy elimination misses. It is implemented as `.pocs_witness()`, and `check_corr_psd()` now decides all 202. Both searches were chosen on these 202 cases. On 546 held-out cases with $p$ up to 20, six (1.1%) stay undecided (Section 4).
3. **Three block-based ideas did not pay off** (Section 5): a witness built from canonical correlations, a `block` localization class, and a block-shrinkage construction.

The block ideas came from reading Royen's proof of the Gaussian correlation inequality. The inequality itself gives `port` nothing; `gci.md` explains why.

---

## 1. The box bound and its penalty

Let $R$ be the reported matrix and $\delta$ the rounding half-width. The rounding box is the set of symmetric matrices $X$ with unit diagonal and $|X_{ij} - R_{ij}| \le \delta$ off the diagonal. For any $v \ne 0$,

$$
\max_{X \in \text{box}} v^\top X\, v \;=\; v^\top R\, v + \delta\big(\|v\|_1^2 - \|v\|_2^2\big),
$$

attained at $X_{ij} = R_{ij} + \delta \operatorname{sign}(v_i v_j)$. If this maximum is negative, no matrix in the box is PSD, and $v$ certifies the table as inconsistent. The bound is the mirror image of the inequality in Rohn's (1994) proof of Theorem 1.

Two properties of the bound matter for the search.

**Sparse directions pay less.** The penalty equals $2\delta \sum_{i<j} |v_i|\,|v_j|$. For a unit vector spread evenly over $k$ coordinates it is $\delta\,(k - 1)$. A direction supported on few variables therefore pays less of it. The price is a less negative $v^\top R\, v$, because the smallest eigenvalue of a principal submatrix $R_S$ is never below that of $R$. The best witness balances the two, and near the boundary it is usually supported on a subset of intermediate size.

**The best direction on a given support is an eigenvector of a shifted matrix.** For a fixed support $S$ and sign pattern $s$, the bound is the quadratic form of $R_S + \delta\,(s s^\top - I)$, and the best direction with that support and those signs is an eigenvector of *that* matrix whose signs match $s$, unless the best is approached only at the edge of the sign pattern, on a smaller support. (These are the vertex matrices $A_c - T_z \Delta T_z$ of Rohn, 1994, restricted to $S$ and with the sign of $\Delta$ reversed.) The bottom eigenvector of $R_S$ itself minimises $v^\top R\, v$ and ignores the penalty. `port`'s polishing iteration alternates between the sign pattern and the eigenvector, so it finds the better direction if it is restricted to $S$. Run on the full matrix, it turns a sparse direction dense again.

---

## 2. Simulation

### Design

The script `special-scripts/sparse-witness-simulation.R` (seed 2026; about 10 minutes, including the held-out run of Section 4 and the runs behind `forensic-value.qmd`) generates near-boundary rounded matrices:

1. a correlation matrix of rank 2 to $p - 2$ plus small diagonal noise,
2. one off-diagonal cell perturbed by up to ±0.08,
3. rounding to 2 decimals ($\delta = 0.005$).

A case counts as **undecided** if `port`'s original witness search finds no certificate and neither the constructions nor POCS certify consistency. On each undecided case the script tries the candidate families below. Every candidate goes through the rigorous box bound, so every count is a genuine certificate.

### Results

| | p = 5 | p = 6 | p = 8 | p = 10 |
|---|---|---|---|---|
| Matrices | 2436 | 2457 | 2473 | 1488 |
| Inconsistent (original search) | 466 | 483 | 507 | 282 |
| Consistent (construction or POCS) | 1960 | 1941 | 1880 | 1133 |
| **Undecided (original search)** | **10** | **33** | **86** | **73** |
| … certified by submatrix eigenvectors, sizes 3–4 | 7 | 26 | 48 | 20 |
| … certified by submatrix eigenvectors, all sizes | 7 | 30 | 76 | 67 |
| … certified by greedy elimination, unpolished | 7 | 30 | 75 | 66 |
| … certified by greedy elimination, winners polished on their support (implemented) | 10 | 33 | 84 | 71 |
| … certified by greedy elimination, every trial polished | 10 | 33 | 86 | 70 |
| … certified by the margin-free POCS direction, at most 1000 iterations | 9 | 32 | 83 | 70 |
| … certified by the margin-free POCS direction, at most 10000 iterations (implemented as a last resort, Section 4) | 10 | 33 | 86 | 73 |
| … certified by the best single vector (enumeration) | 10 | 33 | 86 | 73 |
| … certified by *s*-optimised R² direction (Section 5.1) | 0 | 2 | 0 | 0 |
| … certified by block witness, all bipartitions (Section 5.1) | 3 | 16 | 25 | 21 |
| … certified consistent by cross-block shrinkage (Section 5.3) | 0 | 0 | 0 | 0 |
| **… certified by `check_corr_psd()` now** | **10** | **33** | **86** | **73** |

### Interpretation

What matters is **sparsity**. `port` already used submatrices of size 3, but only those, and only for $p \le 12$. Restricting to sizes 3–4 loses much of the gain as $p$ grows, so the useful subsets are of intermediate size.

The right support is not enough, though. Unpolished, greedy elimination misses a tenth of the certificates it finds when each round's winner is polished on its support (178 against 198). This is the second property from Section 1.

The enumeration row gives the ceiling for any single test vector. Within one cell (a support and a sign pattern) a minimiser of the bound is an eigenvector of $R_S + \delta\,(s s^\top - I)$ whose signs match $s$, so enumerating all $(3^p - 1)/2$ cells finds the global optimum. **All 202 undecided cases have a single-vector certificate.** Every one of them was inconsistent; none was a consistent matrix that POCS failed to reach, and none needed more than one vector.

---

## 3. Implementation

Exhaustive search over subsets costs $O(2^p)$, and the enumeration above $O(3^p)$. `.greedy_subset_witness()` (`R/witness.R`) approximates them by backward elimination:

1. Start from all $p$ variables.
2. For each variable in the current set, form the bottom eigenvector of the submatrix without it and evaluate the rigorous box bound.
3. Polish the best of these directions on its own support.
4. Drop the variable whose removal gave the most negative bound; stop at the first certificate, or at size 3 (size 2 is covered by the coordinate pairs).

This takes $O(p^2)$ eigendecompositions of size up to $p - 1$, so $O(p^5)$ time: about a second at $p = 80$. It found 198 of the 202 certificates (98%). Polishing every direction in step 2, not only the best, found one more in total (199: two more at $p = 8$, one fewer at $p = 10$) at about twice the cost.

Because of that cost, `check_corr_psd()` runs the search only when the existing directions (eigenvectors, pairs, triples, residual directions, polishing) certify nothing *and* no construction (or, for heterogeneous boxes, POCS) certifies consistency. A matrix certified consistent has no witness, so nothing is lost, and matrices that a construction certifies cost what they did before. Verdicts that were already `inconsistent` keep their witnesses. Results that reach the search without a certificate can report a smaller `b_upper` than before. That has a side effect: where the smaller `b_upper` falls below `tau`, a table that used to be presumed consistent without a certificate now goes on to POCS and can leave it with one. Both the uniform path (`.witness_search()`) and the heterogeneous-box path (`.box_inconsistent()`, used by the localizers and the disattenuation checks) use the search.

The localizers' feasibility oracle, `.box_feasible()`, does not defer the search. Bisection asks about many boxes that are only just infeasible, and there a failed POCS run costs far more than the sparse search. With the search deferred until after POCS, `localize_psd_fault()` took 2.4 times as long over eight tables with $p$ from 8 to 30.

With greedy elimination in place, dropping the triples family lost no certificate in any run. The triples are kept anyway as a cheap first pass, since `.box_inconsistent()` sits inside that oracle.

**A corrected regression test.** `tests/testthat/test-escalation.R` contained a 5 × 5 matrix described as "genuinely boundary-only", expected to stay `undecided`. The new search certifies it as inconsistent with a witness supported on variables 1–4. An exact integer-arithmetic check confirms the certificate: the maximum of $v^\top X v$ over the rounding box is $-6.46 \times 10^{-5}\, \|v\|_2^2 < 0$. The original description was wrong; the old search simply could not find the sparse witness. The test now expects `inconsistent` and repeats the integer check.

A second matrix, first chosen to replace it as the `undecided` example, turned out to be inconsistent as well: its witness on variables 2–5 only appears after polishing on that support. It is now a test of exactly that. The `undecided` path is tested with a box that clips to the single point $r = 1$ in a 2 × 2 matrix. The only in-box matrix is singular, so there is no witness and nothing for the Cholesky test to certify, whatever the search does.

---

## 4. The last resort: a direction from alternating projections

Greedy elimination leaves 4 of the 8854 generated matrices undecided (none at $p \le 6$, two each at $p = 8$ and $p = 10$). All four are inconsistent, with a single-vector certificate that it did not reach.

A second direction closes that gap. Alternate projections between the box and the PSD cone without a margin, and take the bottom eigenvector of the box-side iterate. For an inconsistent box the iterates approach a closest pair $(X, Y)$, with $X$ in the box and $Y$ in the cone, and $Y - X$ separates the box from the cone. $Y - X$ is the negative eigen-part of $X$ with its sign flipped, so the bottom eigenvector of $X$ is its leading rank-one term.

Near the boundary the iteration converges slowly, so the number of iterations matters. Capped at 1000 iterations, this direction certified 194 of the 202 cases, fewer than greedy elimination, though the 194 include the four that greedy elimination misses. With up to 10000 iterations it certified all 202 on its own.

`.pocs_witness()` (`R/pocs.R`) implements it. It tests the direction every 100 iterations, stops at the first certificate, and gives up after 10000. That is up to 10000 eigendecompositions, about 17 seconds at $p = 100$ when no certificate turns up. `check_corr_psd()` therefore tries it only where the verdict would otherwise be `undecided`, and the localizers' oracle does not use it. A verdict reached this way has `tier = "pocs"`. With it, none of the 8854 generated matrices stays undecided.

**Held-out matrices.** That last result is in-sample: both searches were chosen on these matrices. The script therefore ends with a run on fresh matrices from the same generator (seed 2027, 2000 draws per $p$, larger $p$ included), on which nothing was tuned. "Open" has the meaning of "undecided" in Section 2: the original search finds no certificate and neither the constructions nor POCS certify consistency.

| | p = 8 | p = 10 | p = 15 | p = 20 |
|---|---|---|---|---|
| Open after the original search | 70 | 106 | 159 | 211 |
| … certified by greedy elimination | 70 | 103 | 152 | 200 |
| … certified by the POCS direction | 0 | 2 | 6 | 7 |
| … … of which after more than 1000 iterations | 0 | 0 | 2 | 1 |
| Presumed consistent (margin above `tau`, so the POCS direction is never tried) | 0 | 0 | 0 | 0 |
| **Still undecided** | **0** | **1** | **1** | **4** |

`check_corr_psd()` decides 540 of the 546 open cases (98.9%), 525 of them by greedy elimination. Six stay undecided, four of them at $p = 20$, so a remainder exists and it grows with $p$. The script does not establish whether these six are inconsistent; the enumeration of Section 2 is out of reach beyond $p \approx 10$.

Single vectors need not always suffice. In general, the box is inconsistent if and only if some PSD matrix $W \ne 0$ satisfies

$$
\max_{X \in \text{box}} \langle W, X \rangle = \langle W, R \rangle + \delta \sum_{i \ne j} |W_{ij}| < 0
$$

(separation of the compact box from the closed PSD cone; uniform unclipped box shown). Single vectors are the rank-one case $W = v v^\top$. A sum of rank-one terms can succeed where each term fails, because the off-diagonal absolute values can cancel. Finding the best $W$ is a semidefinite program. That would make the inconsistency test complete, at the cost of an SDP solver, which `port` deliberately avoids. The matrix $Y - X$ above is such a $W$, so alternating projections could supply matrix certificates without a solver, but `port`'s results hold a single witness vector. None of the 202 cases needed more.

---

## 5. Alternatives tried: block structure

Royen's proof of the Gaussian correlation inequality (`gci.md`) works with a covariance matrix split into two blocks. Three ideas for `port` came from that structure. They rest on Schur complements and canonical correlations (Hotelling, 1936), not on the inequality, and none of them is recommended.

### 5.1 A block witness

**Proposition (block PSD criterion).** Let $R$ be a symmetric matrix partitioned into blocks $R_{11}, R_{12}, R_{21} = R_{12}^\top, R_{22}$, with $R_{11}$ and $R_{22}$ positive definite. Define

$$
K = R_{11}^{-1/2} R_{12} R_{22}^{-1/2}, \qquad \rho_1 = \sigma_{\max}(K).
$$

$\rho_1$ is the largest canonical correlation between the two blocks. Then

$$
R \text{ is PSD} \iff \rho_1 \le 1.
$$

*Proof.*

1. Since $R_{11}$ is positive definite, $R$ is PSD if and only if its Schur complement $R_{22} - R_{21} R_{11}^{-1} R_{12}$ is PSD.
2. Multiplying on both sides by $R_{22}^{-1/2}$ preserves this property. Since $R_{22}^{-1/2} R_{21} R_{11}^{-1} R_{12} R_{22}^{-1/2} = K^\top K$, the condition becomes $I - K^\top K \succeq 0$.
3. The eigenvalues of $K^\top K$ are the squared singular values of $K$, so $I - K^\top K \succeq 0$ holds if and only if $\sigma_{\max}(K) \le 1$. $\square$

**The criterion does not depend on the partition.** The right-hand side of the equivalence does not mention the partition. So if $R$ is not PSD, *every* partition whose diagonal blocks are positive definite has $\rho_1 > 1$. At the reported point, "the canonical correlation between blocks A and B exceeds 1" says nothing about A and B in particular; it only restates that $R$ is not PSD.

**Witness vector.** When $\rho_1 > 1$, we can write down a vector that certifies non-PSDness.

1. Take unit singular vectors $u, w$ with $K w = \rho_1 u$ (so also $K^\top u = \rho_1 w$).
2. Set $a = R_{11}^{-1/2} u$, $b = R_{22}^{-1/2} w$, and $v = (a;\, -b)$.
3. Then
   - $a^\top R_{11} a = u^\top u = 1$,
   - $b^\top R_{22} b = w^\top w = 1$,
   - $a^\top R_{12} b = u^\top R_{11}^{-1/2} R_{12} R_{22}^{-1/2} w = u^\top K w = \rho_1$.
4. Therefore

   $$
   v^\top R\, v = a^\top R_{11} a + b^\top R_{22} b - 2\, a^\top R_{12} b = 2\,(1 - \rho_1).
   $$

If $\rho_1 > 1$, then $v^\top R v < 0$, and $v$ can be used as a test direction in the box bound. It is still a single vector, so it is a certificate of the same kind `port` already produces. Since every partition yields one, its value lies only in which *direction* it supplies. In the simulation the block witnesses of all bipartitions together certified 65 of the 202 undecided cases. They help only because they add directions, and sparse directions do much better.

**Special case: one variable versus the rest.** Let block 1 be a single variable $i$, so $R_{11} = 1$ and $R_{12} = r_i^\top$, the vector of its correlations with the other variables.

1. $K = r_i^\top R_{22}^{-1/2}$ is a row vector, so its only singular value satisfies

   $$
   \rho_1^2 = K K^\top = r_i^\top R_{22}^{-1} r_i = R_i^2,
   $$

   the squared multiple correlation of variable $i$ on the others.
2. The criterion $\rho_1 \le 1$ is therefore exactly van Tilburg and van Tilburg's (2023) $R_i^2 \le 1$ criterion, which `port` uses in its R² localizer.

The witness directions belong to the same family. `port` uses the regression-residual direction $v = (\beta_i;\, -1)$ with $\beta_i = R_{22}^{-1} r_i$, for which $v^\top R v = 1 - R_i^2$. The block witness above gives, after reordering, $v \propto (\beta_i;\, -\rho_1)$. More generally, every direction $v_s = (\beta_i;\, -s)$ with $s > 0$ gives

$$
v_s^\top R\, v_s = \beta_i^\top R_{22}\beta_i - 2s\, r_i^\top \beta_i + s^2 = R_i^2 - 2s R_i^2 + s^2 .
$$

- With $s = 1$ this is $1 - R_i^2$ (`port`'s choice).
- With $s = \rho_1 = R_i$ it is $2R_i^2(1 - R_i)$.

Both are negative exactly when $R_i > 1$, so they certify the same thing when rounding is ignored. The raw values cannot be compared to rank the directions, because the box bound is homogeneous of degree 2 in $v$ and $\|v_s\|$ depends on $s$. With rounding taken into account, the penalty depends on the direction, so a one-dimensional search over $s$ could in principle sharpen verdicts near the boundary. In the simulation it almost never mattered (2 of 202 undecided cases), so it has not been implemented.

**Precondition.** The Proposition requires $R_{11}$ and $R_{22}$ to be positive definite. This matches `port`'s note that clean blame via the R² localizer requires a positive definite complement.

### 5.2 A `block` localization class

`port`'s current verdict classes are `cell`, `cell_tentative`, `triad`, `variable`, `joint`, `diffuse`, and `none`. A `block` class would sit between `variable` and `joint` / `diffuse` and describe this situation:

> Instrument A's items are internally consistent, and instrument B's items are internally consistent, but the reported cross-correlations imply a canonical correlation above 1.

In box terms: the rounding boxes of $R_{11}$ and $R_{22}$ each contain a PSD matrix, but the full box does not.

**Why it adds little.**

1. *Detection is unaffected.* Whether the full box is inconsistent is decided by the witness search, whatever the partition.
2. *The canonical correlation is not needed.* Both premises can be checked by running `check_corr_psd()` on the two diagonal submatrices.
3. *The conclusion is close to a tautology.* If both diagonal blocks are consistent, setting every cross-block cell to 0 gives a PSD block-diagonal matrix. So freeing the cross-block cells always restores consistency, and "the fault lies in the cross-block cells" follows automatically. The existing localizer already tries to narrow the fault further, to cells, triads or variables.
4. *Partitions are costly.* There are $2^{p-1} - 1$ bipartitions, so they would have to come from the paper's design (two scales, predictors versus outcomes, two time points), adding an argument to `localize_psd_fault()`.

**What might be worth keeping.** As a *reporting* device, a sentence such as "the reported correlations imply a canonical correlation of 1.04 between scale A and scale B" is a vivid, effect-size-style description in the spirit of van Tilburg and van Tilburg's R² > 1. If added, it should be descriptive output for a user-supplied partition, not a verdict class. It is low priority.

### 5.3 Block shrinkage for certifying consistency

**Claim.** For $R(\tau) = \tau R + (1 - \tau)\operatorname{diag}(R_{11}, R_{22})$, the interpolation path of Royen's proof,

$$
\lambda_{\min}\big(R(\tau)\big) \;\ge\; \tau\, \lambda_{\min}(R) + (1 - \tau)\, \min\{\lambda_{\min}(R_{11}),\, \lambda_{\min}(R_{22})\}.
$$

*Proof.*

1. $\lambda_{\min}(A) = \min_{\|x\|_2 = 1} x^\top A x$ is a minimum of functions that are linear in $A$, so $\lambda_{\min}$ is concave on symmetric matrices:
   $\lambda_{\min}\big(\tau A + (1-\tau) B\big) \ge \tau\, \lambda_{\min}(A) + (1-\tau)\, \lambda_{\min}(B)$.
2. The eigenvalues of a block-diagonal matrix are those of its blocks, so $\lambda_{\min}\big(\operatorname{diag}(R_{11}, R_{22})\big) = \min\{\lambda_{\min}(R_{11}), \lambda_{\min}(R_{22})\}$. $\square$

**Why it does not support the proposed construction.** The proposal was to shrink only the cross-block cells toward zero, as far as their rounding boxes allow, and verify the result with Rump's (2006) one-sided Cholesky test.

1. Shrinking within the box subtracts up to $\delta$ from each cell's absolute value. $R(\tau)$ instead *scales* every cross-block cell by $\tau$. The two coincide only by accident, so the bound says nothing about the matrix actually constructed.
2. In the interesting case, $\lambda_{\min}(R) < 0$, the bound is weak exactly where it would be needed.
3. Behind the existing constructions, `port` already runs a POCS search for an in-box PSD matrix.

In the simulation, cross-block shrinkage over *all* bipartitions certified none of the 202 undecided cases.

---

## References

- Hotelling, H. (1936). Relations between two sets of variates. *Biometrika*, 28(3/4), 321–377.
- Hussey, I. (2026). *port: Positive-semi-definiteness Of Rounded correlation Tables* [Computer software]. https://github.com/ianhussey/port. doi:10.5281/zenodo.21434236
- Rohn, J. (1994). Positive definiteness and stability of interval matrices. *SIAM Journal on Matrix Analysis and Applications*, 15(1), 175–184.
- Royen, T. (2014). A simple proof of the Gaussian correlation conjecture extended to some multivariate gamma distributions. *Far East Journal of Theoretical Statistics*, 48(2), 139–145. arXiv:1408.1028
- Rump, S. M. (2006). Verification of positive definiteness. *BIT Numerical Mathematics*, 46(2), 433–452.
- van Tilburg, W. A. P., & van Tilburg, L. J. A. (2023). Inconsistent hypotheses and effect-size limits. *Advances in Methods and Practices in Psychological Science*, 6(4).
