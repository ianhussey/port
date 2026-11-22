# The Gaussian correlation inequality, positive semidefiniteness, and `port`

## Summary

The Gaussian correlation inequality (GCI) and positive semidefiniteness (PSD) of correlation matrices are connected, but only in one direction. PSD is a *precondition* of the GCI and the main tool in Royen's (2014) proof. The GCI itself produces no new constraint on correlation matrices, so it offers no new test of the kind `port` performs.

Working through the proof brings out a block-level structure that looks useful for `port`. On closer examination, and in simulation (Section 6), the three proposals it suggests do not pay off:

1. A **block witness** built from canonical correlations (Section 3). It is a valid certificate, but at the reported point it fires for *every* partition, so it carries no partition-specific information. Its extra detection power comes from adding more test directions, not from the block structure, and a simpler family does much better.
2. A **`block` localization class** (Section 4). It is close to a tautology and adds little to the existing localizer.
3. A **block-shrinkage construction** (Section 5). The supporting bound does not describe the construction, and in simulation it certified no case the existing tiers missed.

The simulation did turn up a real improvement, which has been implemented: a **sparse witness search by greedy backward elimination over principal submatrices** (Section 6). It certifies most of the cases the previous witness search left undecided, at polynomial cost.

None of this depends on the GCI. The block material is classical multivariate analysis (Schur complements; Hotelling, 1936), which Royen's proof happens to put on display.

---

## 1. Where PSD enters Royen's proof

### Setting

Let $X \sim \mathcal{N}(0, C)$ be an $n$-dimensional centred Gaussian vector, with coordinates split into two blocks of sizes $n_1 + n_2 = n$:

$$
C = \begin{pmatrix} C_{11} & C_{12} \\ C_{21} & C_{22} \end{pmatrix}.
$$

In the form Royen proves, the GCI states that for any thresholds $t_i > 0$,

$$
P\big(|X_i| \le t_i \;\; \forall i \le n\big) \;\ge\; P\big(|X_i| \le t_i \;\; \forall i \le n_1\big)\; P\big(|X_i| \le t_i \;\; \forall\, n_1 < i \le n\big).
$$

The general statement for symmetric convex sets $K, L$, $\mu(K \cap L) \ge \mu(K)\,\mu(L)$, follows by approximation: every closed symmetric convex set is a countable intersection of symmetric strips $\{x : |\langle x, v\rangle| \le t\}$ (Latała & Matlak, 2017; Panchenko, 2020). Without loss of generality, $C$ is positive definite.

### Step 1: Interpolation

Define, for $\tau \in [0, 1]$,

$$
C(\tau) = \begin{pmatrix} C_{11} & \tau\, C_{12} \\ \tau\, C_{21} & C_{22} \end{pmatrix}.
$$

At $\tau = 0$ the two blocks are independent, so the probability on the left side equals the product on the right. At $\tau = 1$ it equals the joint probability. The proof shows that the probability is **nondecreasing in $\tau$**, which yields the inequality.

For this, $C(\tau)$ must be a valid covariance matrix for every $\tau$. This follows from convexity of the PSD cone:

$$
C(\tau) = \tau\, C + (1 - \tau)\, C_0, \qquad C_0 := \operatorname{diag}(C_{11}, C_{22}).
$$

For any $v = (v_1, v_2)$,

$$
v^\top C(\tau)\, v
= \tau\, v^\top C\, v + (1 - \tau)\big(v_1^\top C_{11} v_1 + v_2^\top C_{22} v_2\big)
\;\ge\; 0.
$$

Each term is nonnegative:

- $v^\top C v \ge 0$ because $C$ is PSD.
- $v_1^\top C_{11} v_1 = (v_1, 0)^\top C\, (v_1, 0) \ge 0$, and likewise for $C_{22}$. Every principal submatrix of a PSD matrix is PSD.

The same argument with $C$ positive definite shows that $C(\tau)$ is positive definite, uniformly in $\tau$.

### Step 2: Laplace transform of the squares

Royen works with the squared coordinates $X_i^2$, whose joint law is a multivariate gamma (chi-square) distribution. For $\lambda_i \ge 0$ and $\Lambda = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$, Lemma 4 of Latała and Matlak (2017) gives

$$
\mathbb{E}\exp\!\left(-\sum_{i=1}^n \lambda_i X_i^2\right) = \big|I_n + 2\Lambda C\big|^{-1/2}.
$$

The determinant can be expanded entirely in terms of **principal minors**:

$$
\big|I_n + 2\Lambda C(\tau)\big| = \sum_{J \subseteq \{1, \ldots, n\}} \left(\prod_{j \in J} 2\lambda_j\right) \big|C(\tau)_J\big|,
\qquad |C(\tau)_\emptyset| := 1.
$$

Here $A_J$ denotes the principal submatrix of $A$ with rows and columns in $J$.

*Justification.* The determinant is multilinear in the rows. Row $j$ of $I + DA$, with $D$ diagonal, is $e_j^\top + d_j\, a_j^\top$. Expanding every row into its two summands, the term that takes the $d_j a_j^\top$ part exactly for $j \in J$ equals $\prod_{j \in J} d_j$ times the determinant of $A_J$.

### Step 3: Monotonicity of the principal minors

There are two cases.

**Case 1: $J$ lies within one block.** Then $C(\tau)_J = C_J$, which does not depend on $\tau$.

**Case 2: $J$ meets both blocks.** Write $J = J_1 \cup J_2$, with $J_1$ in block 1 and $J_2$ in block 2.

1. The Schur-complement factorization of the determinant gives

   $$
   \big|C(\tau)_J\big| = |C_{J_1}| \cdot \big|C_{J_2} - \tau^2\, C_{J_2 J_1} C_{J_1}^{-1} C_{J_1 J_2}\big|.
   $$

2. Factor $C_{J_2}^{1/2}$ out on both sides of the second determinant:

   $$
   \big|C(\tau)_J\big| = |C_{J_1}|\,|C_{J_2}| \cdot \big|I - \tau^2 M\big|,
   \qquad M := C_{J_2}^{-1/2} C_{J_2 J_1} C_{J_1}^{-1} C_{J_1 J_2} C_{J_2}^{-1/2}.
   $$

3. $M = B^\top B$ with $B = C_{J_1}^{-1/2} C_{J_1 J_2} C_{J_2}^{-1/2}$, so $M$ is PSD. Its eigenvalues $\rho_k^2$ are the squared **canonical correlations** between the variable groups $J_1$ and $J_2$ (Hotelling, 1936). Because $C$ is PSD, $\rho_k \le 1$ (see the Proposition in Section 3). Hence

   $$
   \big|C(\tau)_J\big| = |C_{J_1}|\,|C_{J_2}| \prod_k \big(1 - \tau^2 \rho_k^2\big).
   $$

4. Every factor $1 - \tau^2 \rho_k^2$ lies in $[0, 1]$ and is nonincreasing in $\tau \in [0, 1]$. So $|C(\tau)_J|$ is **nonincreasing in $\tau$**, that is, $a_J(\tau) := -\partial_\tau |C(\tau)_J| \ge 0$.

**Interpretation.** Royen's interpolation turns down every canonical correlation between the two blocks by the factor $\tau$. At $\tau = 0$ they are all zero (independent blocks); at $\tau = 1$ they take their actual values.

### Step 4: From sign facts to probabilities

This step is where the multivariate gamma structure, rather than PSD alone, does the work. Write $f(\cdot, \tau)$ for the density of $Z(\tau) = \tfrac12 X(\tau)^2$ (coordinatewise).

1. Differentiating the Laplace transform $|I + \Lambda C(\tau)|^{-1/2}$ in $\tau$ gives

   $$
   \tfrac12 \sum_{J \ne \emptyset} a_J(\tau)\, \big|I + \Lambda C(\tau)\big|^{-3/2} \prod_{j \in J} \lambda_j .
   $$

2. $|I + \Lambda C(\tau)|^{-3/2}$ is itself the Laplace transform of a **different** density $h$: the multivariate gamma density with shape $3/2$ and the same $C(\tau)$. Multiplying a Laplace transform by $\prod_{j \in J} \lambda_j$ corresponds to taking the mixed partial derivative $\partial^{|J|} h / \partial x_J$ (Latała & Matlak, 2017, Lemmas 7–8). Hence

   $$
   \frac{\partial f}{\partial \tau}(x, \tau) = \tfrac12 \sum_{J \ne \emptyset} a_J(\tau)\, \frac{\partial^{|J|} h}{\partial x_J}(x).
   $$

3. Integrating over the rectangle $\prod_i [0, s_i]$ turns each mixed partial into an integral of $h$ over a face of the rectangle, because the boundary terms at $x_i = 0$ vanish. These integrals are nonnegative, and so are the $a_J(\tau)$, so the probability is nondecreasing in $\tau$.

### What PSD does in the proof

PSD has three jobs:

1. It guarantees that the Gaussian $\mathcal{N}(0, C)$ exists.
2. It keeps the interpolation path $C(\tau)$ inside the PSD cone.
3. It makes the principal minors nonnegative and monotone along that path (via Schur complements and $\rho_k \le 1$).

---

## 2. Why the GCI gives no new forensic test

### Direction of implication

The logical chain is:

$$
C \text{ is PSD} \;\Longrightarrow\; \mathcal{N}(0, C) \text{ exists} \;\Longrightarrow\; \text{the GCI holds for } \mathcal{N}(0, C).
$$

The contrapositive, "if the GCI fails, then $C$ is not PSD," cannot be used. For a non-PSD reported matrix $R$, the distribution $\mathcal{N}(0, R)$ does not exist. There are therefore no probabilities that could fail the inequality. The GCI cannot detect non-PSD matrices; it presupposes PSD.

### Exact identity versus distributional property

`port` works because PSD is an **exact algebraic identity** for the complete-case Pearson correlation matrix of *every* finite dataset, whatever its distribution.

1. Let $Z$ be the column-centred $n \times p$ data matrix. The sample covariance matrix is $S = Z^\top Z / (n - 1)$.
2. The sample correlation matrix is $R = D^{-1/2} S D^{-1/2}$, with $D = \operatorname{diag}(S)$.
3. For any $v \in \mathbb{R}^p$, set $w = D^{-1/2} v$. Then

   $$
   v^\top R\, v = w^\top S\, w = \frac{\|Z w\|_2^2}{n - 1} \;\ge\; 0.
   $$

This requires no distributional assumption. It does require that every correlation was computed on the same cases. Pairwise deletion, polychoric estimation, meta-analytic pooling, or disattenuation can each break it, which is why `port` asks the analyst to confirm these before reading an inconsistency forensically.

The GCI, by contrast, is a statement about **population Gaussian measures**. A finite sample's relative frequencies can violate it by chance, since the empirical distribution is not Gaussian.

### A hypothetical GCI-based check

Suppose a paper reported the proportions of participants within ±k SD of the mean on several variables, both separately and jointly. Under normality, the joint proportion should be at least the product of the marginal ones. Checking this would be a **probabilistic test resting on a normality assumption**, not a deterministic consistency check in the GRIM / `port` sense.

Moreover, for intervals taken coordinate by coordinate, the inequality

$$
P\big(|X_i| \le t_i \;\forall i\big) \;\ge\; \prod_i P\big(|X_i| \le t_i\big)
$$

is already the **Khatri–Šidák inequality** (Khatri, 1967; Šidák, 1967). Royen's contribution is the version where *each side is a block* of several coordinates, or more generally any pair of symmetric convex sets.

---

## 3. A block witness

### Proposition (block PSD criterion)

Let $R$ be a symmetric matrix partitioned into blocks $R_{11}, R_{12}, R_{21} = R_{12}^\top, R_{22}$, with $R_{11}$ and $R_{22}$ positive definite. Define

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

**Consequence: the criterion does not depend on the partition.** The right-hand side of the equivalence does not mention the partition. So if $R$ is not PSD, *every* partition whose diagonal blocks are positive definite has $\rho_1 > 1$. At the reported point, "the canonical correlation between blocks A and B exceeds 1" says nothing about A and B in particular; it only restates that $R$ is not PSD.

### Witness vector

When $\rho_1 > 1$, we can write down a vector that certifies non-PSDness.

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

If $\rho_1 > 1$, then $v^\top R v < 0$. This $v$ can be used directly as a test direction in `port`'s rigorous box bound, which flags the matrix as inconsistent when

$$
v^\top R\, v + \delta\big(\|v\|_1^2 - \|v\|_2^2\big) < 0.
$$

It is still a single vector, so it is a certificate of the same kind `port` already produces. Since every partition yields one, its value lies only in which *direction* it supplies to the box bound. Section 6 shows that this is worth something, but much less than sparse directions are.

### Special case: one variable versus the rest

Let block 1 be a single variable $i$, so $R_{11} = 1$ and $R_{12} = r_i^\top$, the vector of its correlations with the other variables.

1. $K = r_i^\top R_{22}^{-1/2}$ is a row vector, so its only singular value satisfies

   $$
   \rho_1^2 = K K^\top = r_i^\top R_{22}^{-1} r_i = R_i^2,
   $$

   the squared multiple correlation of variable $i$ on the others.
2. The criterion $\rho_1 \le 1$ is therefore exactly van Tilburg and van Tilburg's (2023) $R_i^2 \le 1$ criterion, which `port` uses in its R² localizer.

**The witness directions belong to the same family.** `port` uses the regression-residual direction $v = (\beta_i;\, -1)$ with $\beta_i = R_{22}^{-1} r_i$, for which $v^\top R v = 1 - R_i^2$. The block witness above gives, after reordering, $v \propto (\beta_i;\, -\rho_1)$.

More generally, every direction $v_s = (\beta_i;\, -s)$ with $s > 0$ gives

$$
v_s^\top R\, v_s = \beta_i^\top R_{22}\beta_i - 2s\, r_i^\top \beta_i + s^2 = R_i^2 - 2s R_i^2 + s^2 .
$$

- With $s = 1$ this is $1 - R_i^2$ (`port`'s choice).
- With $s = \rho_1 = R_i$ it is $2R_i^2(1 - R_i)$.

Both are negative exactly when $R_i > 1$, so they certify the same thing when rounding is ignored. The raw values cannot be compared to rank the directions, because the box bound is homogeneous of degree 2 in $v$ and $\|v_s\|$ depends on $s$. With rounding taken into account, the penalty $\delta\big(\|v\|_1^2 - \|v\|_2^2\big)$ depends on the direction, so a one-dimensional search over $s$ could in principle sharpen verdicts near the boundary. In simulation it almost never mattered (2 of 202 undecided cases; Section 6), so it has not been implemented.

### Precondition

The Proposition requires $R_{11}$ and $R_{22}$ to be positive definite. This matches `port`'s note that clean blame via the R² localizer requires a positive definite complement.

---

## 4. A `block` localization class (not recommended)

`port`'s current verdict classes are `cell`, `cell_tentative`, `triad`, `variable`, `joint`, `diffuse`, and `none`. A `block` class would sit between `variable` and `joint` / `diffuse` and describe this situation:

> Instrument A's items are internally consistent, and instrument B's items are internally consistent, but the reported cross-correlations imply a canonical correlation above 1.

In box terms: the rounding boxes of $R_{11}$ and $R_{22}$ each contain a PSD matrix, but the full box does not.

**Why it adds little.**

1. *Detection is unaffected.* Whether the full box is inconsistent is decided by the witness search, whatever the partition.
2. *The canonical correlation is not needed.* Both premises can be checked by running `check_corr_psd()` on the two diagonal submatrices.
3. *The conclusion is close to a tautology.* If both diagonal blocks are consistent, setting every cross-block cell to 0 gives a PSD block-diagonal matrix. So freeing the cross-block cells always restores consistency, and "the fault lies in the cross-block cells" follows automatically. The existing localizer already tries to narrow the fault further, to cells, triads or variables.
4. *Partitions are costly.* There are $2^{p-1} - 1$ bipartitions, so they would have to come from the paper's design (two scales, predictors versus outcomes, two time points), adding an argument to `localize_psd_fault()`.

**What might be worth keeping.** As a *reporting* device, a sentence such as "the reported correlations imply a canonical correlation of 1.04 between scale A and scale B" is a vivid, effect-size-style description in the spirit of van Tilburg and van Tilburg's R² > 1. If added, it should be descriptive output for a user-supplied partition, not a verdict class. It is low priority.

---

## 5. Block shrinkage for certifying consistency (not recommended)

**Claim.** For $R(\tau) = \tau R + (1 - \tau)\operatorname{diag}(R_{11}, R_{22})$,

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

In simulation, cross-block shrinkage over *all* bipartitions certified none of the 202 undecided cases (Section 6).

---

## 6. Simulation, and what was implemented

### Design

The script `special-scripts/gci-psd-simulation.R` (seed 2026, about 5 minutes) generates near-boundary rounded matrices:

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
| … certified by *s*-optimised R² direction | 0 | 2 | 0 | 0 |
| … certified by block witness, all bipartitions | 3 | 16 | 25 | – |
| … certified consistent by cross-block shrinkage | 0 | 0 | 0 | 0 |
| … certified by submatrix eigenvectors, sizes 3–4 | 7 | 26 | 48 | 20 |
| … certified by submatrix eigenvectors, all sizes | 7 | 30 | 76 | 67 |
| … certified by greedy elimination, unpolished | 7 | 30 | 75 | 66 |
| … certified by greedy elimination, winners polished on their support (implemented) | 10 | 33 | 84 | 71 |
| … certified by the best single vector (enumeration) | 10 | 33 | 86 | 73 |

The block witness was not run at $p = 10$ (511 bipartitions per case).

### Interpretation

The block witness helps only because it adds directions. What matters is **sparsity**. The rounding penalty $\delta\big(\|v\|_1^2 - \|v\|_2^2\big)$ grows with the spread of $v$ across coordinates, so directions supported on few variables pay less of it. The bottom eigenvectors of principal submatrices are such directions. `port` already used submatrices of size 3, but only those, and only for $p \le 12$. Restricting to sizes 3–4 loses much of the gain as $p$ grows, so the useful subsets are of intermediate size.

The right support is not enough, though. The bottom eigenvector of $R_S$ minimises $v^\top R v$ and ignores the penalty. For a fixed support $S$ and sign pattern $s$, the box bound is the quadratic form of $R_S + \delta\,(s s^\top - I)$, and the best direction is the bottom eigenvector of *that* matrix. (These are the vertex matrices $A_c - T_z \Delta T_z$ of Rohn, 1994, restricted to $S$ and with the sign of $\Delta$ reversed; the box bound itself is the mirror image of the inequality in his proof of Theorem 1.) The polishing iteration finds it, provided it is restricted to $S$. Run on the full matrix, as `port` did before, it turns a sparse direction dense again and the gain is lost. Unpolished, greedy elimination misses a tenth of the certificates it finds with polishing (178 against 198).

The last row gives the ceiling for any single test vector. Within one cell (a support and a sign pattern) a minimiser of the bound is an eigenvector of the matrix above whose signs match $s$, so enumerating all $(3^p - 1)/2$ cells finds the global optimum. **All 202 undecided cases have a single-vector certificate.** Every one of them was inconsistent; none was a consistent matrix that POCS failed to reach, and none needed more than one vector.

### Implementation

Exhaustive search over subsets costs $O(2^p)$, and the enumeration above $O(3^p)$. `.greedy_subset_witness()` (`R/witness.R`) approximates them by backward elimination:

1. Start from all $p$ variables.
2. For each variable in the current set, form the bottom eigenvector of the submatrix without it and evaluate the rigorous box bound.
3. Polish the best of these directions on its own support.
4. Drop the variable whose removal gave the most negative bound; stop at the first certificate, or at size 3 (size 2 is covered by the coordinate pairs).

This takes $O(p^2)$ eigendecompositions of size up to $p - 1$, so $O(p^5)$ time: about a second at $p = 80$. It found 198 of the 202 certificates (98%). Polishing every direction in step 2, not only the best, found one more at three times the cost.

Because of that cost, `check_corr_psd()` runs the search only when the existing directions (eigenvectors, pairs, triples, residual directions, polishing) certify nothing *and* no construction (or, for heterogeneous boxes, POCS) certifies consistency. A matrix certified consistent has no witness, so nothing is lost, and ordinary consistent matrices cost what they did before. Verdicts that were already `inconsistent` keep their witnesses. Results that reach the search without a certificate can report a smaller `b_upper` than before. Both the uniform path (`.witness_search()`) and the heterogeneous-box path (`.box_inconsistent()`, used by the localizers and the disattenuation checks) use the search.

With greedy elimination in place, dropping the triples family lost no certificate in any run. The triples are kept anyway as a cheap first pass, since `.box_inconsistent()` sits inside the localizers' feasibility oracle.

**A corrected regression test.** `tests/testthat/test-escalation.R` contained a 5 × 5 matrix described as "genuinely boundary-only", expected to stay `undecided`. The new search certifies it as inconsistent with a witness supported on variables 1–4. An exact integer-arithmetic check confirms the certificate: the maximum of $v^\top X v$ over the rounding box is $-2.27 \times 10^{-5}\, \|v\|_2^2 < 0$. The original description was wrong; the old search simply could not find the sparse witness. The test now expects `inconsistent` and repeats the integer check.

A second matrix, first chosen to replace it as the `undecided` example, turned out to be inconsistent as well: its witness on variables 2–5 only appears after polishing on that support. It is now a test of exactly that. The `undecided` path is tested with a box that clips to the single point $r = 1$ in a 2 × 2 matrix. The only in-box matrix is singular, so there is no witness and nothing for the Cholesky test to certify, whatever the search does.

### What remains undecided

After the change, 4 of the 8854 generated matrices stayed undecided (none at $p \le 6$, two each at $p = 8$ and $p = 10$). All four are inconsistent with a single-vector certificate that greedy elimination did not reach, so in this simulation the remaining gap is one of search, not of principle. One cheap direction that was not implemented: alternate projections between the box and the PSD cone without a margin, and take the bottom eigenvector of the box-side limit point. In a side experiment it certified 201 of the 202 cases.

Single vectors need not always suffice. In general, the box is inconsistent if and only if some PSD matrix $W \ne 0$ satisfies

$$
\max_{X \in \text{box}} \langle W, X \rangle = \langle W, R \rangle + \delta \sum_{i \ne j} |W_{ij}| < 0
$$

(separation of the compact box from the closed PSD cone; uniform unclipped box shown). Single vectors are the rank-one case $W = v v^\top$. A sum of rank-one terms can succeed where each term fails, because the off-diagonal absolute values can cancel. Finding the best $W$ is a semidefinite program. That would make the inconsistency test complete, at the cost of an SDP solver, which `port` deliberately avoids. No case in the simulation needed it.

---

## References

- Hotelling, H. (1936). Relations between two sets of variates. *Biometrika*, 28(3/4), 321–377.
- Hussey, I. (2026). *port: Positive-semi-definiteness Of Rounded correlation Tables* [Computer software]. https://github.com/ianhussey/port. doi:10.5281/zenodo.21434236
- Khatri, C. G. (1967). On certain inequalities for normal distributions and their applications to simultaneous confidence bounds. *Annals of Mathematical Statistics*, 38(6), 1853–1867.
- Latała, R., & Matlak, D. (2017). Royen's proof of the Gaussian correlation inequality. In B. Klartag & E. Milman (Eds.), *Geometric Aspects of Functional Analysis* (Lecture Notes in Mathematics, Vol. 2169, pp. 265–275). Springer. Preprint: arXiv:1512.08776
- Panchenko, D. (2020). *Gaussian correlation inequality* [Lecture notes]. https://sites.stat.washington.edu/people/mdperlma/542%20Multivariate%20Analysis%20Material/Panchenko%20re%20Royen.pdf
- Rohn, J. (1994). Positive definiteness and stability of interval matrices. *SIAM Journal on Matrix Analysis and Applications*, 15(1), 175–184.
- Royen, T. (2014). A simple proof of the Gaussian correlation conjecture extended to some multivariate gamma distributions. *Far East Journal of Theoretical Statistics*, 48(2), 139–145. arXiv:1408.1028
- Rump, S. M. (2006). Verification of positive definiteness. *BIT Numerical Mathematics*, 46(2), 433–452.
- Šidák, Z. (1967). Rectangular confidence regions for the means of multivariate normal distributions. *Journal of the American Statistical Association*, 62(318), 626–633.
- van Tilburg, W. A. P., & van Tilburg, L. J. A. (2023). Inconsistent hypotheses and effect-size limits. *Advances in Methods and Practices in Psychological Science*, 6(4).
