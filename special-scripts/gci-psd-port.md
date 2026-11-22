# The Gaussian correlation inequality, positive semidefiniteness, and `port`

## Summary

The Gaussian correlation inequality (GCI) and positive semidefiniteness (PSD) of correlation matrices are connected, but only in one direction. PSD is a *precondition* of the GCI and the main tool in Royen's (2014) proof. The GCI itself produces no new constraint on correlation matrices, so it offers no new test of the kind `port` performs.

Working through the proof does bring out a block-level structure that is useful for `port`:

1. A **block witness** built from canonical correlations. Its single-variable special case gives the same criterion as `port`'s R² localizer.
2. A possible **`block` localization class**.
3. A cheap **block-shrinkage construction** for certifying consistency.

None of these depend on the GCI. They are classical multivariate analysis (Schur complements; Hotelling, 1936), which Royen's proof happens to put on display.

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

The general statement for symmetric convex sets $K, L$, $\mu(K \cap L) \ge \mu(K)\,\mu(L)$, follows by approximation: every closed symmetric convex set is a countable intersection of symmetric strips $\{x : |\langle x, v\rangle| \le t\}$ (Latała & Matlak, 2017). Without loss of generality, $C$ is positive definite.

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

4. Every factor $1 - \tau^2 \rho_k^2$ lies in $[0, 1]$ and is nonincreasing in $\tau \in [0, 1]$. So $|C(\tau)_J|$ is **nonincreasing in $\tau$**.

**Interpretation.** Royen's interpolation turns down every canonical correlation between the two blocks by the factor $\tau$. At $\tau = 0$ they are all zero (independent blocks); at $\tau = 1$ they take their actual values.

Latała and Matlak (2017) use the same $B^\top B \succeq 0$ structure. The remaining steps of the proof turn these sign facts into a formula for $\partial f / \partial \tau$, the $\tau$-derivative of the multivariate gamma density. That formula is a sum of mixed partial derivatives $\partial^{|J|} / \partial x_J$ with nonnegative coefficients $a_J(\tau) \ge 0$. Integrating over the rectangle then shows that the probability is nondecreasing in $\tau$.

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

is already **Šidák's (1967) inequality**. Royen's contribution is the version where *each side is a block* of several coordinates, or more generally any pair of symmetric convex sets.

---

## 3. Where the overlap is useful for `port`: a block witness

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

Both are negative exactly when $R_i > 1$, so they certify the same thing when rounding is ignored. With rounding taken into account, the penalty $\delta\big(\|v\|_1^2 - \|v\|_2^2\big)$ depends on the direction. Choosing $s$ to make the bound as negative as possible could therefore give slightly sharper verdicts near the boundary. This is a cheap one-dimensional optimization.

### Precondition

The Proposition requires $R_{11}$ and $R_{22}$ to be positive definite. This matches `port`'s note that clean blame via the R² localizer requires a positive definite complement.

---

## 4. A proposed `block` localization class

`port`'s current verdict classes are `cell`, `cell_tentative`, `triad`, `variable`, `joint`, `diffuse`, and `none`. A `block` class would sit naturally between `variable` and `joint` / `diffuse`. It would describe the following situation:

> Instrument A's items are internally consistent, and instrument B's items are internally consistent, but the reported cross-correlations imply a canonical correlation above 1.

In terms of the Proposition: $R_{11}$ and $R_{22}$ are each PSD within their rounding boxes, but the block witness certifies $\rho_1 > 1$ in the full rounding box.

**Choosing the partition.** Searching over all partitions grows exponentially with the number of variables, so the partitions would come from the paper's design, for example:

- two questionnaires or scales (items of A versus items of B),
- predictors versus outcomes,
- time point 1 versus time point 2 in a longitudinal design.

`localize_psd_fault()` could accept an optional variable partition, or a list of partitions, and report a block-level fault when the within-block checks pass but the block witness certifies inconsistency.

---

## 5. Block shrinkage for certifying consistency

The interpolation identity from Section 1 also gives a bound useful on the consistency side.

**Claim.** For $R(\tau) = \tau R + (1 - \tau)\operatorname{diag}(R_{11}, R_{22})$,

$$
\lambda_{\min}\big(R(\tau)\big) \;\ge\; \tau\, \lambda_{\min}(R) + (1 - \tau)\, \min\{\lambda_{\min}(R_{11}),\, \lambda_{\min}(R_{22})\}.
$$

*Proof.*

1. $\lambda_{\min}(A) = \min_{\|x\|_2 = 1} x^\top A x$ is a minimum of functions that are linear in $A$, so $\lambda_{\min}$ is concave on symmetric matrices:
   $\lambda_{\min}\big(\tau A + (1-\tau) B\big) \ge \tau\, \lambda_{\min}(A) + (1-\tau)\, \lambda_{\min}(B)$.
2. The eigenvalues of a block-diagonal matrix are those of its blocks, so $\lambda_{\min}\big(\operatorname{diag}(R_{11}, R_{22})\big) = \min\{\lambda_{\min}(R_{11}), \lambda_{\min}(R_{22})\}$. $\square$

**Use.** When the diagonal blocks are well conditioned, shrinking only the **cross-block** cells toward zero, as far as their rounding boxes allow, is a cheap way to construct a candidate PSD matrix inside the box. The candidate can then be verified with Rump's (2006) one-sided Cholesky test, as in `port`'s existing certification path.

`port`'s existing construction, which shrinks the box-edge matrix toward the diagonal, is the special case where every block is a single variable. Block shrinkage keeps the within-block entries fixed, so it searches a different family of candidates.

---

## References

- Hotelling, H. (1936). Relations between two sets of variates. *Biometrika*, 28(3/4), 321–377.
- Hussey, I. (2026). *port: Positive-semi-definiteness Of Rounded correlation Tables* [Computer software]. https://github.com/ianhussey/port. doi:10.5281/zenodo.21434236
- Latała, R., & Matlak, D. (2017). Royen's proof of the Gaussian correlation inequality. In B. Klartag & E. Milman (Eds.), *Geometric Aspects of Functional Analysis* (Lecture Notes in Mathematics, Vol. 2169, pp. 265–275). Springer. Preprint: arXiv:1512.08776
- Royen, T. (2014). A simple proof of the Gaussian correlation conjecture extended to some multivariate gamma distributions. *Far East Journal of Theoretical Statistics*, 48(2), 139–145. arXiv:1408.1028
- Rump, S. M. (2006). Verification of positive definiteness. *BIT Numerical Mathematics*, 46(2), 433–452.
- Šidák, Z. (1967). Rectangular confidence regions for the means of multivariate normal distributions. *Journal of the American Statistical Association*, 62(318), 626–633.
- van Tilburg, W. A. P., & van Tilburg, L. J. A. (2023). Inconsistent hypotheses and effect-size limits. *Advances in Methods and Practices in Psychological Science*, 6(4).
