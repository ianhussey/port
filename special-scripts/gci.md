# The Gaussian correlation inequality gives `port` no new test

## Summary

The Gaussian correlation inequality (GCI) and positive semidefiniteness (PSD) of correlation matrices are connected, but only in one direction. PSD is a *precondition* of the GCI and a tool in Royen's (2014) proof. The GCI itself produces no new constraint on correlation matrices, so it offers no new test of the kind `port` performs.

Royen's proof does put the block structure of a covariance matrix on display (Schur complements and canonical correlations). Three ideas for `port` were taken from it. None of them paid off, and none depends on the GCI; they are recorded in `sparse-witness.md`, Section 5, next to the improvement that the same investigation did lead to.

---

## 1. The inequality

Let $X \sim \mathcal{N}(0, C)$ be an $n$-dimensional centred Gaussian vector, with coordinates split into two blocks of sizes $n_1 + n_2 = n$:

$$
C = \begin{pmatrix} C_{11} & C_{12} \\ C_{21} & C_{22} \end{pmatrix}.
$$

In the form Royen proves, the GCI states that for any thresholds $t_i > 0$,

$$
P\big(|X_i| \le t_i \;\; \forall i \le n\big) \;\ge\; P\big(|X_i| \le t_i \;\; \forall i \le n_1\big)\; P\big(|X_i| \le t_i \;\; \forall\, n_1 < i \le n\big).
$$

The general statement for symmetric convex sets $K, L$, $\mu(K \cap L) \ge \mu(K)\,\mu(L)$, follows by approximation: every closed symmetric convex set is a countable intersection of symmetric strips $\{x : |\langle x, v\rangle| \le t\}$ (Latała & Matlak, 2017; Panchenko, 2020).

---

## 2. What PSD does in Royen's proof

Latała and Matlak (2017) give the proof in full. Only two of its steps involve PSD.

**The interpolation path stays in the PSD cone.** Define, for $\tau \in [0, 1]$,

$$
C(\tau) = \begin{pmatrix} C_{11} & \tau\, C_{12} \\ \tau\, C_{21} & C_{22} \end{pmatrix}
= \tau\, C + (1 - \tau) \operatorname{diag}(C_{11}, C_{22}).
$$

At $\tau = 0$ the blocks are independent, so the joint probability equals the product. At $\tau = 1$ it is the joint probability under $C$. The proof shows that the probability is nondecreasing in $\tau$. For this, every $C(\tau)$ must be a covariance matrix. It is one because the PSD cone is convex, and $\operatorname{diag}(C_{11}, C_{22})$ is PSD since every principal submatrix of a PSD matrix is PSD.

**The principal minors decrease along the path.** The Laplace transform of the squared coordinates is $\mathbb{E}\exp\!\big(-\sum_i \lambda_i X_i^2\big) = |I + 2\Lambda C|^{-1/2}$ with $\Lambda = \operatorname{diag}(\lambda_i)$, and the determinant expands in the principal minors $|C(\tau)_J|$. A minor whose index set $J$ lies within one block does not depend on $\tau$. If $J = J_1 \cup J_2$ meets both blocks, the Schur complement gives (taking $C$ positive definite, without loss of generality)

$$
\big|C(\tau)_J\big| = |C_{J_1}|\,|C_{J_2}| \prod_k \big(1 - \tau^2 \rho_k^2\big),
$$

where the $\rho_k$ are the **canonical correlations** between the variable groups $J_1$ and $J_2$ (Hotelling, 1936). Because $C$ is PSD, $\rho_k \le 1$, so every factor lies in $[0, 1]$ and is nonincreasing in $\tau$. Royen's interpolation thus turns down every canonical correlation between the two blocks by the factor $\tau$.

The rest of the proof turns these sign facts into a monotone probability. That step rests on the multivariate gamma distribution of the squares, not on PSD.

PSD therefore has three jobs in the proof:

1. It guarantees that the Gaussian $\mathcal{N}(0, C)$ exists.
2. It keeps the interpolation path $C(\tau)$ inside the PSD cone.
3. It makes the principal minors nonnegative and monotone along that path.

---

## 3. Why the GCI gives no forensic test

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

## References

- Hotelling, H. (1936). Relations between two sets of variates. *Biometrika*, 28(3/4), 321–377.
- Khatri, C. G. (1967). On certain inequalities for normal distributions and their applications to simultaneous confidence bounds. *Annals of Mathematical Statistics*, 38(6), 1853–1867.
- Latała, R., & Matlak, D. (2017). Royen's proof of the Gaussian correlation inequality. In B. Klartag & E. Milman (Eds.), *Geometric Aspects of Functional Analysis* (Lecture Notes in Mathematics, Vol. 2169, pp. 265–275). Springer. Preprint: arXiv:1512.08776
- Panchenko, D. (2020). *Gaussian correlation inequality* [Lecture notes]. https://sites.stat.washington.edu/people/mdperlma/542%20Multivariate%20Analysis%20Material/Panchenko%20re%20Royen.pdf
- Royen, T. (2014). A simple proof of the Gaussian correlation conjecture extended to some multivariate gamma distributions. *Far East Journal of Theoretical Statistics*, 48(2), 139–145. arXiv:1408.1028
- Šidák, Z. (1967). Rectangular confidence regions for the means of multivariate normal distributions. *Journal of the American Statistical Association*, 62(318), 626–633.
