# -----------------------------------------------------------------------------
# Verified witness-vector bound (primary impossibility path, all p).
#
# For any test direction v != 0, the most PSD-favourable in-box matrix along v
# gives, exactly,
#
#     max_{X in box} v' X v = v'Rv + delta * (||v||_1^2 - ||v||_2^2).
#
# Derivation: X has X_ii = 1 exactly and X_ij in [R_ij - delta, R_ij + delta]
# for i != j. Then
#     v'Xv = sum_i v_i^2 + 2 sum_{i<j} v_i v_j X_ij,
# and each off-diagonal is chosen to maximise 2 v_i v_j X_ij, i.e.
# X_ij = R_ij + delta*sign(v_i v_j), giving the extra term
#     2*delta*sum_{i<j}|v_i||v_j| = delta*(||v||_1^2 - ||v||_2^2).
# Clipping the box to [-1, 1] can only shrink it, so ignoring that clip keeps
# this an UPPER bound on the achievable quadratic form -- sound for
# impossibility. If this maximum is < 0 then no in-box matrix is PSD (a PSD
# matrix would need v'Xv >= 0), and v certifies impossibility.
#
# v need not be an exact eigenvector: it is only a test direction, so its
# inaccuracy never invalidates the certificate. Only the *evaluation* of the
# bound must be rigorous, which is what .witness_bound() guarantees.
# -----------------------------------------------------------------------------

# Evaluate the witness bound for a single direction v, together with a rigorous
# a-priori forward-error bound `E` such that the exact real value M satisfies
# M <= M_hat + E =: B_upper. Hence B_upper < 0 rigorously implies impossibility.
#
# v may have any (nonzero) norm: we use the general form with ||v||_2^2 rather
# than assuming a normalized v, which is strictly more rigorous.
.witness_bound <- function(R, v, delta) {
  p <- length(v)
  u <- .unit_roundoff()

  Rv <- as.numeric(R %*% v)
  a  <- sum(v * Rv)          # a = v'Rv (quadratic form)
  av <- abs(v)
  # P = |v|' |R| |v|: the magnitude scaffold that bounds the error of a.
  P    <- sum(av * as.numeric(abs(R) %*% av))
  L1   <- sum(av)            # ||v||_1
  n2   <- sum(v * v)         # ||v||_2^2
  L1sq <- L1 * L1
  Dterm <- L1sq - n2         # ||v||_1^2 - ||v||_2^2  (>= 0)
  M_hat <- a + delta * Dterm

  # Rigorous forward-error bound. Each coefficient below strictly dominates the
  # true accumulated error of the corresponding sub-computation (Higham Ch. 3):
  #   * v'Rv computed as (R %*% v) then a dot product: error <= 2*gamma_p * P;
  #     we use gamma_{2p+2} * P (dominates) and add u*P to also cover the <= u
  #     representation error of the stored box centres R_ij.
  #   * L1^2, n2: squarings/sums of p terms, dominated by gamma_{2p+2}, gamma_{p+2}.
  #   * the final few combining flops (subtract, scale by delta, add): gamma_8.
  Ea    <- (.gamma(2L * p + 2L) + u) * P
  EL1sq <- .gamma(2L * p + 2L) * L1sq
  En2   <- .gamma(p + 2L) * n2
  Ecomb <- .gamma(8L) * (abs(a) + delta * (L1sq + n2))
  E <- Ea + delta * (EL1sq + En2) + Ecomb
  # Inflate to cover the (tiny) rounding incurred while computing E itself. The
  # 1e-6 relative bump is ~9 orders of magnitude larger than that rounding, so
  # fl(E * 1.000001) rigorously exceeds the true error bound.
  E <- E * 1.000001

  list(v = v, M_hat = M_hat, E = E, B_upper = M_hat + E,
       a = a, L1 = L1, n2 = n2)
}

# Fast O(p^2) precheck: if any reported off-diagonal is out of range even at the
# nearest box edge (|R_ij| - delta > 1), no in-box value can be a valid
# correlation, so the matrix is impossible. Returns NULL, or a one-row summary.
.precheck_range <- function(R, delta) {
  p <- nrow(R)
  off <- which(upper.tri(R), arr.ind = TRUE)
  if (nrow(off) == 0L) return(NULL)
  vals <- R[off]
  viol <- abs(vals) - delta > 1
  if (!any(viol)) return(NULL)
  k <- which(viol)[1L]
  list(i = as.integer(off[k, 1L]), j = as.integer(off[k, 2L]),
       value = as.numeric(vals[k]))
}

# Search a family of sound test directions and return the one giving the most
# negative (most impossibility-favourable) B_upper. Every direction yields a
# valid certificate, so combining several only increases detection power.
#
# Directions:
#   1. the bottom few eigenvectors of R (minimise v'Rv);
#   2. all coordinate pairs v = (e_i - sign(R_ij) e_j)/sqrt(2), which reproduce
#      the pairwise range condition and subsume the precheck;
#   3. for small p, the smallest eigenvector of every 3x3 principal submatrix
#      (connects to the exact 3x3 determinant condition).
.witness_search <- function(R, delta, triples_max_p = 12L) {
  p <- nrow(R)
  cands <- vector("list", 0L)

  eg <- eigen(R, symmetric = TRUE)
  k <- min(p, 4L)
  for (j in seq_len(k)) {
    cands[[length(cands) + 1L]] <- eg$vectors[, p - j + 1L]
  }

  for (i in seq_len(p - 1L)) {
    for (j in (i + 1L):p) {
      s <- if (R[i, j] >= 0) 1 else -1
      v <- numeric(p)
      v[i] <- 1 / sqrt(2)
      v[j] <- -s / sqrt(2)
      cands[[length(cands) + 1L]] <- v
    }
  }

  if (p >= 3L && p <= triples_max_p) {
    combs <- utils::combn(p, 3L)
    for (c in seq_len(ncol(combs))) {
      S <- combs[, c]
      es <- eigen(R[S, S, drop = FALSE], symmetric = TRUE)
      v <- numeric(p)
      v[S] <- es$vectors[, 3L]
      cands[[length(cands) + 1L]] <- v
    }
  }

  best <- NULL
  for (v in cands) {
    if (all(v == 0)) next
    wb <- .witness_bound(R, v, delta)
    if (is.null(best) || wb$B_upper < best$B_upper) best <- wb
  }
  best
}
