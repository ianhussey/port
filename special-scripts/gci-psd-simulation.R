# Simulation behind Section 6 of gci-psd-port.md.
#
# Question: on rounded correlation matrices that port's original witness search
# leaves undecided, which extra families of test directions (or constructions)
# turn the verdict into a certificate?
#
# Generator: a low-rank-plus-small-noise correlation matrix (near the PSD
# boundary), one cell perturbed by up to +-0.08, rounded to 2 decimals
# (delta = 0.005). "Undecided" = the ORIGINAL witness search finds no
# certificate (B_upper >= 0) and neither the constructions nor POCS certify
# consistency. The original search is .witness_search(sparse = FALSE).
#
# Run from the package root:  Rscript special-scripts/gci-psd-simulation.R
# Takes about 5 minutes.

suppressMessages(pkgload::load_all(".", quiet = TRUE))

# Best single test direction, to numerical precision. Within one cell (support
# T, sign pattern s) the box bound is the quadratic form of X_s = R + d * s s'
# (unit diagonal) on T, so a minimiser is an eigenvector of X_s[T, T] whose
# signs match s. Enumerating all (3^p - 1) / 2 cells therefore finds the global
# minimum; the result still goes through the rigorous bound like any direction.
global_dir <- function(R, d) {
  p <- nrow(R)
  best <- Inf
  best_v <- NULL
  for (m in seq_len(2^p - 1)) {
    S <- which(bitwAnd(m, 2^(0:(p - 1))) > 0)
    k <- length(S)
    if (k < 2) {
      next
    }
    for (sg in 0:(2^(k - 1) - 1)) {
      s <- c(1, ifelse(bitwAnd(sg, 2^(seq_len(k - 1) - 1)) > 0, -1, 1))
      X <- R[S, S] + d * outer(s, s)
      diag(X) <- 1
      e <- eigen(X, symmetric = TRUE)
      for (j in k:1) {
        if (e$values[j] >= best) {
          break
        }
        v <- e$vectors[, j]
        if (v[1] < 0) {
          v <- -v
        }
        if (all(v * s > -1e-12)) {
          best <- e$values[j]
          best_v <- replace(numeric(p), S, v)
          break
        }
      }
    }
  }
  best_v
}

msqrt_inv <- function(A) {
  e <- eigen(A, symmetric = TRUE)
  e$vectors %*% (1 / sqrt(e$values) * t(e$vectors))
}

# Block (canonical-correlation) witness for the partition S | rest; NULL unless
# both diagonal blocks are positive definite.
block_dir <- function(R, S) {
  Tt <- setdiff(seq_len(nrow(R)), S)
  A <- R[S, S, drop = FALSE]
  B <- R[Tt, Tt, drop = FALSE]
  pd <- function(X) min(eigen(X, TRUE, TRUE)$values) > 1e-10
  if (!pd(A) || !pd(B)) {
    return(NULL)
  }
  Ai <- msqrt_inv(A)
  Bi <- msqrt_inv(B)
  sv <- svd(Ai %*% R[S, Tt, drop = FALSE] %*% Bi)
  v <- numeric(nrow(R))
  v[S] <- Ai %*% sv$u[, 1]
  v[Tt] <- -Bi %*% sv$v[, 1]
  v
}

bipartitions <- function(p) {
  lapply(seq_len(2^(p - 1) - 1), function(m) {
    which(bitwAnd(m, 2^(0:(p - 1))) > 0)
  })
}
subsets <- function(p, max_size = p) {
  unlist(
    lapply(3:max_size, function(m) utils::combn(p, m, simplify = FALSE)),
    recursive = FALSE
  )
}

cross_shrink <- function(R, S, d) {
  X <- R
  Tt <- setdiff(seq_len(nrow(R)), S)
  X[S, Tt] <- sign(R[S, Tt]) * pmax(0, abs(R[S, Tt]) - d)
  X[Tt, S] <- t(X[S, Tt])
  X
}

simulate <- function(p, n, d = 0.005) {
  verdicts <- character(0)
  undec <- NULL
  triples_lost <- 0L
  example <- NULL
  for (k in seq_len(n)) {
    r <- sample(2:max(2, p - 2), 1)
    L <- matrix(rnorm(p * r), p)
    S0 <- cov2cor(tcrossprod(L) + diag(runif(p, 0.001, 0.05)))
    ij <- sample(p, 2)
    S0[ij[1], ij[2]] <- S0[ij[2], ij[1]] <- S0[ij[1], ij[2]] +
      runif(1, -0.08, 0.08)
    if (any(abs(S0[upper.tri(S0)]) > 1)) {
      next
    }
    R <- round(S0, 2)

    old_inc <- .witness_search(R, d, sparse = FALSE)$B_upper < 0
    if (old_inc) {
      verdicts <- c(verdicts, "inconsistent")
      # Is the triples family redundant once greedy elimination is added?
      if (.witness_search(R, d, triples_max_p = 0L)$B_upper >= 0) {
        triples_lost <- triples_lost + 1L
      }
      next
    }
    if (!is.null(.construct_consistent(R, d))) {
      verdicts <- c(verdicts, "consistent (construction)")
      next
    }
    if (!is.null(.pocs_consistent(R, d))) {
      verdicts <- c(verdicts, "consistent (POCS)")
      next
    }
    verdicts <- c(verdicts, "undecided")

    neg <- function(v) !is.null(v) && .witness_bound(R, v, d)$B_upper < 0
    s_opt <- any(vapply(
      seq_len(p),
      function(j) {
        dd <- .rsquared_direction(R, j)
        if (is.null(dd)) {
          return(FALSE)
        }
        v0 <- dd$v / -dd$v[j]
        f <- function(s) {
          v <- v0
          v[j] <- -s
          .witness_bound(R, v, d)$B_upper
        }
        stats::optimize(f, c(0.01, 3))$objective < 0
      },
      logical(1)
    ))
    block <- if (p <= 8) {
      any(vapply(bipartitions(p), function(S) neg(block_dir(R, S)), logical(1)))
    } else {
      NA
    }
    shrink <- any(vapply(
      bipartitions(p),
      function(S) {
        .verify_psd(cross_shrink(R, S, d))
      },
      logical(1)
    ))
    sub_le4 <- any(vapply(
      subsets(p, 4),
      function(S) neg(.subset_bottom_eigvec(R, S)),
      logical(1)
    ))
    sub_all <- any(vapply(
      subsets(p),
      function(S) neg(.subset_bottom_eigvec(R, S)),
      logical(1)
    ))
    # greedy elimination without polishing, then as implemented (each round's
    # winner polished on its support), then the single-vector optimum
    raw <- function(v) .witness_bound(R, v, d)
    pol <- function(v) .polish_witness(R, d, v, on_support = TRUE)
    greedy_raw <- isTRUE(.greedy_subset_witness(R, raw, raw)$B_upper < 0)
    greedy <- isTRUE(.greedy_subset_witness(R, raw, pol)$B_upper < 0)
    global <- neg(global_dir(R, d))
    new_inc <- check_corr_psd(R, delta = d)$verdict == "inconsistent"
    if (new_inc && is.null(example)) {
      example <- R
    }

    undec <- rbind(
      undec,
      data.frame(
        s_opt,
        block,
        shrink,
        sub_le4,
        sub_all,
        greedy_raw,
        greedy,
        global,
        new_inc
      )
    )
  }
  list(
    verdicts = table(verdicts),
    undec = undec,
    triples_lost = triples_lost,
    example = example
  )
}

set.seed(2026)
for (p in c(5, 6, 8, 10)) {
  res <- simulate(p, n = if (p == 10) 1500 else 2500)
  cat("\n=== p =", p, "===\n")
  print(res$verdicts)
  cat("Undecided cases certified inconsistent by:\n")
  print(colSums(res$undec, na.rm = TRUE))
  cat("(block = NA when p > 8; shrink counts consistency certificates)\n")
  cat(
    "Original certificates lost if triples are dropped (greedy kept):",
    res$triples_lost,
    "\n"
  )
  if (p == 6 && !is.null(res$example)) {
    cat("Example matrix undecided before, inconsistent now:\n")
    dput(res$example)
  }
}
