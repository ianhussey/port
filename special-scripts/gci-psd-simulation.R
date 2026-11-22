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
# consistency. The original search is recovered by stubbing out
# greedy_subset_witness(); it is restored afterwards.
#
# Run from the package root:  Rscript special-scripts/gci-psd-simulation.R
# Takes about 2-4 minutes.

suppressMessages(pkgload::load_all(".", quiet = TRUE))

greedy_real <- greedy_subset_witness
use_greedy <- function(on) {
  f <- if (on) greedy_real else function(...) NULL
  utils::assignInNamespace("greedy_subset_witness", f, "port")
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

    use_greedy(FALSE)
    old_inc <- .witness_search(R, d)$B_upper < 0
    use_greedy(TRUE)
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
      function(S) neg(subset_bottom_eigvec(R, S)),
      logical(1)
    ))
    sub_all <- any(vapply(
      subsets(p),
      function(S) neg(subset_bottom_eigvec(R, S)),
      logical(1)
    ))
    greedy <- isTRUE(
      greedy_subset_witness(R, function(v) .witness_bound(R, v, d))$B_upper < 0
    )
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
        greedy,
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
