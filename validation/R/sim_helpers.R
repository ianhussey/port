# =============================================================================
# Shared helpers for the port validation simulations.
#
# Contents:
#   1. Matrix generators (each returns a p x p correlation matrix, unit diagonal)
#   2. Contamination (round off-diagonals only; perturbations for Sim 3)
#   3. The INDEPENDENT ORACLE for the box-feasibility question (ground truth)
#   4. Small analyse() wrappers around the package under test
#   5. Auditing / reproducibility utilities
#
# The oracle is deliberately NOT the tool under test and NOT a point lambda_min
# tolerance. Ground truth answers the BOX question: "does a PSD matrix exist
# inside the rounding box around the reported entries?" -- see box_is_feasible().
# =============================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# -----------------------------------------------------------------------------
# 1. Generators
# -----------------------------------------------------------------------------

# Onion method (clusterGeneration): principled uniform-ish coverage of valid
# correlation matrices. eta = 1 spreads over the valid set; small eta piles mass
# on the PSD boundary (near-singular) -- the boundary stress source.
gen_onion <- function(p, eta = 1) {
  S <- clusterGeneration::genPositiveDefMat(
    dim = p,
    covMethod = "onion",
    eta = eta
  )$Sigma
  stats::cov2cor(S)
}

# Factor model: the primary plausibility generator. A positive general factor
# plus k-1 group factors gives a positive-manifold, heterogeneous, moderate-
# magnitude, mildly-redundant structure resembling real psychology matrices.
# Stronger group loadings (loading_sd) raise communalities and move toward the
# PSD boundary. Diagonal is exactly 1 by construction.
gen_factor_model <- function(
  p,
  k = 3,
  loading_sd = 0.35,
  general_min = 0.30,
  general_max = 0.70,
  comm_cap = 0.98
) {
  L <- matrix(stats::rnorm(p * k, 0, loading_sd), p, k)
  L[, 1] <- stats::runif(p, general_min, general_max) # positive general factor
  comm <- rowSums(L^2)
  # rescale any row whose communality would leave no uniqueness
  s <- ifelse(comm >= comm_cap, sqrt(comm_cap / comm), 1)
  L <- L * s
  comm <- rowSums(L^2)
  R <- tcrossprod(L) + diag(1 - comm, p) # PD, unit diagonal
  stats::cov2cor(R)
}

# Equicorrelation: closed-form edge-case oracle (eigenvalues 1+(p-1)rho and
# 1-rho, the latter with multiplicity p-1). Use in tests, not as a main DGP.
gen_equicorrelation <- function(p, rho) {
  R <- matrix(rho, p, p)
  diag(R) <- 1
  R
}

# Independent sampling: trivial, well-conditioned. Sanity row only.
gen_independent <- function(p, n = 250) {
  X <- matrix(stats::rnorm(n * p), n, p)
  stats::cor(X)
}

# Dispatch a generator by name with a per-row seed (so the same base matrix is
# reused across downstream rounding / perturbation / tier conditions).
generate_data <- function(
  generator,
  p,
  seed,
  eta = 1,
  k = 3,
  loading_sd = 0.35,
  rho = 0.3,
  n_independent = 250
) {
  set.seed(seed)
  switch(
    generator,
    onion = gen_onion(p, eta = eta),
    factor_model = gen_factor_model(p, k = k, loading_sd = loading_sd),
    equicorrelation = gen_equicorrelation(p, rho = rho),
    independent = gen_independent(p, n = n_independent),
    stop("unknown generator: ", generator)
  )
}

# -----------------------------------------------------------------------------
# 2. Contamination (Monticelli-style: round ONLY off-diagonals, mirror to keep
#    symmetry, diagonal fixed at 1). Perturbations for Sim 3 come first.
# -----------------------------------------------------------------------------

.rounder <- function(method, decimals) {
  m <- 10^decimals
  switch(
    method,
    nearest = function(x) round(x, decimals), # round-half-to-even: error <= 0.5*10^-d
    truncate = function(x) trunc(x * m) / m, # toward zero: error < 10^-d
    floor = function(x) floor(x * m) / m,
    ceiling = function(x) ceiling(x * m) / m,
    stop("unknown rounding method: ", method)
  )
}

round_offdiag <- function(R, decimals, method = "nearest") {
  p <- nrow(R)
  f <- .rounder(method, decimals)
  Rt <- R
  ut <- upper.tri(R)
  Rt[ut] <- f(R[ut])
  Rt[lower.tri(Rt)] <- t(Rt)[lower.tri(Rt)] # mirror upper -> lower
  diag(Rt) <- 1
  Rt
}

# Perturbations matched to the localization taxonomy. Returns the perturbed
# matrix and the "true culprit" (cell or variable) where one is well defined.
perturb <- function(
  R,
  type,
  magnitude,
  i = NULL,
  j = NULL,
  var = NULL,
  seed = NULL
) {
  p <- nrow(R)
  Rp <- R
  clamp <- function(x) pmax(pmin(x, 1), -1)
  culprit_cell <- NULL
  culprit_var <- NULL
  if (type == "single_cell") {
    push <- if (R[i, j] >= 0) 1 else -1 # push away from 0 -> inconsistency
    Rp[i, j] <- Rp[j, i] <- clamp(R[i, j] + push * magnitude)
    culprit_cell <- c(i, j)
  } else if (type == "sign_flip") {
    Rp[i, j] <- Rp[j, i] <- -R[i, j]
    culprit_cell <- c(i, j)
  } else if (type == "variable_column") {
    others <- setdiff(seq_len(p), var)
    push <- ifelse(R[var, others] >= 0, 1, -1)
    Rp[var, others] <- clamp(R[var, others] + push * magnitude)
    Rp[others, var] <- Rp[var, others]
    culprit_var <- var
  } else if (type == "jitter") {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    ut <- upper.tri(R)
    noise <- stats::runif(sum(ut), -magnitude, magnitude)
    Rp[ut] <- clamp(R[ut] + noise)
    Rp[lower.tri(Rp)] <- t(Rp)[lower.tri(Rp)]
  } else {
    stop("unknown perturbation type: ", type)
  }
  diag(Rp) <- 1
  list(R = Rp, culprit_cell = culprit_cell, culprit_var = culprit_var)
}

# One-call contamination used by the sims: optionally perturb, then round.
contaminate_data <- function(
  R,
  decimals,
  rounding_method = "nearest",
  perturbation = NULL
) {
  culprit_cell <- NULL
  culprit_var <- NULL
  if (!is.null(perturbation)) {
    pr <- do.call(perturb, c(list(R = R), perturbation))
    R <- pr$R
    culprit_cell <- pr$culprit_cell
    culprit_var <- pr$culprit_var
  }
  list(
    R_tilde = round_offdiag(R, decimals, rounding_method),
    culprit_cell = culprit_cell,
    culprit_var = culprit_var
  )
}

# -----------------------------------------------------------------------------
# 3. The INDEPENDENT ORACLE
# -----------------------------------------------------------------------------

# Frobenius distance from a symmetric matrix to the PSD cone: sqrt of the sum of
# squares of its negative eigenvalues. Pure eigendecomposition -- independent of
# the tool under test.
dist_to_psd <- function(M) {
  ev <- eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values
  sqrt(sum(pmin(ev, 0)^2))
}

# Max Frobenius distance from the box centre to any in-box matrix. Each of the
# p(p-1) off-diagonal entries moves at most delta (clipping to [-1,1] only
# shrinks the box; the diagonal is fixed), so this bounds ||B - centre||_F.
box_radius_frob <- function(p, delta) delta * sqrt(p * (p - 1))

# Is a specific matrix W inside the rounding box around `center`?
in_box <- function(W, center, delta, tol = 1e-10) {
  if (any(abs(diag(W) - 1) > 1e-8)) {
    return(FALSE)
  }
  off <- upper.tri(center)
  lo <- pmax(center[off] - delta, -1) - tol
  hi <- pmin(center[off] + delta, 1) + tol
  w <- W[off]
  all(w >= lo & w <= hi)
}

is_psd <- function(M, tol = 1e-10) {
  min(eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values) >=
    -tol
}

# Independent nearest-correlation-matrix (Higham 2002, Dykstra-corrected
# alternating projection) used ONLY to generate a feasibility witness. Its output
# is independently verified PSD (eigen) and in-box (arithmetic); it is written
# here from scratch rather than calling the package's POCS, so the oracle does
# not depend on the tool under test.
oracle_nearcorr <- function(R, max_iter = 300L, tol = 1e-10) {
  n <- nrow(R)
  Y <- R
  Dz <- matrix(0, n, n)
  for (it in seq_len(max_iter)) {
    Rk <- Y - Dz
    ev <- eigen((Rk + t(Rk)) / 2, symmetric = TRUE)
    X <- ev$vectors %*% (pmax(ev$values, 0) * t(ev$vectors))
    Dz <- X - Rk
    Yn <- X
    diag(Yn) <- 1
    Yn <- pmin(pmax(Yn, -1), 1)
    diag(Yn) <- 1
    if (max(abs(Yn - Y)) < tol) {
      Y <- Yn
      break
    }
    Y <- Yn
  }
  (Y + t(Y)) / 2
}

# Ground-truth answer to the BOX question. Returns a list with `status` in
# {"feasible", "infeasible", "uncertain"} plus the certificate and diagnostics.
#
#  * infeasible: rigorous -- dist(center, cone) exceeds the box radius by a margin
#    that dominates the eigenvalue rounding error, so NO in-box matrix is PSD.
#  * feasible: an exhibited matrix (a supplied known-PSD `witness`, or the
#    independent nearest-correlation matrix) is verified PSD and inside the box.
#  * uncertain: neither certificate fired (a thin near-boundary band); excluded
#    from rate denominators and reported.
box_is_feasible <- function(
  R_reported,
  delta,
  witness = NULL,
  margin = 1e-6,
  try_nearcorr = TRUE
) {
  p <- nrow(R_reported)
  radius <- box_radius_frob(p, delta)
  dist <- dist_to_psd(R_reported)
  if (dist > radius + margin) {
    return(list(
      status = "infeasible",
      cert = "distance>radius",
      dist = dist,
      radius = radius
    ))
  }
  if (
    !is.null(witness) && in_box(witness, R_reported, delta) && is_psd(witness)
  ) {
    return(list(
      status = "feasible",
      cert = "known_psd_in_box",
      dist = dist,
      radius = radius
    ))
  }
  if (try_nearcorr) {
    W <- oracle_nearcorr(R_reported)
    if (in_box(W, R_reported, delta) && is_psd(W)) {
      return(list(
        status = "feasible",
        cert = "nearcorr_in_box",
        dist = dist,
        radius = radius
      ))
    }
  }
  list(status = "uncertain", cert = NA_character_, dist = dist, radius = radius)
}

# -----------------------------------------------------------------------------
# 3b. Disattenuation DGP + heterogeneous-box oracle (Sim 4)
# -----------------------------------------------------------------------------

# Attenuate a construct matrix by per-variable reliabilities (the inverse of the
# package's disattenuate): R_ij = C_ij * sqrt(rho_i * rho_j), diagonal 1. The
# result is ALWAYS a valid correlation matrix -- attenuation shrinks correlations
# toward the PSD interior -- so the "observed" matrix is valid by construction.
attenuate <- function(C, reliability) {
  p <- nrow(C)
  rel <- if (length(reliability) == 1L) rep(reliability, p) else reliability
  s <- sqrt(rel)
  R <- outer(s, s) * C
  diag(R) <- 1
  R
}

# Draw per-variable reliabilities from a regime (low includes the poor-measure
# range the tool is meant to catch).
gen_reliability <- function(p, regime = c("low", "moderate", "high")) {
  regime <- match.arg(regime)
  rng <- switch(
    regime,
    low = c(0.30, 0.60),
    moderate = c(0.60, 0.85),
    high = c(0.85, 0.98)
  )
  stats::runif(p, rng[1], rng[2])
}

# Independent construction of the disattenuated box (per-cell interval division
# of the correlation box by the reliability box). Written from scratch so the
# oracle does not depend on the package internals.
disatten_box <- function(R, reliability, delta_R, delta_rel, eps = 1e-6) {
  p <- nrow(R)
  rel <- if (length(reliability) == 1L) rep(reliability, p) else reliability
  rel_lo <- pmax(rel - delta_rel, eps)
  rel_hi <- pmin(rel + delta_rel, 1)
  s_lo <- 1 / sqrt(rel_hi)
  s_hi <- 1 / sqrt(rel_lo)
  lo <- matrix(0, p, p)
  hi <- matrix(0, p, p)
  for (i in seq_len(p - 1L)) {
    for (j in (i + 1L):p) {
      rlo <- R[i, j] - delta_R
      rhi <- R[i, j] + delta_R
      pr <- c(
        s_lo[i] * s_lo[j] * rlo,
        s_lo[i] * s_lo[j] * rhi,
        s_hi[i] * s_hi[j] * rlo,
        s_hi[i] * s_hi[j] * rhi
      )
      lo[i, j] <- lo[j, i] <- max(min(pr), -1)
      hi[i, j] <- hi[j, i] <- min(max(pr), 1)
    }
  }
  diag(lo) <- 1
  diag(hi) <- 1
  list(lo = lo, hi = hi)
}

# Heterogeneous-box feasibility oracle: does a PSD matrix exist in the box
# [lo, hi]? Same certificates as box_is_feasible(), generalised to per-cell
# radii: infeasible if the box centre's distance to the PSD cone exceeds the box
# radius sqrt(sum of squared per-cell radii); feasible via an exhibited PSD
# witness (a known construct matrix, or the nearest-correlation matrix).
box_is_feasible_hetero <- function(
  lo,
  hi,
  witness = NULL,
  margin = 1e-6,
  try_nearcorr = TRUE
) {
  off <- upper.tri(lo) | lower.tri(lo)
  center <- (lo + hi) / 2
  diag(center) <- 1
  r <- (hi - lo) / 2
  box_radius <- sqrt(sum((r[off])^2))
  dist <- dist_to_psd(center)
  if (dist > box_radius + margin) {
    return(list(status = "infeasible", dist = dist, radius = box_radius))
  }
  in_box_h <- function(X) {
    if (any(abs(diag(X) - 1) > 1e-8)) {
      return(FALSE)
    }
    all(X[off] >= lo[off] - 1e-9 & X[off] <= hi[off] + 1e-9)
  }
  if (!is.null(witness) && in_box_h(witness) && is_psd(witness)) {
    return(list(
      status = "feasible",
      cert = "known_psd",
      dist = dist,
      radius = box_radius
    ))
  }
  if (try_nearcorr) {
    W <- oracle_nearcorr(center)
    if (in_box_h(W) && is_psd(W)) {
      return(list(
        status = "feasible",
        cert = "nearcorr",
        dist = dist,
        radius = box_radius
      ))
    }
  }
  list(status = "uncertain", dist = dist, radius = box_radius)
}

# Independent critical common-reliability oracle: bisect rho for the point
# disattenuated matrix's feasibility (min eigenvalue >= 0), combined with the
# range bound max|R_ij|. Validates the package's closed form rho* = 1 - lambda_min
# WITHOUT using that formula.
critical_rho_oracle <- function(R, tol = 1e-6, max_it = 60L) {
  rho_range <- max(abs(R[upper.tri(R)]))
  feas_psd <- function(rho) {
    D <- R / rho
    diag(D) <- 1
    min(eigen(D, symmetric = TRUE, only.values = TRUE)$values) >= 0
  }
  if (!feas_psd(1)) {
    return(list(rho_star = NA_real_, rho_psd = NA_real_, rho_range = rho_range))
  }
  lo <- 1e-4
  hi <- 1
  for (it in seq_len(max_it)) {
    if (hi - lo < tol) {
      break
    }
    m <- (lo + hi) / 2
    if (feas_psd(m)) hi <- m else lo <- m
  }
  list(rho_star = max(hi, rho_range), rho_psd = hi, rho_range = rho_range)
}

# -----------------------------------------------------------------------------
# 4. analyse() wrappers around the package under test
# -----------------------------------------------------------------------------

analyse_check <- function(R, decimals = NULL, delta = NULL) {
  res <- tryCatch(
    if (!is.null(delta)) {
      port::check_corr_psd(R, delta = delta)
    } else {
      port::check_corr_psd(R, decimals = decimals)
    },
    error = function(e) NULL
  )
  if (is.null(res)) {
    return(tibble::tibble(
      verdict = "error",
      tier = NA_character_,
      margin = NA_real_
    ))
  }
  tibble::tibble(
    verdict = res$verdict,
    tier = res$tier %||% NA_character_,
    margin = res$margin %||% NA_real_
  )
}

# Localization wrapper: returns the localization verdict, a flat string of the
# implicated cells/variable, whether the known culprit was recovered, and the
# severity. `culprit_cell` / `culprit_var` encode the planted fault (or NULL).
analyse_localize <- function(
  R,
  decimals,
  culprit_cell = NULL,
  culprit_var = NULL
) {
  lf <- tryCatch(
    port::localize_psd_fault(R, decimals = decimals),
    error = function(e) NULL
  )
  if (is.null(lf)) {
    return(tibble::tibble(
      localization_verdict = "error",
      implicated = NA_character_,
      hit_cell = NA,
      hit_variable = NA,
      severity_max = NA_real_
    ))
  }
  cells <- lf$implicated$cells
  var <- lf$implicated$variable
  cell_str <- if (!is.null(cells) && nrow(cells) > 0L) {
    paste(
      apply(cells, 1L, function(r) sprintf("(%d,%d)", r[1], r[2])),
      collapse = " "
    )
  } else {
    NA_character_
  }
  imp_str <- paste(
    c(if (!is.null(var)) sprintf("var%d", var), cell_str %||% NULL),
    collapse = "; "
  )

  # Did localization recover the planted culprit?
  hit_cell <- NA
  hit_variable <- NA
  cell_set <- if (!is.null(cells) && nrow(cells) > 0L) {
    apply(cells, 1L, function(r) paste(sort(r), collapse = "-"))
  } else {
    character(0)
  }
  if (!is.null(culprit_cell)) {
    key <- paste(sort(culprit_cell), collapse = "-")
    endpoints <- culprit_cell
    hit_cell <- (key %in% cell_set) ||
      (!is.null(var) && var %in% endpoints) # a cell fault may read as a variable
  }
  if (!is.null(culprit_var)) {
    in_col <- length(cell_set) > 0L &&
      all(vapply(
        strsplit(cell_set, "-"),
        function(z) culprit_var %in% as.integer(z),
        logical(1)
      ))
    hit_variable <- (!is.null(var) && var == culprit_var) || in_col
  }
  tibble::tibble(
    localization_verdict = lf$localization_verdict,
    implicated = imp_str,
    hit_cell = hit_cell,
    hit_variable = hit_variable,
    severity_max = lf$severity$severity_max %||% NA_real_
  )
}

# Wrapper for check_disattenuated_psd(): returns the disattenuation verdict.
analyse_disattenuated <- function(
  R,
  reliability,
  decimals,
  reliability_decimals,
  max_plausible_r = 1
) {
  res <- tryCatch(
    port::check_disattenuated_psd(
      R,
      reliability = reliability,
      decimals = decimals,
      reliability_decimals = reliability_decimals,
      max_plausible_r = max_plausible_r
    ),
    error = function(e) NULL
  )
  if (is.null(res)) {
    return(tibble::tibble(verdict = "error"))
  }
  tibble::tibble(verdict = res$verdict)
}

# -----------------------------------------------------------------------------
# 5. Auditing / reproducibility utilities
# -----------------------------------------------------------------------------

# One-row summary of the off-diagonal correlation distribution, so generator
# plausibility is auditable in Results (PSD alone is too thin a description).
offdiag_summary <- function(R) {
  x <- R[upper.tri(R)]
  q <- stats::quantile(x, c(.05, .25, .5, .75, .95), names = FALSE)
  tibble::tibble(
    offdiag_mean = mean(x),
    offdiag_sd = stats::sd(x),
    offdiag_min = min(x),
    offdiag_max = max(x),
    offdiag_q05 = q[1],
    offdiag_q25 = q[2],
    offdiag_q50 = q[3],
    offdiag_q75 = q[4],
    offdiag_q95 = q[5]
  )
}

# Deterministic per-replication seed from the base-matrix identity only
# (generator, p, iteration) so downstream conditions reuse the same base matrix.
GENERATOR_INDEX <- c(
  onion = 1L,
  factor_model = 2L,
  equicorrelation = 3L,
  independent = 4L
)
seed_from <- function(generator, p, iteration, salt = 0L) {
  g <- GENERATOR_INDEX[[generator]] %||% 9L
  as.integer((g * 1e7 + p * 1e4 + iteration + salt) %% .Machine$integer.max)
}

# Wilson score interval for a proportion (robust when phat is at 0 or 1, where
# the Wald/binomial MCSE collapses and is misleading).
wilson_ci <- function(x, n, conf = 0.95) {
  if (n == 0L) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  z <- stats::qnorm(1 - (1 - conf) / 2)
  phat <- x / n
  denom <- 1 + z^2 / n
  centre <- (phat + z^2 / (2 * n)) / denom
  half <- z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2)) / denom
  c(lower = max(0, centre - half), upper = min(1, centre + half))
}
