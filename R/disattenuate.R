# -----------------------------------------------------------------------------
# Reliability disattenuation, box-aware and rounding-aware.
#
# Spearman's correction for attenuation inflates each observed correlation by the
# geometric mean of the two variables' reliabilities:
#     D_ij = R_ij / sqrt(rho_i * rho_j),   D_ii = 1.
# Because reliabilities are in (0, 1], this only ever moves the matrix AWAY from
# the PSD cone, so a valid observed matrix can become an inconsistent *construct*
# matrix under the reliabilities the authors themselves report. This module asks
# whether that disattenuated matrix is inconsistent / consistent but implausible / consistent,
# accounting for the rounding of both the correlations and the reliabilities.
#
# Key facts used (single common reliability rho, exact values):
#   * D = I + (1/rho) * (R - I), so the disattenuated matrix loses positive
#     semidefiniteness exactly below rho*_psd = 1 - lambda_min(R);
#   * some |D_ij| exceeds 1 below rho*_range = max|R_ij|;
#   * some |D_ij| exceeds a plausibility cutoff c below max|R_ij| / c.
# The forward verdict at reported reliabilities is the rigorous, box-aware one
# (it reuses the package's witness / Rump box machinery via .box_feasible());
# the reported rho* thresholds are the point closed forms, for the narrative.
# -----------------------------------------------------------------------------

#' Disattenuate a correlation matrix for measurement (un)reliability
#'
#' Apply Spearman's correction for attenuation:
#' `D_ij = R_ij / sqrt(rho_i * rho_j)` with the diagonal held at 1.
#'
#' @param R A symmetric correlation matrix with unit diagonal.
#' @param reliability A single reliability applied to all variables, or a
#'   length-`p` vector of per-variable reliabilities. Each must be in `(0, 1]`.
#' @return The disattenuated matrix (unit diagonal). Off-diagonals may exceed 1
#'   in magnitude, which is itself diagnostic (see [check_disattenuated_psd()]).
#' @examples
#' R <- matrix(c(1, 0.6, 0.5, 0.6, 1, 0.55, 0.5, 0.55, 1), 3, 3)
#' disattenuate(R, 0.7)
#' @seealso [check_disattenuated_psd()]
#' @export
disattenuate <- function(R, reliability) {
  R <- .validate_corr(R)
  p <- nrow(R)
  rel <- .check_reliability(reliability, p)
  s <- 1 / sqrt(rel)
  D <- outer(s, s) * R
  diag(D) <- 1
  D
}

.check_reliability <- function(reliability, p) {
  if (!is.numeric(reliability) || anyNA(reliability)) {
    stop("`reliability` must be numeric and non-missing.", call. = FALSE)
  }
  if (!(length(reliability) %in% c(1L, p))) {
    stop(
      sprintf("`reliability` must have length 1 or p = %d.", p),
      call. = FALSE
    )
  }
  if (any(reliability <= 0) || any(reliability > 1)) {
    stop("`reliability` values must lie in (0, 1].", call. = FALSE)
  }
  if (length(reliability) == 1L) rep(reliability, p) else reliability
}

# Per-cell disattenuated intervals from the boxed correlations and boxed
# reliabilities (sound interval arithmetic; denominators are strictly positive).
# Returns clipped [lo, hi] for the PSD box, the unclipped minimum achievable
# magnitude per cell (for the range / plausibility flags), and the point centre.
.disattenuated_intervals <- function(
  R,
  reliability,
  delta_R,
  delta_rel,
  eps = 1e-6
) {
  p <- nrow(R)
  rel <- .check_reliability(reliability, p)
  rel_lo <- pmax(rel - delta_rel, eps)
  rel_hi <- pmin(rel + delta_rel, 1)
  s_lo <- 1 / sqrt(rel_hi) # s = 1/sqrt(rho) decreases in rho
  s_hi <- 1 / sqrt(rel_lo)

  off <- upper.tri(R) | lower.tri(R)
  loD <- matrix(0, p, p)
  hiD <- matrix(0, p, p)
  min_abs <- matrix(0, p, p)
  for (i in seq_len(p - 1L)) {
    for (j in (i + 1L):p) {
      ss_lo <- s_lo[i] * s_lo[j]
      ss_hi <- s_hi[i] * s_hi[j] # positive
      rlo <- R[i, j] - delta_R
      rhi <- R[i, j] + delta_R
      prods <- c(ss_lo * rlo, ss_lo * rhi, ss_hi * rlo, ss_hi * rhi)
      dlo <- min(prods)
      dhi <- max(prods)
      ma <- if (dlo <= 0 && dhi >= 0) 0 else min(abs(dlo), abs(dhi))
      loD[i, j] <- loD[j, i] <- max(dlo, -1)
      hiD[i, j] <- hiD[j, i] <- min(dhi, 1)
      min_abs[i, j] <- min_abs[j, i] <- ma
    }
  }
  diag(loD) <- 1
  diag(hiD) <- 1
  list(
    lo = loD,
    hi = hiD,
    off = off,
    min_abs = min_abs,
    center = disattenuate(R, reliability)
  )
}

# Forward verdict on a disattenuated box: inconsistent / consistent but implausible / consistent /
# undecided, plus the offending cells.
.disatten_forward <- function(iv, cutoff) {
  off <- iv$off
  ut <- which(upper.tri(iv$min_abs), arr.ind = TRUE)
  ma <- iv$min_abs[upper.tri(iv$min_abs)]
  # cells certainly out of range (|D| > 1) or certainly consistent but implausible (|D| > cutoff)
  imp_rows <- ut[ma > 1, , drop = FALSE]
  impl_rows <- ut[ma > cutoff & ma <= 1, , drop = FALSE]

  range_inconsistent <- nrow(imp_rows) > 0L
  fe <- .box_feasible(iv$lo, iv$hi, off)
  psd_inconsistent <- identical(fe$status, "infeasible")
  inconsistent <- range_inconsistent || psd_inconsistent

  any_implausible <- any(ma > cutoff) # includes >1 cells, but those imply inconsistent
  verdict <- if (inconsistent) {
    "inconsistent"
  } else if (any_implausible) {
    "consistent but implausible"
  } else if (identical(fe$status, "feasible")) {
    "consistent"
  } else {
    "undecided"
  }

  list(
    verdict = verdict,
    psd_status = fe$status,
    range_inconsistent = range_inconsistent,
    psd_inconsistent = psd_inconsistent,
    inconsistent_cells = imp_rows,
    implausible_cells = impl_rows
  )
}

.new_disattenuation_check <- function(...) {
  structure(list(...), class = "disattenuation_check")
}

# Common-rho point thresholds (the interpretable narrative numbers, computed at
# the exact reported values; the box-sound counterpart is .rho_inconsistent_box).
.rho_thresholds <- function(R, cutoff) {
  lam_min <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
  max_r <- max(abs(R[upper.tri(R)]))
  rho_inconsistent <- max(1 - lam_min, max_r) # below -> inconsistent
  rho_plausible <- if (cutoff < 1) {
    max(1 - lam_min, max_r / cutoff)
  } else {
    rho_inconsistent
  }
  list(
    lambda_min = lam_min,
    max_r = max_r,
    rho_inconsistent = rho_inconsistent,
    rho_plausible = rho_plausible,
    # which constraint sets each boundary
    inconsistent_binds = if ((1 - lam_min) >= max_r) "PSD" else "range",
    plausible_binds = if (cutoff < 1 && (max_r / cutoff) > (1 - lam_min)) {
      "range"
    } else {
      "PSD"
    }
  )
}

# Box-sound critical reliability: the largest common rho at which the
# disattenuated BOX (correlations rounded to delta_R; reliabilities treated as
# exact, since rho is the variable being solved for) is still certifiably
# inconsistent via the rigorous box witness. The claim "the disattenuation is
# inconsistent for any common reliability <= rho*_box" is therefore sound even
# allowing for rounding of the reported correlations -- unlike the point closed
# form, which ignores the rounding box. Located by bisection (the certified-
# inconsistent region is an interval (0, rho*]).
# Returns 0 when no reliability >= floor_rho is certifiably inconsistent
# ("vacuous"), and NA when the observed box is itself inconsistent at rho = 1.
.rho_inconsistent_box <- function(
  R,
  delta_R,
  observed_inconsistent,
  floor_rho = 0.01,
  tol = 1e-4,
  max_it = 30L
) {
  if (isTRUE(observed_inconsistent)) {
    return(NA_real_)
  }
  cert_imp <- function(rho) {
    iv <- .disattenuated_intervals(R, rho, delta_R, delta_rel = 0)
    if (any(iv$min_abs[upper.tri(iv$min_abs)] > 1)) {
      return(TRUE)
    } # range precheck
    isTRUE(.box_inconsistent(iv$lo, iv$hi, iv$off)$inconsistent)
  }
  if (!cert_imp(floor_rho)) {
    return(0)
  }
  lo <- floor_rho
  hi <- 1
  for (it in seq_len(max_it)) {
    if (hi - lo < tol) {
      break
    }
    m <- (lo + hi) / 2
    if (cert_imp(m)) lo <- m else hi <- m
  }
  lo
}

#' Check whether a disattenuated correlation matrix can be valid
#'
#' Given a reported (rounded) correlation matrix and measurement reliabilities,
#' decide whether the reliability-disattenuated construct matrix is
#' `inconsistent` (no PSD matrix fits, or a corrected correlation exceeds 1),
#' `consistent but implausible` (valid, but a corrected correlation exceeds `max_plausible_r`),
#' `consistent`, or `undecided`. The inconsistency side is rigorous and box-aware
#' (it reuses the package's witness / Rump machinery); rounding of both the
#' correlations and the reliabilities is taken into account.
#'
#' With `reliability = NULL` the function instead reports the *critical
#' reliability* thresholds (the reliability below which disattenuation is
#' inconsistent / consistent but implausible) -- useful when reliabilities were not reported.
#'
#' @param R A reported correlation matrix (rounded), symmetric, unit diagonal.
#' @param reliability A single common reliability, a length-`p` vector of
#'   per-variable reliabilities, or `NULL` to report only the critical thresholds.
#'   Per-variable thresholds are a later addition; with a vector, the forward
#'   verdict is returned but the `rho*` thresholds refer to the common-reliability
#'   model.
#' @param decimals Decimals the correlations were rounded to (sets the box on R).
#' @param reliability_decimals Decimals the reliabilities were rounded to (sets
#'   the box on the reliabilities). Use `Inf` to treat reliabilities as exact.
#' @param max_plausible_r Plausibility cutoff on corrected correlations. A
#'   corrected `|D_ij|` above this (but at most 1) is flagged `consistent but implausible`.
#'   Default `1` disables the plausibility level (only the mathematical bound of
#'   1 flags `inconsistent`); set e.g. `0.9` to activate it.
#' @param plausible_floor Optional lowest reliability the measures could
#'   plausibly attain, used only for the narrative when `reliability = NULL`
#'   (a reachability note, never a filter; no default).
#'
#' @return An S3 `disattenuation_check` object; see the print method.
#' @examples
#' R <- matrix(c(1, 0.60, 0.55, 0.60, 1, 0.50, 0.55, 0.50, 1), 3, 3)
#' # reported reliabilities:
#' check_disattenuated_psd(R, reliability = 0.6, decimals = 2, max_plausible_r = 0.9)
#' # reliabilities not reported -> critical thresholds:
#' check_disattenuated_psd(R, reliability = NULL, decimals = 2, max_plausible_r = 0.9)
#' @seealso [disattenuate()], [check_corr_psd()]
#' @export
check_disattenuated_psd <- function(
  R,
  reliability = NULL,
  decimals = 2,
  reliability_decimals = 2,
  max_plausible_r = 1,
  plausible_floor = NULL
) {
  R <- .validate_corr(R)
  p <- nrow(R)
  if (
    !is.numeric(max_plausible_r) ||
      length(max_plausible_r) != 1L ||
      max_plausible_r <= 0 ||
      max_plausible_r > 1
  ) {
    stop("`max_plausible_r` must be a single number in (0, 1].", call. = FALSE)
  }
  delta_R <- 0.5 * 10^(-decimals)
  delta_rel <- if (is.infinite(reliability_decimals)) {
    0
  } else {
    0.5 * 10^(-reliability_decimals)
  }
  cutoff <- max_plausible_r

  observed <- check_corr_psd(R, decimals = decimals)$verdict
  thr <- .rho_thresholds(R, cutoff)
  thr$rho_inconsistent_box <- .rho_inconsistent_box(
    R,
    delta_R,
    observed_inconsistent = identical(observed, "inconsistent")
  )

  common <- list(
    observed_verdict = observed,
    thresholds = thr,
    max_plausible_r = cutoff,
    decimals = decimals,
    reliability_decimals = reliability_decimals,
    plausible_floor = plausible_floor,
    p = p,
    R = R
  )

  if (is.null(reliability)) {
    return(do.call(
      .new_disattenuation_check,
      c(
        list(
          mode = "critical",
          reliability = NULL,
          verdict = NA_character_,
          disattenuated = NULL,
          forward = NULL,
          headroom = NA_real_
        ),
        common
      )
    ))
  }

  rel <- .check_reliability(reliability, p)
  iv <- .disattenuated_intervals(R, rel, delta_R, delta_rel)
  fwd <- .disatten_forward(iv, cutoff)
  reported_rho <- if (length(reliability) == 1L) reliability else NA_real_
  # headroom is measured against the box-sound boundary (falls back to the point
  # closed form only when the box bisection is unavailable)
  rho_anchor <- if (
    is.finite(thr$rho_inconsistent_box %||% NA_real_) &&
      (thr$rho_inconsistent_box %||% 0) > 0
  ) {
    thr$rho_inconsistent_box
  } else {
    thr$rho_inconsistent
  }
  headroom <- if (!is.na(reported_rho)) reported_rho - rho_anchor else NA_real_
  max_disattenuated <- max(abs(iv$center[upper.tri(iv$center)]))

  do.call(
    .new_disattenuation_check,
    c(
      list(
        mode = if (length(reliability) == 1L) "common" else "per_variable",
        reliability = reliability,
        verdict = fwd$verdict,
        disattenuated = iv$center,
        forward = fwd,
        headroom = headroom,
        max_disattenuated = max_disattenuated,
        reported_rho = reported_rho
      ),
      common
    )
  )
}
