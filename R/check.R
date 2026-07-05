# -----------------------------------------------------------------------------
# Main entry point: check_corr_psd().
# -----------------------------------------------------------------------------

# S3 constructor for the result object.
.new_corr_psd_check <- function(verdict, tier, delta, p,
                                witness = NULL, margin = NA_real_,
                                b_upper = NA_real_, certified_matrix = NULL,
                                detail = NULL, note = NULL,
                                certified = NA, rounding = "nearest") {
  structure(
    list(
      verdict = verdict,                 # "inconsistent" / "consistent" / "undecided"
      tier = tier,                       # "precheck" / "witness" / "pocs"
      certified = certified,             # TRUE: verdict carries a certificate
      witness = witness,                 # certificate vector v (or NULL)
      margin = margin,                   # B = v'Rv + delta*(||v||_1^2 - ||v||_2^2)
      b_upper = b_upper,                 # margin + rigorous FP slack
      delta = delta,                     # (max) half-width of the box
      p = p,
      rounding = rounding,               # rounding rule the box was built for
      certified_matrix = certified_matrix, # in-box PSD witness matrix (or NULL)
      detail = detail,                   # tier-specific detail (e.g. precheck entry)
      note = note
    ),
    class = "corr_psd_check"
  )
}

#' Decide whether a rounded correlation matrix can be positive semidefinite
#'
#' Given a reported correlation matrix `R` whose off-diagonal entries are
#' rounded to `decimals` decimal places, decide whether *any* positive
#' semidefinite (PSD) matrix is consistent with the induced rounding box. The
#' box is `X_ii = 1` exactly and, for `i != j`,
#' `X_ij` in `[R_ij - delta, R_ij + delta]` intersected with `[-1, 1]`, where
#' `delta = 0.5 * 10^(-decimals)` unless a `delta` is supplied directly.
#'
#' The method is a verified witness-vector bound. For any test direction `v`,
#' the most PSD-favourable in-box matrix along `v` attains
#' `max_{X in box} v'Xv = v'Rv + delta*(||v||_1^2 - ||v||_2^2)`. If this maximum
#' is negative, no in-box matrix is PSD and `v` certifies inconsistency. The
#' evaluation of this quantity is bounded by a rigorous a-priori floating-point
#' error model (Higham 2002, Ch. 3), so the inconsistency verdict does not
#' depend on any solver tolerance.
#'
#' @section Verdicts, tiers and certification:
#' \describe{
#'   \item{`inconsistent`}{No PSD matrix exists in the rounding box. Reached at
#'     tier `precheck` (a reported entry is out of range even at its nearest box
#'     edge) or `witness` (a test vector `v` gives a rigorously-negative upper
#'     bound `B_upper < 0`). Always certified.}
#'   \item{`consistent`}{`certified = TRUE` when an explicit in-box matrix was
#'     verified PSD (Rump 2006) -- a rigorous existence certificate.
#'     `certified = FALSE` marks a *presumed* consistent: the witness margin is
#'     comfortably positive but no in-box PSD matrix was exhibited, so the
#'     verdict means only "not shown inconsistent".}
#'   \item{`undecided`}{The witness margin sits in the precision-limited
#'     ambiguous zone `[0, tau]` and no in-box PSD matrix could be certified,
#'     including by the alternating-projection escalation.}
#' }
#'
#' @section Rigor guarantee (and what it does not prove):
#' If `verdict == "inconsistent"`, then no PSD matrix exists in the rounding box,
#' subject to the stated floating-point error model (IEEE double precision,
#' round-to-nearest, unit roundoff `u = .Machine$double.eps / 2`, with Higham
#' `gamma_n` backward/forward error bounds applied conservatively). The verdict
#' does not rely on solver tolerances. Equivalently: the reported values cannot
#' be the **complete-data, single-sample Pearson** correlation matrix of any
#' dataset, up to the stated rounding. The inference from there to a reporting
#' error requires ruling out legitimate generators of non-PSD reported matrices
#' -- pairwise deletion (per-cell n), polychoric/tetrachoric estimation,
#' meta-analytically assembled cells, or upstream disattenuation; large
#' severities (see [localize_psd_fault()]) are hard to explain by any of these.
#' The complementary guarantee also holds: because the bound is a sound
#' over-estimate, a rounding box that contains any PSD matrix can never be
#' reported `inconsistent`.
#'
#' @param R A symmetric numeric `p x p` matrix with unit diagonal (the reported
#'   correlation matrix). Small asymmetries are symmetrized with a warning;
#'   large ones error. `p >= 2` is required. Off-diagonal `NA` entries
#'   (symmetric) are allowed and are freed to `[-1, 1]`: verdicts then hold for
#'   *every* possible value of the missing cells (an `inconsistent` is the
#'   stronger claim "inconsistent whatever the unreported values were").
#' @param decimals Decimal places the off-diagonals were rounded to: a single
#'   integer, or a `p x p` symmetric matrix for mixed-precision tables.
#' @param delta Optional absolute half-width(s) of the rounding box: a scalar or
#'   a `p x p` symmetric matrix. Overrides `decimals`; requires
#'   `rounding = "nearest"`.
#' @param tau Non-negative width of the precision-limited ambiguous zone. When
#'   the witness upper bound lies in `[0, tau]` and no in-box PSD matrix is
#'   certified by the cheap constructions, the result escalates to the
#'   self-contained alternating-projection search, and is reported `undecided`
#'   if that too fails. Defaults to `10 *` the (maximum) box half-width.
#' @param rounding The rule that produced the reported values: `"nearest"`
#'   (default; symmetric box of half-width `0.5 * 10^(-decimals)`),
#'   `"truncate"` (toward zero), `"floor"`, or `"ceiling"` (asymmetric boxes of
#'   width `10^(-decimals)`). Using the rule that actually generated the values
#'   is a soundness requirement: a mismatched symmetric box can exclude the true
#'   value and mislabel a valid matrix.
#'
#' @return An object of class `corr_psd_check` (a list) with elements
#'   `verdict`, `tier`, `certified`, `witness`, `margin`, `b_upper`, `delta`
#'   (maximum half-width), `p`, `rounding`, `certified_matrix`, `detail` and
#'   `note`. See [certificate()] and the print method.
#'
#' @examples
#' # Strongly inconsistent 3x3: r12 = r13 = 0.9, r23 = -0.9.
#' R <- matrix(c(1, 0.9, 0.9,
#'               0.9, 1, -0.9,
#'               0.9, -0.9, 1), 3, 3)
#' check_corr_psd(R, decimals = 2)
#'
#' # A genuine correlation matrix rounded to 2dp is never inconsistent.
#' G <- round(cor(matrix(rnorm(200), ncol = 4)), 2)
#' check_corr_psd(G, decimals = 2)
#'
#' @seealso [certificate()], [check_corr_psd_batch()]
#' @export
check_corr_psd <- function(R, decimals = 2, delta = NULL, tau = NULL,
                           rounding = c("nearest", "truncate", "floor", "ceiling")) {
  rounding <- match.arg(rounding)
  if (is.data.frame(R)) R <- as.matrix(R)
  has_na <- is.matrix(R) && is.numeric(R) && anyNA(R)
  R <- if (has_na) .validate_corr_na(R) else .validate_corr(R)
  p <- nrow(R)

  # --- precision arguments ---------------------------------------------------
  dec_ok <- function(x) is.numeric(x) && !anyNA(x) && all(is.finite(x)) &&
    all(x >= 0) && all(x == round(x))
  if (is.null(delta)) {
    if (is.matrix(decimals)) {
      if (!all(dim(decimals) == p) || !dec_ok(decimals) ||
          max(abs(decimals - t(decimals))) > 0) {
        stop("A `decimals` matrix must be p x p, symmetric, with non-negative ",
             "integer entries.", call. = FALSE)
      }
    } else if (!(length(decimals) == 1L && dec_ok(decimals))) {
      stop("`decimals` must be a single non-negative integer (or a p x p matrix).",
           call. = FALSE)
    }
  } else {
    if (!identical(rounding, "nearest")) {
      stop("Supply `decimals` (not `delta`) when `rounding` is not \"nearest\": ",
           "asymmetric boxes are derived from the decimal width.", call. = FALSE)
    }
    if (is.matrix(delta)) {
      if (!all(dim(delta) == p) || anyNA(delta) || any(!is.finite(delta)) ||
          any(delta < 0) || max(abs(delta - t(delta))) > 0) {
        stop("A `delta` matrix must be p x p, symmetric, non-negative and finite.",
             call. = FALSE)
      }
    } else if (!(is.numeric(delta) && length(delta) == 1L && is.finite(delta) &&
                 delta >= 0)) {
      stop("`delta` must be a single non-negative finite number (or a p x p matrix).",
           call. = FALSE)
    }
  }
  chk_tau <- function(tau) {
    if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau < 0) {
      stop("`tau` must be a single non-negative finite number.", call. = FALSE)
    }
    tau
  }

  uniform <- !has_na && identical(rounding, "nearest") &&
    !is.matrix(decimals) && !is.matrix(delta)

  if (!uniform) {
    bx <- .reported_box(R, decimals, delta, rounding)
    if (is.null(tau)) tau <- 10 * bx$half_max
    res <- .check_box(R, bx, chk_tau(tau), rounding)
    attr(res, "R") <- R
    attr(res, "uniform_box") <- FALSE
    return(res)
  }

  if (is.null(delta)) delta <- 0.5 * 10^(-decimals)
  if (is.null(tau)) tau <- 10 * delta
  tau <- chk_tau(tau)

  # The tiered logic lives in build(); we attach R to the result at a single
  # exit point so downstream tools (e.g. localize_psd_fault()) can recover it.
  build <- function() {
    # Tier 1: precheck (out-of-range reported entry).
    pc <- .precheck_range(R, delta)
    if (!is.null(pc)) {
      return(.new_corr_psd_check(
        verdict = "inconsistent", tier = "precheck", delta = delta, p = p,
        certified = TRUE, detail = pc,
        note = sprintf(
          "Reported entry R[%d,%d] = %s is out of range: |%s| - delta = %.4g > 1.",
          pc$i, pc$j, format(pc$value), format(pc$value), abs(pc$value) - delta)))
    }

    # Tier 2: witness-vector bound (primary inconsistency path).
    w <- .witness_search(R, delta)
    if (!is.null(w) && w$B_upper < 0) {
      return(.new_corr_psd_check(
        verdict = "inconsistent", tier = "witness", delta = delta, p = p,
        certified = TRUE, witness = w$v, margin = w$M_hat, b_upper = w$B_upper))
    }

    b_upper <- if (is.null(w)) NA_real_ else w$B_upper
    margin  <- if (is.null(w)) NA_real_ else w$M_hat

    # Consistency: try to exhibit and verify an in-box PSD matrix.
    cons <- .construct_consistent(R, delta)
    if (!is.null(cons)) {
      return(.new_corr_psd_check(
        verdict = "consistent", tier = "witness", delta = delta, p = p,
        certified = TRUE, margin = margin, b_upper = b_upper,
        certified_matrix = cons$X, note = cons$how))
    }

    # No construction certified. If the witness margin is comfortably positive we
    # report a PRESUMED "consistent" (not shown inconsistent); otherwise we are in
    # the ambiguous zone and escalate.
    if (!is.na(b_upper) && b_upper > tau) {
      return(.new_corr_psd_check(
        verdict = "consistent", tier = "witness", delta = delta, p = p,
        certified = FALSE, margin = margin, b_upper = b_upper,
        note = paste("presumed: not shown inconsistent (witness margin above tau);",
                     "no in-box PSD matrix was certified")))
    }

    # Ambiguous zone: escalate to the self-contained POCS consistency search.
    pocs <- .pocs_consistent(R, delta)
    if (!is.null(pocs)) {
      return(.new_corr_psd_check(
        verdict = "consistent", tier = "pocs", delta = delta, p = p,
        certified = TRUE, margin = margin, b_upper = b_upper,
        certified_matrix = pocs$X, note = pocs$how))
    }

    .new_corr_psd_check(
      verdict = "undecided", tier = "pocs", delta = delta, p = p,
      certified = FALSE, margin = margin, b_upper = b_upper,
      note = paste("witness margin within the precision-limited zone [0, tau]",
                   "and no in-box PSD matrix could be constructed or found by",
                   "alternating projections"))
  }

  res <- build()
  attr(res, "R") <- R
  attr(res, "uniform_box") <- TRUE
  res
}

# General (heterogeneous) box path: asymmetric rounding rules, per-cell
# precision, and NA cells freed to [-1, 1]. Mirrors the uniform tiers on the
# generalized box machinery; all certificates are as rigorous as the uniform
# path (witness bound + Rump-verified POCS point).
.check_box <- function(R, bx, tau, rounding) {
  p <- nrow(R)
  n_na <- sum(bx$na_mask & upper.tri(R))
  na_note <- if (n_na > 0L) {
    sprintf(paste0("%d missing cell(s) freed to [-1, 1]; the verdict holds for ",
                   "every possible value of the missing entries."), n_na)
  } else NULL

  # Tier 1: precheck -- a reported interval entirely outside [-1, 1].
  if (any(bx$empty)) {
    k <- which(bx$empty & upper.tri(R), arr.ind = TRUE)[1L, ]
    i <- as.integer(k[1]); j <- as.integer(k[2])
    return(.new_corr_psd_check(
      verdict = "inconsistent", tier = "precheck", delta = bx$half_max, p = p,
      certified = TRUE, rounding = rounding,
      detail = list(i = i, j = j, value = R[i, j]),
      note = c(sprintf(
        "Reported entry R[%d,%d] = %s cannot arise from any value in [-1, 1] under the '%s' rule.",
        i, j, format(R[i, j]), rounding), na_note)))
  }

  # Tier 2: box witness (primary inconsistency path).
  imp <- .box_inconsistent(bx$lo, bx$hi, bx$off)
  if (isTRUE(imp$inconsistent)) {
    return(.new_corr_psd_check(
      verdict = "inconsistent", tier = "witness", delta = bx$half_max, p = p,
      certified = TRUE, rounding = rounding,
      witness = imp$witness, margin = imp$margin, b_upper = imp$b_upper,
      note = na_note))
  }

  # Consistency: POCS search + Rump verification on the heterogeneous box.
  hit <- .pocs_feasible(bx$lo, bx$hi, bx$off)
  if (!is.null(hit)) {
    return(.new_corr_psd_check(
      verdict = "consistent", tier = "pocs", delta = bx$half_max, p = p,
      certified = TRUE, rounding = rounding,
      margin = imp$margin, b_upper = imp$b_upper, certified_matrix = hit$X,
      note = c("alternating projections found an in-box matrix, independently verified PSD",
               na_note)))
  }

  if (!is.na(imp$b_upper) && imp$b_upper > tau) {
    return(.new_corr_psd_check(
      verdict = "consistent", tier = "witness", delta = bx$half_max, p = p,
      certified = FALSE, rounding = rounding,
      margin = imp$margin, b_upper = imp$b_upper,
      note = c(paste("presumed: not shown inconsistent (witness margin above tau);",
                     "no in-box PSD matrix was certified"), na_note)))
  }

  .new_corr_psd_check(
    verdict = "undecided", tier = "pocs", delta = bx$half_max, p = p,
    certified = FALSE, rounding = rounding,
    margin = imp$margin, b_upper = imp$b_upper,
    note = c(paste("witness margin within the precision-limited zone [0, tau]",
                   "and no in-box PSD matrix could be found by alternating",
                   "projections"), na_note))
}
