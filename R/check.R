# -----------------------------------------------------------------------------
# Main entry point: check_corr_psd().
# -----------------------------------------------------------------------------

# S3 constructor for the result object.
.new_corr_psd_check <- function(verdict, tier, delta, p,
                                witness = NULL, margin = NA_real_,
                                b_upper = NA_real_, certified_matrix = NULL,
                                detail = NULL, note = NULL) {
  structure(
    list(
      verdict = verdict,                 # "impossible" / "possible" / "undecided"
      tier = tier,                       # "precheck" / "witness" / "pocs"
      witness = witness,                 # certificate vector v (or NULL)
      margin = margin,                   # B = v'Rv + delta*(||v||_1^2 - ||v||_2^2)
      b_upper = b_upper,                 # margin + rigorous FP slack
      delta = delta,
      p = p,
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
#' is negative, no in-box matrix is PSD and `v` certifies impossibility. The
#' evaluation of this quantity is bounded by a rigorous a-priori floating-point
#' error model (Higham 2002, Ch. 3), so the impossibility verdict does not
#' depend on any solver tolerance.
#'
#' @section Verdicts and tiers:
#' \describe{
#'   \item{`impossible`}{No PSD matrix exists in the rounding box. Reached at
#'     tier `precheck` (a reported entry is out of range even at its nearest box
#'     edge) or `witness` (a test vector `v` gives a rigorously-negative upper
#'     bound `B_upper < 0`).}
#'   \item{`possible`}{No impossibility certificate was found. When
#'     `certified_matrix` is non-`NULL`, an explicit in-box matrix was verified
#'     PSD (Rump 2006), a strong possibility certificate; otherwise the verdict
#'     means only "not shown impossible".}
#'   \item{`undecided`}{The witness margin sits in the precision-limited
#'     ambiguous zone `[0, tau]` and no in-box PSD matrix could be certified,
#'     including by the alternating-projection escalation.}
#' }
#'
#' @section Rigor guarantee:
#' If `verdict == "impossible"`, then no PSD matrix exists in the rounding box,
#' subject to the stated floating-point error model (IEEE double precision,
#' round-to-nearest, unit roundoff `u = .Machine$double.eps / 2`, with Higham
#' `gamma_n` backward/forward error bounds applied conservatively). The verdict
#' does not rely on solver tolerances. The complementary guarantee also holds:
#' because the bound is a sound over-estimate, a rounding box that contains any
#' PSD matrix can never be reported `impossible`.
#'
#' @param R A symmetric numeric `p x p` matrix with unit diagonal (the reported
#'   correlation matrix). Small asymmetries are symmetrized with a warning;
#'   large ones error. `p >= 2` is required.
#' @param decimals Number of decimal places the off-diagonals were rounded to.
#'   Used to derive `delta = 0.5 * 10^(-decimals)` when `delta` is `NULL`.
#' @param delta Optional absolute half-width of the rounding box. Overrides
#'   `decimals` when supplied.
#' @param tau Non-negative width of the precision-limited ambiguous zone. When
#'   the witness upper bound lies in `[0, tau]` and no in-box PSD matrix is
#'   certified by the cheap constructions, the result escalates to the
#'   self-contained alternating-projection search, and is reported `undecided`
#'   if that too fails. Defaults to `10 * delta`.
#'
#' @return An object of class `corr_psd_check` (a list) with elements
#'   `verdict`, `tier`, `witness`, `margin`, `b_upper`, `delta`, `p`,
#'   `certified_matrix`, `detail` and `note`. See [certificate()] and the print
#'   method.
#'
#' @examples
#' # Strongly inconsistent 3x3: r12 = r13 = 0.9, r23 = -0.9.
#' R <- matrix(c(1, 0.9, 0.9,
#'               0.9, 1, -0.9,
#'               0.9, -0.9, 1), 3, 3)
#' check_corr_psd(R, decimals = 2)
#'
#' # A genuine correlation matrix rounded to 2dp is never impossible.
#' G <- round(cor(matrix(rnorm(200), ncol = 4)), 2)
#' check_corr_psd(G, decimals = 2)
#'
#' @seealso [certificate()], [check_corr_psd_batch()]
#' @export
check_corr_psd <- function(R, decimals = 2, delta = NULL, tau = NULL) {
  R <- .validate_corr(R)
  p <- nrow(R)

  if (is.null(delta)) {
    if (!is.numeric(decimals) || length(decimals) != 1L || decimals < 0 ||
        decimals != round(decimals)) {
      stop("`decimals` must be a single non-negative integer.", call. = FALSE)
    }
    delta <- 0.5 * 10^(-decimals)
  }
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta < 0) {
    stop("`delta` must be a single non-negative finite number.", call. = FALSE)
  }
  if (is.null(tau)) tau <- 10 * delta
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau < 0) {
    stop("`tau` must be a single non-negative finite number.", call. = FALSE)
  }

  # The tiered logic lives in build(); we attach R to the result at a single
  # exit point so downstream tools (e.g. localize_psd_fault()) can recover it.
  build <- function() {
    # Tier 1: precheck (out-of-range reported entry).
    pc <- .precheck_range(R, delta)
    if (!is.null(pc)) {
      return(.new_corr_psd_check(
        verdict = "impossible", tier = "precheck", delta = delta, p = p,
        detail = pc,
        note = sprintf(
          "Reported entry R[%d,%d] = %s is out of range: |%s| - delta = %.4g > 1.",
          pc$i, pc$j, format(pc$value), format(pc$value), abs(pc$value) - delta)))
    }

    # Tier 2: witness-vector bound (primary impossibility path).
    w <- .witness_search(R, delta)
    if (!is.null(w) && w$B_upper < 0) {
      return(.new_corr_psd_check(
        verdict = "impossible", tier = "witness", delta = delta, p = p,
        witness = w$v, margin = w$M_hat, b_upper = w$B_upper))
    }

    b_upper <- if (is.null(w)) NA_real_ else w$B_upper
    margin  <- if (is.null(w)) NA_real_ else w$M_hat

    # Possibility: try to exhibit and verify an in-box PSD matrix.
    cons <- .construct_possible(R, delta)
    if (!is.null(cons)) {
      return(.new_corr_psd_check(
        verdict = "possible", tier = "witness", delta = delta, p = p,
        margin = margin, b_upper = b_upper,
        certified_matrix = cons$X, note = cons$how))
    }

    # No construction certified. If the witness margin is comfortably positive we
    # report "possible" (not shown impossible); otherwise we are in the ambiguous
    # zone and escalate.
    if (!is.na(b_upper) && b_upper > tau) {
      return(.new_corr_psd_check(
        verdict = "possible", tier = "witness", delta = delta, p = p,
        margin = margin, b_upper = b_upper,
        note = paste("not shown impossible (witness margin above tau);",
                     "no in-box PSD matrix was certified")))
    }

    # Ambiguous zone: escalate to the self-contained POCS possibility search.
    pocs <- .pocs_possible(R, delta)
    if (!is.null(pocs)) {
      return(.new_corr_psd_check(
        verdict = "possible", tier = "pocs", delta = delta, p = p,
        margin = margin, b_upper = b_upper,
        certified_matrix = pocs$X, note = pocs$how))
    }

    .new_corr_psd_check(
      verdict = "undecided", tier = "pocs", delta = delta, p = p,
      margin = margin, b_upper = b_upper,
      note = paste("witness margin within the precision-limited zone [0, tau]",
                   "and no in-box PSD matrix could be constructed or found by",
                   "alternating projections"))
  }

  res <- build()
  attr(res, "R") <- R
  res
}
