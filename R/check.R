# -----------------------------------------------------------------------------
# Main entry point: check_corr_psd().
# -----------------------------------------------------------------------------

# S3 constructor for the result object.
.new_corr_psd_check <- function(
  verdict,
  tier,
  delta,
  p,
  witness = NULL,
  margin = NA_real_,
  b_upper = NA_real_,
  certified_matrix = NULL,
  detail = NULL,
  note = NULL,
  certified = NA,
  rounding = "nearest"
) {
  structure(
    list(
      verdict = verdict, # "inconsistent" / "consistent" / "undecided"
      tier = tier, # "precheck" / "witness" / "pocs"
      certified = certified, # TRUE: verdict carries a certificate
      witness = witness, # certificate vector v (or NULL)
      margin = margin, # B = v'Rv + delta*(||v||_1^2 - ||v||_2^2)
      b_upper = b_upper, # margin + rigorous FP slack
      delta = delta, # (max) half-width of the box
      p = p,
      rounding = rounding, # rounding rule the box was built for
      certified_matrix = certified_matrix, # in-box PSD witness matrix (or NULL)
      detail = detail, # tier-specific detail (e.g. precheck entry)
      note = note
    ),
    class = "corr_psd_check"
  )
}

# Summary of the reported off-diagonal correlations: min / max / mean / SD over
# the unique upper-triangle entries, with NA (freed) cells excluded. The SD is
# the sample SD and is NA when fewer than two correlations are present.
.offdiag_stats <- function(R) {
  v <- R[upper.tri(R)]
  v <- v[!is.na(v)]
  n <- length(v)
  list(
    r_min = if (n >= 1L) min(v) else NA_real_,
    r_max = if (n >= 1L) max(v) else NA_real_,
    r_mean = if (n >= 1L) mean(v) else NA_real_,
    r_sd = if (n >= 2L) sqrt(sum((v - mean(v))^2) / (n - 1L)) else NA_real_
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
#'     edge), `witness` (a test vector `v` gives a rigorously-negative upper
#'     bound `B_upper < 0`), or `pocs` (the same bound, for a test vector that
#'     only the alternating-projection escalation found). Always certified.}
#'   \item{`consistent`}{`certified = TRUE` when an explicit in-box matrix was
#'     verified PSD (Rump 2006) -- a rigorous existence certificate.
#'     `certified = FALSE` marks a *presumed* consistent: the witness margin is
#'     comfortably positive but no in-box PSD matrix was exhibited, so the
#'     verdict means only "not shown inconsistent".}
#'   \item{`undecided`}{Neither proof was found: the witness margin sits in
#'     the ambiguous zone `[0, tau]`, no in-box PSD matrix could be certified,
#'     and the alternating-projection escalation produced neither such a matrix
#'     nor a witness. The table may be inconsistent with a witness the search
#'     missed, or consistent only on the boundary of the PSD cone.}
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
#' @section Running time:
#' Most tables take well under a second. When neither the dense witness
#' directions nor a cheap construction settle the verdict, a sparse witness
#' search runs whose cost grows like `p^5`: about 1.5 s at `p = 80`, 4 s at
#' `p = 100` and 20 s at `p = 150`. If a verdict is still `undecided` after
#' that, the alternating-projection witness search can add up to about 17 s at
#' `p = 100`.
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
#' @param tau Non-negative width of the ambiguous zone. When the witness upper
#'   bound lies in `[0, tau]` and no in-box PSD matrix is certified by the cheap
#'   constructions, the result escalates to the self-contained
#'   alternating-projection search, which looks for an in-box PSD matrix and
#'   then for a witness direction, and is reported `undecided` if both fail.
#'   Defaults to `10 *` the (maximum) box half-width.
#' @param rounding The rule that produced the reported values: `"nearest"`
#'   (default; symmetric box of half-width `0.5 * 10^(-decimals)`),
#'   `"truncate"` (toward zero), `"floor"`, or `"ceiling"` (asymmetric boxes of
#'   width `10^(-decimals)`). Using the rule that actually generated the values
#'   is a soundness requirement: a mismatched symmetric box can exclude the true
#'   value and mislabel a valid matrix.
#'
#' @return An object of class `corr_psd_check` (a list) with elements
#'   `verdict`, `tier`, `certified`, `witness`, `margin`, `b_upper`, `delta`
#'   (maximum half-width), `p`, `rounding`, `certified_matrix`, `detail`,
#'   `note`, and a summary of the reported off-diagonal correlations `r_min`,
#'   `r_max`, `r_mean`, `r_sd` (over the unique off-diagonal entries, excluding
#'   the unit diagonal and any `NA` cells; `r_sd` is the sample SD, `NA` when
#'   fewer than two correlations are present). See [certificate()] and the print
#'   method.
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
check_corr_psd <- function(
  R,
  decimals = 2,
  delta = NULL,
  tau = NULL,
  rounding = c("nearest", "truncate", "floor", "ceiling")
) {
  rounding <- match.arg(rounding)
  if (is.data.frame(R)) {
    R <- as.matrix(R)
  }
  has_na <- is.matrix(R) && is.numeric(R) && anyNA(R)
  R <- if (has_na) .validate_corr_na(R) else .validate_corr(R)
  p <- nrow(R)

  .validate_precision(decimals, delta, rounding, p)

  chk_tau <- function(tau) {
    if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau < 0) {
      stop("`tau` must be a single non-negative finite number.", call. = FALSE)
    }
    tau
  }

  uniform <- .is_uniform_box(has_na, decimals, delta, rounding)

  if (!uniform) {
    bx <- .reported_box(R, decimals, delta, rounding)
    if (is.null(tau)) {
      tau <- 10 * bx$half_max
    }
    res <- .check_box(R, bx, chk_tau(tau), rounding)
    res[c("r_min", "r_max", "r_mean", "r_sd")] <- .offdiag_stats(R)
    attr(res, "R") <- R
    attr(res, "uniform_box") <- FALSE
    return(res)
  }

  if (is.null(delta)) {
    delta <- 0.5 * 10^(-decimals)
  }
  if (is.null(tau)) {
    tau <- 10 * delta
  }
  tau <- chk_tau(tau)

  # The tiered logic lives in build(); we attach R to the result at a single
  # exit point so downstream tools (e.g. localize_psd_fault()) can recover it.
  build <- function() {
    # Tier 1: precheck (out-of-range reported entry).
    pc <- .precheck_range(R, delta)
    if (!is.null(pc)) {
      return(.new_corr_psd_check(
        verdict = "inconsistent",
        tier = "precheck",
        delta = delta,
        p = p,
        certified = TRUE,
        detail = pc,
        note = sprintf(
          "Reported entry R[%d,%d] = %s is out of range: |%s| - delta = %.4g > 1.",
          pc$i,
          pc$j,
          format(pc$value),
          format(pc$value),
          abs(pc$value) - delta
        )
      ))
    }

    # Tier 2: witness-vector bound (primary inconsistency path). The dense
    # directions go first. The O(p^5) sparse search runs only if they certify
    # nothing and no construction certifies consistency (in which case no
    # witness can exist).
    w <- .witness_search(R, delta, sparse = FALSE)
    cons <- NULL
    if (is.null(w) || w$B_upper >= 0) {
      cons <- .construct_consistent(R, delta)
      # ponytail: repeats the dense search, negligible next to the sparse one;
      # split .witness_search() if that ever shows up in a profile.
      if (is.null(cons)) w <- .witness_search(R, delta)
    }
    if (!is.null(w) && w$B_upper < 0) {
      return(.new_corr_psd_check(
        verdict = "inconsistent",
        tier = "witness",
        delta = delta,
        p = p,
        certified = TRUE,
        witness = w$v,
        margin = w$M_hat,
        b_upper = w$B_upper
      ))
    }

    b_upper <- if (is.null(w)) NA_real_ else w$B_upper
    margin <- if (is.null(w)) NA_real_ else w$M_hat

    # Consistency: an in-box PSD matrix was exhibited and verified.
    if (!is.null(cons)) {
      return(.new_corr_psd_check(
        verdict = "consistent",
        tier = "witness",
        delta = delta,
        p = p,
        certified = TRUE,
        margin = margin,
        b_upper = b_upper,
        certified_matrix = cons$X,
        note = cons$how
      ))
    }

    # No construction certified. If the witness margin is comfortably positive we
    # report a PRESUMED "consistent" (not shown inconsistent); otherwise we are in
    # the ambiguous zone and escalate.
    if (!is.na(b_upper) && b_upper > tau) {
      return(.new_corr_psd_check(
        verdict = "consistent",
        tier = "witness",
        delta = delta,
        p = p,
        certified = FALSE,
        margin = margin,
        b_upper = b_upper,
        note = paste(
          "presumed: not shown inconsistent (witness margin above tau);",
          "no in-box PSD matrix was certified"
        )
      ))
    }

    # Ambiguous zone: escalate to the self-contained POCS consistency search.
    pocs <- .pocs_consistent(R, delta)
    if (!is.null(pocs)) {
      return(.new_corr_psd_check(
        verdict = "consistent",
        tier = "pocs",
        delta = delta,
        p = p,
        certified = TRUE,
        margin = margin,
        b_upper = b_upper,
        certified_matrix = pocs$X,
        note = pocs$how
      ))
    }

    # Last resort before "undecided": if the projections cannot bring the box
    # and the PSD cone together, try the direction that separates them.
    bb <- .box_bounds(R, delta)
    wp <- .pocs_witness(bb$lo, bb$hi, bb$off, function(v) {
      .witness_bound(R, v, delta)
    })
    if (wp$B_upper < 0) {
      return(.new_corr_psd_check(
        verdict = "inconsistent",
        tier = "pocs",
        delta = delta,
        p = p,
        certified = TRUE,
        witness = wp$v,
        margin = wp$M_hat,
        b_upper = wp$B_upper,
        note = "witness direction found by alternating projections"
      ))
    }

    .new_corr_psd_check(
      verdict = "undecided",
      tier = "pocs",
      delta = delta,
      p = p,
      certified = FALSE,
      margin = margin,
      b_upper = b_upper,
      note = paste(
        "witness margin within the ambiguous zone [0, tau]; neither an",
        "in-box PSD matrix nor a witness was found, including by",
        "alternating projections"
      )
    )
  }

  res <- build()
  res[c("r_min", "r_max", "r_mean", "r_sd")] <- .offdiag_stats(R)
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
    sprintf(
      paste0(
        "%d missing cell(s) freed to [-1, 1]; the verdict holds for ",
        "every possible value of the missing entries."
      ),
      n_na
    )
  } else {
    NULL
  }

  # Tier 1: precheck -- a reported interval entirely outside [-1, 1].
  if (any(bx$empty)) {
    k <- which(bx$empty & upper.tri(R), arr.ind = TRUE)[1L, ]
    i <- as.integer(k[1])
    j <- as.integer(k[2])
    return(.new_corr_psd_check(
      verdict = "inconsistent",
      tier = "precheck",
      delta = bx$half_max,
      p = p,
      certified = TRUE,
      rounding = rounding,
      detail = list(i = i, j = j, value = R[i, j]),
      note = c(
        sprintf(
          "Reported entry R[%d,%d] = %s cannot arise from any value in [-1, 1] under the '%s' rule.",
          i,
          j,
          format(R[i, j]),
          rounding
        ),
        na_note
      )
    ))
  }

  # Tier 2: box witness (primary inconsistency path). As in the uniform path,
  # the sparse search is deferred until POCS has failed to certify consistency.
  imp <- .box_inconsistent(bx$lo, bx$hi, bx$off, sparse = FALSE)
  hit <- NULL
  if (!isTRUE(imp$inconsistent)) {
    hit <- .pocs_feasible(bx$lo, bx$hi, bx$off)
    if (is.null(hit)) imp <- .box_sparse(bx$lo, bx$hi, bx$off, imp)
  }
  if (isTRUE(imp$inconsistent)) {
    return(.new_corr_psd_check(
      verdict = "inconsistent",
      tier = "witness",
      delta = bx$half_max,
      p = p,
      certified = TRUE,
      rounding = rounding,
      witness = imp$witness,
      margin = imp$margin,
      b_upper = imp$b_upper,
      note = na_note
    ))
  }

  # Consistency: POCS found an in-box matrix, Rump-verified.
  if (!is.null(hit)) {
    return(.new_corr_psd_check(
      verdict = "consistent",
      tier = "pocs",
      delta = bx$half_max,
      p = p,
      certified = TRUE,
      rounding = rounding,
      margin = imp$margin,
      b_upper = imp$b_upper,
      certified_matrix = hit$X,
      note = c(
        "alternating projections found an in-box matrix, independently verified PSD",
        na_note
      )
    ))
  }

  if (!is.na(imp$b_upper) && imp$b_upper > tau) {
    return(.new_corr_psd_check(
      verdict = "consistent",
      tier = "witness",
      delta = bx$half_max,
      p = p,
      certified = FALSE,
      rounding = rounding,
      margin = imp$margin,
      b_upper = imp$b_upper,
      note = c(
        paste(
          "presumed: not shown inconsistent (witness margin above tau);",
          "no in-box PSD matrix was certified"
        ),
        na_note
      )
    ))
  }

  # Last resort before "undecided", as in the uniform path.
  wp <- .pocs_witness(bx$lo, bx$hi, bx$off, function(v) {
    .witness_box_bound(bx$lo, bx$hi, bx$off, v)
  })
  if (wp$B_upper < 0) {
    return(.new_corr_psd_check(
      verdict = "inconsistent",
      tier = "pocs",
      delta = bx$half_max,
      p = p,
      certified = TRUE,
      rounding = rounding,
      witness = wp$v,
      margin = wp$M_hat,
      b_upper = wp$B_upper,
      note = c("witness direction found by alternating projections", na_note)
    ))
  }

  .new_corr_psd_check(
    verdict = "undecided",
    tier = "pocs",
    delta = bx$half_max,
    p = p,
    certified = FALSE,
    rounding = rounding,
    margin = imp$margin,
    b_upper = imp$b_upper,
    note = c(
      paste(
        "witness margin within the ambiguous zone [0, tau]; neither an",
        "in-box PSD matrix nor a witness was found, including by",
        "alternating projections"
      ),
      na_note
    )
  )
}

# Whether the rounding box is the uniform one (nearest rounding, one precision
# for every cell, no missing cells), which check_corr_psd() and check_witness()
# treat with the uniform bound rather than the general box.
.is_uniform_box <- function(has_na, decimals, delta, rounding) {
  !has_na &&
    identical(rounding, "nearest") &&
    !is.matrix(decimals) &&
    !is.matrix(delta)
}

#' Check a witness vector against a reported correlation table
#'
#' Evaluates the bound behind an `inconsistent` verdict for one given vector of
#' weights, without any search. Use it to check the witness that
#' [check_corr_psd()] reports (see [certificate()]), a witness from another
#' source, or one constructed by hand.
#'
#' For weights `v`, the bound is the largest value that `v'Xv` can take over
#' all matrices `X` whose entries lie within the rounding ranges of the
#' reported table, plus a rigorous allowance for floating-point error. `v'Xv`
#' is the variance of the weighted sum of the standardized variables, so a
#' negative bound proves that no dataset can have produced the table. A
#' non-negative bound proves nothing: another vector may still give a negative
#' one.
#'
#' The rounding ranges are those of [check_corr_psd()], and so is the bound: for
#' the witness reported by [check_corr_psd()], `b_upper` is the same number.
#' Reported entries outside `[-1, 1]` even at the nearest end of their range are
#' not examined here; [check_corr_psd()] reports them at tier `precheck`.
#'
#' @inheritParams check_corr_psd
#' @param witness A numeric vector of length `p` with finite entries, not all
#'   zero.
#' @return A one-row tibble with columns `inconsistent` (`TRUE` if the witness
#'   proves the table inconsistent), `b_upper` (the rigorous upper bound;
#'   negative proves inconsistency) and `margin` (the same maximum without the
#'   floating-point allowance).
#' @examples
#' R <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
#'
#' # The witness found by check_corr_psd():
#' result <- check_corr_psd(R, decimals = 2)
#' check_witness(R, certificate(result)$witness, decimals = 2)
#'
#' # A witness by hand: variable 1 minus variables 2 and 3.
#' check_witness(R, c(1, -1, -1), decimals = 2)
#'
#' @seealso [check_corr_psd()], [certificate()]
#' @export
check_witness <- function(
  R,
  witness,
  decimals = 2,
  delta = NULL,
  rounding = c("nearest", "truncate", "floor", "ceiling")
) {
  rounding <- match.arg(rounding)
  if (is.data.frame(R)) {
    R <- as.matrix(R)
  }
  has_na <- is.matrix(R) && is.numeric(R) && anyNA(R)
  R <- if (has_na) .validate_corr_na(R) else .validate_corr(R)
  p <- nrow(R)
  .validate_precision(decimals, delta, rounding, p)
  if (
    !is.numeric(witness) ||
      length(witness) != p ||
      !all(is.finite(witness)) ||
      all(witness == 0)
  ) {
    stop(
      sprintf(
        "`witness` must be a numeric vector of length %d, finite and not all zero.",
        p
      ),
      call. = FALSE
    )
  }
  witness <- as.numeric(witness)

  bound <- if (.is_uniform_box(has_na, decimals, delta, rounding)) {
    if (is.null(delta)) {
      delta <- 0.5 * 10^(-decimals)
    }
    .witness_bound(R, witness, delta)
  } else {
    box <- .reported_box(R, decimals, delta, rounding)
    .witness_box_bound(box$lo, box$hi, box$off, witness)
  }
  tibble::tibble(
    inconsistent = bound$B_upper < 0,
    b_upper = bound$B_upper,
    margin = bound$M_hat
  )
}
