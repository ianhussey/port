# -----------------------------------------------------------------------------
# Component B: rounding-robust inconsistent-triple scan (primary cheap localizer).
#
# For a triple (i,j,k) with a = X_ij, b = X_ik, c = X_jk, the 3x3 PSD condition
# reduces to det = 1 + 2abc - a^2 - b^2 - c^2 >= 0. We need a SOUND upper bound
# on max det over the three rounding boxes: if that bound is < 0, the triple is
# inconsistent even accounting for rounding. This path is pure interval arithmetic
# -- no solver, no CVXR.
#
# Sound bound. Use the exact identities
#   det = (1 - b^2)(1 - c^2) - (a - b c)^2
#       = (1 - a^2)(1 - c^2) - (b - a c)^2
#       = (1 - a^2)(1 - b^2) - (c - a b)^2 .
# For each form, upper-bound the (nonnegative) product term by placing each
# factor's variable at the point of its box nearest 0 (maximising 1 - x^2), and
# lower-bound the squared term by the squared interval gap between the freed
# cell's box and the product interval of the other two. Both are sound over-/
# under-estimates, so the difference is a sound upper bound on det. We take the
# tightest (minimum) over the three forms.
# -----------------------------------------------------------------------------

# Point of interval [lo, hi] nearest 0 (maximises 1 - x^2).
.nearest_zero <- function(lo, hi) min(hi, max(lo, 0))

# Interval product [lo1,hi1] * [lo2,hi2].
.interval_prod <- function(lo1, hi1, lo2, hi2) {
  p <- c(lo1 * lo2, lo1 * hi2, hi1 * lo2, hi1 * hi2)
  c(lo = min(p), hi = max(p))
}

# Gap between intervals [alo,ahi] and [plo,phi] (0 if they overlap).
.interval_gap <- function(alo, ahi, plo, phi) max(0, alo - phi, plo - ahi)

# Sound upper bound on max det over the three boxes; box_i = c(lo, hi).
.triple_det_ub <- function(bx_a, bx_b, bx_c) {
  a0 <- .nearest_zero(bx_a[1], bx_a[2]); fa <- 1 - a0 * a0
  b0 <- .nearest_zero(bx_b[1], bx_b[2]); fb <- 1 - b0 * b0
  c0 <- .nearest_zero(bx_c[1], bx_c[2]); fc <- 1 - c0 * c0

  pbc <- .interval_prod(bx_b[1], bx_b[2], bx_c[1], bx_c[2])
  pac <- .interval_prod(bx_a[1], bx_a[2], bx_c[1], bx_c[2])
  pab <- .interval_prod(bx_a[1], bx_a[2], bx_b[1], bx_b[2])

  ub_a <- fb * fc - .interval_gap(bx_a[1], bx_a[2], pbc[1], pbc[2])^2
  ub_b <- fa * fc - .interval_gap(bx_b[1], bx_b[2], pac[1], pac[2])^2
  ub_c <- fa * fb - .interval_gap(bx_c[1], bx_c[2], pab[1], pab[2])^2
  min(ub_a, ub_b, ub_c)
}

# Scan all triples; return a list of inconsistent ones (box_max_det < 0), each a
# list(vars = c(i,j,k), cells = 3x2 matrix of (row,col), box_max_det).
.scan_inconsistent_triples <- function(R, delta) {
  p <- nrow(R)
  bnd <- .box_bounds(R, delta)
  lo <- bnd$lo; hi <- bnd$hi
  out <- list()
  if (p < 3L) return(out)
  combs <- utils::combn(p, 3L)
  for (m in seq_len(ncol(combs))) {
    i <- combs[1, m]; j <- combs[2, m]; k <- combs[3, m]
    bx_a <- c(lo[i, j], hi[i, j])   # X_ij
    bx_b <- c(lo[i, k], hi[i, k])   # X_ik
    bx_c <- c(lo[j, k], hi[j, k])   # X_jk
    ub <- .triple_det_ub(bx_a, bx_b, bx_c)
    if (ub < 0) {
      out[[length(out) + 1L]] <- list(
        vars = c(i, j, k),
        cells = rbind(c(i, j), c(i, k), c(j, k)),
        box_max_det = ub)
    }
  }
  out
}

#' Scan for rounding-robust inconsistent triples
#'
#' For every triple of variables, test whether the induced 3x3 correlation
#' submatrix can be positive semidefinite for *some* assignment within the
#' rounding boxes. Uses a sound interval bound on the 3x3 determinant, so a
#' flagged triple is inconsistent even after accounting for rounding. Pure
#' arithmetic: no solver required. This is component B of [localize_psd_fault()]
#' exposed on its own.
#'
#' @param R A symmetric numeric correlation matrix with unit diagonal.
#' @param decimals,delta Rounding precision, as in [check_corr_psd()].
#' @return A [tibble][tibble::tibble] with one row per inconsistent triple: columns
#'   `i`, `j`, `k`, and `box_max_det` (the sound upper bound on the box-max
#'   determinant; negative means inconsistent). Empty if none.
#' @examples
#' R <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
#' inconsistent_triples(R, decimals = 2)
#' @seealso [localize_psd_fault()]
#' @export
inconsistent_triples <- function(R, decimals = 2, delta = NULL) {
  R <- .validate_corr(R)
  if (is.null(delta)) delta <- 0.5 * 10^(-decimals)
  trs <- .scan_inconsistent_triples(R, delta)
  if (length(trs) == 0L) {
    return(tibble::tibble(i = integer(), j = integer(), k = integer(),
                          box_max_det = numeric()))
  }
  tibble::tibble(
    i = vapply(trs, function(t) t$vars[1], integer(1)),
    j = vapply(trs, function(t) t$vars[2], integer(1)),
    k = vapply(trs, function(t) t$vars[3], integer(1)),
    box_max_det = vapply(trs, function(t) t$box_max_det, numeric(1)))
}
