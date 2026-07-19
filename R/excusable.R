# -----------------------------------------------------------------------------
# The excusable-imprecision threshold: one headline number in precision units.
# -----------------------------------------------------------------------------

#' Smallest reporting imprecision that could excuse a matrix
#'
#' The excusable delta of a reported correlation matrix is the smallest uniform
#' rounding half-width `w` at which the box `[R_ij - w, R_ij + w]` (unit
#' diagonal, clipped to `[-1, 1]`) admits a positive semidefinite matrix. It
#' unifies the verdict and the severity into a single interpretable statistic in
#' the units of reporting precision:
#'
#' * the matrix is inconsistent at reporting precision `delta` **iff**
#'   `excusable_delta(R) > delta`;
#' * a value like `0.04` reads as "the reported values are inconsistent with any
#'   correlation matrix unless each entry were misstated by more than +/-0.04 --
#'   not excusable even by rounding to 1 decimal place (+/-0.05 would be needed)".
#'
#' For a matrix that is already consistent at fine precision the threshold is ~0.
#' The value is located by bisection over the box feasibility oracle (witness
#' refutation + Rump-verified POCS), to tolerance `tol`.
#'
#' @param R A symmetric numeric correlation matrix with unit diagonal.
#' @param tol Bisection tolerance for the threshold.
#' @return A single number: the smallest excusing half-width `w`.
#' @examples
#' R <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
#' excusable_delta(R)   # far above 0.005: 2dp rounding cannot excuse it
#' excusable_delta(diag(3))   # ~0
#' @seealso [check_corr_psd()], [localize_psd_fault()] (whose severity carries
#'   the same quantity as `excusable_delta`)
#' @export
excusable_delta <- function(R, tol = 1e-6) {
  R <- .validate_corr(R)
  .severity_max(R, delta = 0, tol = tol)
}

# Translate an excusable half-width into the coarsest decimal precision whose
# nearest-rounding could excuse it: the largest d with 0.5*10^(-d) >= w.
# Returns -Inf when even 0-decimal rounding (+-0.5) could not excuse it.
.excusable_decimals <- function(w) {
  if (!is.finite(w) || w <= 0) {
    return(Inf)
  }
  d <- floor(-log10(2 * w) + 1e-12)
  if (0.5 * 10^(-d) < w) {
    d <- d - 1
  } # guard against boundary rounding
  d
}
