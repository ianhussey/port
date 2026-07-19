# -----------------------------------------------------------------------------
# Input validation and symmetrization.
# -----------------------------------------------------------------------------

# Validate a reported correlation matrix and return a cleaned, symmetric copy.
#
# Rules (see ?check_corr_psd):
#   * must be a numeric, square matrix with p >= 2;
#   * no NA / non-finite entries;
#   * unit diagonal (within a tight tolerance);
#   * non-symmetric input is symmetrized (with a warning) if the asymmetry is
#     small, and errors if the asymmetry is large.
#
# Returns the symmetric numeric matrix.
.validate_corr <- function(
  R,
  sym_warn = 1e-8,
  sym_error = 1e-4,
  diag_tol = 1e-8
) {
  if (is.data.frame(R)) {
    R <- as.matrix(R)
  }
  if (!is.matrix(R) || !is.numeric(R)) {
    stop(
      "`R` must be a numeric matrix (or a data frame coercible to one).",
      call. = FALSE
    )
  }
  if (nrow(R) != ncol(R)) {
    stop(
      sprintf("`R` must be square; got %d x %d.", nrow(R), ncol(R)),
      call. = FALSE
    )
  }
  p <- nrow(R)
  if (p < 2L) {
    stop("`R` must have at least 2 variables (p >= 2).", call. = FALSE)
  }
  if (anyNA(R) || any(!is.finite(R))) {
    stop(
      "`R` contains NA or non-finite entries; please supply a complete, ",
      "finite matrix.",
      call. = FALSE
    )
  }

  # Diagonal must be (numerically) 1.
  d <- diag(R)
  if (any(abs(d - 1) > diag_tol)) {
    bad <- which(abs(d - 1) > diag_tol)
    stop(
      sprintf(
        "Diagonal entries must equal 1; entry/entries %s = %s.",
        paste(bad, collapse = ", "),
        paste(format(d[bad]), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Symmetry.
  asym <- max(abs(R - t(R)))
  if (asym > sym_error) {
    stop(
      sprintf(
        paste0(
          "`R` is not symmetric (max |R - t(R)| = %.3g > %.3g). ",
          "Refusing to symmetrize a strongly asymmetric matrix."
        ),
        asym,
        sym_error
      ),
      call. = FALSE
    )
  }
  if (asym > sym_warn) {
    warning(
      sprintf(
        paste0(
          "`R` is not exactly symmetric (max |R - t(R)| = %.3g); ",
          "symmetrizing as (R + t(R)) / 2."
        ),
        asym
      ),
      call. = FALSE
    )
  }
  R <- (R + t(R)) / 2
  diag(R) <- 1 # pin diagonal to exactly 1 after averaging
  dimnames(R) <- NULL
  R
}
