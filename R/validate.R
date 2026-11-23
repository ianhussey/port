# -----------------------------------------------------------------------------
# Input validation and symmetrization.
# -----------------------------------------------------------------------------

# Validate a reported correlation matrix and return a cleaned, symmetric copy.
#
# Rules (see ?check_corr_psd):
#   * must be a numeric, square matrix with p >= 2;
#   * no infinite entries; no NA entries unless `allow_na`, and then only
#     off the diagonal and mirrored across it;
#   * unit diagonal (within a tight tolerance);
#   * non-symmetric input is symmetrized (with a warning) if the asymmetry is
#     small, and errors if the asymmetry is large.
#
# Returns the symmetric numeric matrix, without dimnames.
.validate_corr <- function(
  R,
  allow_na = FALSE,
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
  if (any(is.infinite(R))) {
    stop("`R` contains infinite entries.", call. = FALSE)
  }
  if (!allow_na && anyNA(R)) {
    stop(
      "`R` contains NA entries; please supply a complete matrix.",
      call. = FALSE
    )
  }
  if (allow_na && !all(is.na(R) == is.na(t(R)))) {
    stop("Missing (NA) entries must be symmetric.", call. = FALSE)
  }

  # Diagonal must be (numerically) 1.
  d <- diag(R)
  if (anyNA(d) || any(abs(d - 1) > diag_tol)) {
    bad <- which(is.na(d) | abs(d - 1) > diag_tol)
    stop(
      sprintf(
        "Diagonal entries must equal 1 (and be non-missing); entry/entries %s = %s.",
        paste(bad, collapse = ", "),
        paste(format(d[bad]), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Symmetry (over the reported cells).
  asym <- max(0, abs(R - t(R)), na.rm = TRUE)
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

# Validate a single `decimals` or `delta` for the functions that support only
# a uniform nearest-rounding box, and return the half-width.
.uniform_delta <- function(decimals, delta, p) {
  if (is.matrix(decimals) || is.matrix(delta)) {
    stop(
      "Only a single `decimals` or `delta` is supported here (no per-cell ",
      "precision).",
      call. = FALSE
    )
  }
  .validate_precision(decimals, delta, "nearest", p)
  if (is.null(delta)) 0.5 * 10^(-decimals) else delta
}

# Validate the precision arguments shared by check_corr_psd() and
# check_witness(): `decimals` or `delta` (scalar or p x p matrix) and the
# rounding rule. Errors on invalid input; returns nothing.
.validate_precision <- function(decimals, delta, rounding, p) {
  dec_ok <- function(x) {
    is.numeric(x) &&
      !anyNA(x) &&
      all(is.finite(x)) &&
      all(x >= 0) &&
      all(x == round(x))
  }
  if (is.null(delta)) {
    if (is.matrix(decimals)) {
      if (
        !all(dim(decimals) == p) ||
          !dec_ok(decimals) ||
          max(abs(decimals - t(decimals))) > 0
      ) {
        stop(
          "A `decimals` matrix must be p x p, symmetric, with non-negative ",
          "integer entries.",
          call. = FALSE
        )
      }
    } else if (!(length(decimals) == 1L && dec_ok(decimals))) {
      stop(
        "`decimals` must be a single non-negative integer (or a p x p matrix).",
        call. = FALSE
      )
    }
  } else {
    if (!identical(rounding, "nearest")) {
      stop(
        "Supply `decimals` (not `delta`) when `rounding` is not \"nearest\": ",
        "asymmetric boxes are derived from the decimal width.",
        call. = FALSE
      )
    }
    if (is.matrix(delta)) {
      if (
        !all(dim(delta) == p) ||
          anyNA(delta) ||
          any(!is.finite(delta)) ||
          any(delta < 0) ||
          max(abs(delta - t(delta))) > 0
      ) {
        stop(
          "A `delta` matrix must be p x p, symmetric, non-negative and finite.",
          call. = FALSE
        )
      }
    } else if (
      !(is.numeric(delta) &&
        length(delta) == 1L &&
        is.finite(delta) &&
        delta >= 0)
    ) {
      stop(
        "`delta` must be a single non-negative finite number (or a p x p matrix).",
        call. = FALSE
      )
    }
  }
  invisible()
}
