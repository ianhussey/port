# -----------------------------------------------------------------------------
# Closed-form per-cell implied intervals (V3) and single-missing-cell imputation
# (V8), from van Tilburg & van Tilburg (2023).
#
# Freeing one entry X_ij and requiring the matrix PSD, with the OTHER cells held
# at their reported points and the complementary submatrix M_{-{i,j}} PD, gives a
# closed-form interval. Taking the Schur complement onto the (i,j) block,
#   X_ij must satisfy  (X_ij - g)^2 <= (1 - alpha)(1 - beta),
# where, with A = M_{-{i,j}}^{-1} and c_i, c_j the two variables' correlations to
# the rest,  alpha = c_i' A c_i,  beta = c_j' A c_j,  g = c_i' A c_j. Hence
#   X_ij in [g - sqrt((1-alpha)(1-beta)), g + sqrt((1-alpha)(1-beta))] intersect [-1, 1].
# This is a fast, point-conditioned (partial-correlation) reading. It is NOT the
# same as the box-aware interval (component A), which lets the other cells roam
# their rounding boxes; that stays the rigorous localization path. Here the
# closed form supplements it and drives missing-cell imputation.
# -----------------------------------------------------------------------------

# Point-conditioned implied interval for cell (i,j); others at their reported
# values. Returns c(lo, hi), or NA endpoints when the complement is not usable
# (M_{-{i,j}} not invertible, or i/j already over-determined by the rest alone).
.implied_interval_points <- function(R, i, j) {
  p <- nrow(R)
  others <- setdiff(seq_len(p), c(i, j))
  if (length(others) == 0L) {
    return(c(lo = -1, hi = 1))
  } # p = 2
  M <- R[others, others, drop = FALSE]
  if (anyNA(M)) {
    return(c(lo = NA_real_, hi = NA_real_))
  } # another cell missing
  Ainv <- tryCatch(solve(M), error = function(e) NULL)
  if (is.null(Ainv)) {
    return(c(lo = NA_real_, hi = NA_real_))
  }
  ci <- R[others, i]
  cj <- R[others, j]
  if (anyNA(ci) || anyNA(cj)) {
    return(c(lo = NA_real_, hi = NA_real_))
  }
  alpha <- sum(ci * (Ainv %*% ci))
  beta <- sum(cj * (Ainv %*% cj))
  g <- sum(ci * (Ainv %*% cj))
  if (1 - alpha < 0 || 1 - beta < 0) {
    return(c(lo = NA_real_, hi = NA_real_))
  } # empty
  rad <- sqrt((1 - alpha) * (1 - beta))
  c(lo = max(-1, g - rad), hi = min(1, g + rad))
}

# Relaxed validation: like .validate_corr but permits NA off-diagonals.
.validate_corr_na <- function(R) {
  if (is.data.frame(R)) {
    R <- as.matrix(R)
  }
  if (!is.matrix(R) || !is.numeric(R)) {
    stop("`R` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(R) != ncol(R)) {
    stop("`R` must be square.", call. = FALSE)
  }
  if (nrow(R) < 2L) {
    stop("`R` must have at least 2 variables.", call. = FALSE)
  }
  d <- diag(R)
  if (anyNA(d) || any(abs(d - 1) > 1e-8)) {
    stop("Diagonal entries must equal 1 (and be non-missing).", call. = FALSE)
  }
  asym <- max(abs(R - t(R)), na.rm = TRUE)
  if (is.finite(asym) && asym > 1e-6) {
    stop("`R` is not symmetric.", call. = FALSE)
  }
  # NA must be mirrored across the diagonal
  na_ok <- all(is.na(R) == is.na(t(R)))
  if (!na_ok) {
    stop("Missing (NA) entries must be symmetric.", call. = FALSE)
  }
  R
}

#' Implied interval for a correlation cell (and single-cell imputation)
#'
#' Compute the interval a correlation cell could occupy for the matrix to be a
#' valid (PSD) correlation matrix, given the other cells. Two uses:
#' interrogating a *reported* cell (how far outside its rounding box it would
#' have to move), and imputing a genuinely *missing* (`NA`) cell (the interval it
#' must have occupied -- van Tilburg & van Tilburg 2023, Use 3).
#'
#' @param R A symmetric numeric correlation matrix with unit diagonal. Up to one
#'   off-diagonal entry per interrogated cell may be `NA` (missing); all other
#'   off-diagonals must be finite.
#' @param cells Optional 2-column matrix (or length-2 vector) of `(i, j)` cells
#'   to interrogate. Defaults to the `NA` off-diagonal cell(s) if any.
#' @param decimals,delta Rounding precision (as in [check_corr_psd()]); used for
#'   the `required_edit` of reported cells and for the box path.
#' @param hold `"points"` (default) holds the other cells at their reported
#'   values -- a fast closed form and the cleanest partial-correlation reading;
#'   `"box"` lets the other cells roam their rounding boxes (rigorous, reuses the
#'   box feasibility oracle).
#' @param verify Passed to the box oracle when `hold = "box"`.
#'
#' @return A [tibble][tibble::tibble] with columns `i`, `j`, `reported`, `lo`,
#'   `hi`, `required_edit` (signed distance from the rounding box to `[lo, hi]`;
#'   `NA` for missing cells), `hold`, and `status` (`"feasible"`, `"empty"`, or
#'   `"na"`).
#' @examples
#' R <- matrix(c(1, 0.5, 0.5, 0.82,
#'               0.5, 1, 0.5, 0.82,
#'               0.5, 0.5, 1, 0.82,
#'               0.82, 0.82, 0.82, 1), 4, 4)
#' R[1, 4] <- R[4, 1] <- NA           # pretend this correlation was unreported
#' implied_interval(R)
#' @seealso [localize_psd_fault()]
#' @export
implied_interval <- function(
  R,
  cells = NULL,
  decimals = 2,
  delta = NULL,
  hold = c("points", "box"),
  verify = TRUE
) {
  hold <- match.arg(hold)
  R <- .validate_corr_na(R)
  p <- nrow(R)
  if (is.null(delta)) {
    delta <- 0.5 * 10^(-decimals)
  }

  if (is.null(cells)) {
    na_idx <- which(is.na(R) & upper.tri(R), arr.ind = TRUE)
    if (nrow(na_idx) == 0L) {
      stop(
        "No missing (NA) cells found; supply `cells` to interrogate ",
        "reported entries.",
        call. = FALSE
      )
    }
    cells <- na_idx
  }
  if (is.null(dim(cells))) {
    cells <- matrix(cells, nrow = 1)
  }
  storage.mode(cells) <- "integer"

  rows <- vector("list", nrow(cells))
  for (r in seq_len(nrow(cells))) {
    i <- cells[r, 1]
    j <- cells[r, 2]
    if (i == j) {
      stop("Cells must be off-diagonal.", call. = FALSE)
    }
    # every other off-diagonal must be finite for this interrogation
    Rtest <- R
    Rtest[i, j] <- Rtest[j, i] <- 0
    if (anyNA(Rtest)) {
      stop(
        sprintf(
          "Interrogating cell (%d,%d) requires all other cells to be ",
          i,
          j
        ),
        "non-missing.",
        call. = FALSE
      )
    }
    reported <- R[i, j]
    if (identical(hold, "points")) {
      iv <- .implied_interval_points(R, i, j)
    } else {
      bnd <- .box_bounds(Rtest, delta)
      fb <- .cell_free_box(bnd, i, j)
      fo <- .box_feasible(fb$lo, fb$hi, fb$off, verify = verify)
      if (identical(fo$status, "feasible")) {
        xf <- fo$X[i, j]
        lo_ij <- if (
          identical(.feasible_status_at(bnd, i, j, -1, verify), "feasible")
        ) {
          -1
        } else {
          .bisect_boundary(bnd, i, j, -1, xf, verify)
        }
        hi_ij <- if (
          identical(.feasible_status_at(bnd, i, j, 1, verify), "feasible")
        ) {
          1
        } else {
          .bisect_boundary(bnd, i, j, 1, xf, verify)
        }
        iv <- c(lo = lo_ij, hi = hi_ij)
      } else {
        iv <- c(lo = NA_real_, hi = NA_real_)
      }
    }
    status <- if (any(is.na(iv))) "empty" else "feasible"
    if (is.na(reported)) {
      edit <- NA_real_
      if (status == "feasible") status <- "feasible" # imputation
    } else {
      edit <- if (status == "feasible") {
        .signed_gap(reported - delta, reported + delta, iv["lo"], iv["hi"])
      } else {
        NA_real_
      }
    }
    rows[[r]] <- data.frame(
      i = i,
      j = j,
      reported = reported,
      lo = unname(iv["lo"]),
      hi = unname(iv["hi"]),
      required_edit = unname(edit),
      hold = hold,
      status = status,
      stringsAsFactors = FALSE
    )
  }
  tibble::as_tibble(do.call(rbind, rows))
}
