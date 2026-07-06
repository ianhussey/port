# -----------------------------------------------------------------------------
# S3 methods and accessors for corr_psd_check.
# -----------------------------------------------------------------------------

#' @export
print.corr_psd_check <- function(x, digits = 4, ...) {
  head <- switch(x$verdict,
    inconsistent = "INCONSISTENT  (no PSD matrix fits the rounding box)",
    consistent   = if (isTRUE(x$certified))
                   "CONSISTENT    (certified: in-box PSD matrix exhibited)"
                 else "CONSISTENT    (presumed: not shown inconsistent)",
    undecided    = "UNDECIDED     (ambiguous at this precision)",
    x$verdict)
  cat("<corr_psd_check>\n")
  cat(sprintf("  verdict : %s\n", head))
  cat(sprintf("  tier    : %s\n", x$tier))
  cat(sprintf("  p       : %d      delta : %g", x$p, x$delta))
  if (!identical(x$rounding %||% "nearest", "nearest")) {
    cat(sprintf("      rounding : %s", x$rounding))
  }
  cat("\n")
  if (!is.null(x$r_min)) {
    cat(sprintf("  corr    : min %.*g  max %.*g  mean %.*g  sd %.*g  (off-diagonal)\n",
                digits, x$r_min, digits, x$r_max, digits, x$r_mean, digits, x$r_sd))
  }

  if (identical(x$verdict, "inconsistent")) {
    if (identical(x$tier, "precheck")) {
      cat("  reason  : ", paste(x$note, collapse = " "), "\n", sep = "")
    } else {
      v <- x$witness
      cat(sprintf("  margin  : B = %.*g  (B_upper = %.*g < 0)\n",
                  digits, x$margin, digits, x$b_upper))
      cat("  witness : v =", .fmt_vec(v, digits), "\n")
      cat(sprintf("  certificate: v'Rv + delta*(||v||_1^2 - ||v||_2^2) = %.*g < 0\n",
                  digits, x$margin))
      if (length(x$note)) cat("  note    : ", paste(x$note, collapse = " "), "\n", sep = "")
    }
    cat("  caution : pairwise deletion, polychoric estimation, meta-analytic\n",
        "            assembly, or upstream disattenuation can legitimately\n",
        "            produce non-PSD reported matrices; see README.\n", sep = "")
  } else {
    if (!is.na(x$b_upper)) {
      cat(sprintf("  margin  : B_upper = %.*g  (>= 0: no inconsistency certificate)\n",
                  digits, x$b_upper))
    }
    if (!is.null(x$certified_matrix)) {
      cat("  evidence: ", paste(x$note %||% "an in-box PSD matrix was verified",
                                collapse = " "), "\n", sep = "")
    } else if (!is.null(x$note)) {
      cat("  note    : ", paste(x$note, collapse = " "), "\n", sep = "")
    }
  }
  invisible(x)
}

# Compact vector formatter for printing certificates.
.fmt_vec <- function(v, digits = 4) {
  paste0("(", paste(formatC(v, digits = digits, format = "g"), collapse = ", "), ")")
}

#' Extract the inconsistency certificate from a `corr_psd_check`
#'
#' @param x A `corr_psd_check` object.
#' @return A list with `witness` (the certificate vector `v`, or `NULL` when the
#'   verdict is not `inconsistent` via the witness tier), `margin` (the value
#'   `B = v'Rv + delta*(||v||_1^2 - ||v||_2^2)`), `b_upper` (the rigorous upper
#'   bound `B_upper`), `verdict` and `tier`.
#' @examples
#' R <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
#' certificate(check_corr_psd(R, decimals = 2))
#' @export
certificate <- function(x) {
  UseMethod("certificate")
}

#' @export
certificate.corr_psd_check <- function(x) {
  list(
    witness = x$witness,
    margin = x$margin,
    b_upper = x$b_upper,
    verdict = x$verdict,
    tier = x$tier
  )
}

#' @export
certificate.default = function(x) {
  stop("`certificate()` expects a <corr_psd_check> object.", call. = FALSE)
}
