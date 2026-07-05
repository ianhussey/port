# -----------------------------------------------------------------------------
# Print method for disattenuation_check: the plain-language forensic conclusion.
# -----------------------------------------------------------------------------

.fmt_cell <- function(cells, D) {
  if (is.null(cells) || nrow(cells) == 0L) return(character(0))
  vapply(seq_len(nrow(cells)), function(k) {
    i <- cells[k, 1]; j <- cells[k, 2]
    sprintf("(%d,%d) disattenuates to %.2f", i, j, D[i, j])
  }, character(1))
}

#' @export
print.disattenuation_check <- function(x, digits = 2, ...) {
  cat("<disattenuation_check>  mode:", x$mode, "\n")
  thr <- x$thresholds
  cut <- x$max_plausible_r
  cat(sprintf("  Observed matrix: %s given rounding.\n", x$observed_verdict))

  bind_txt <- function(which_bind) {
    if (identical(which_bind, "PSD")) "loses positive semidefiniteness"
    else "produces a corrected correlation above 1"
  }
  # Box-sound boundary when available; the point closed form as reference.
  rho_box <- thr$rho_inconsistent_box
  rho_hdl <- if (is.finite(rho_box %||% NA_real_) && (rho_box %||% 0) > 0) rho_box
             else thr$rho_inconsistent
  box_txt <- if (identical(rho_hdl, rho_box)) {
    sprintf(" (rounding-aware; boundary at exact reported values: %.*f)",
            digits, thr$rho_inconsistent)
  } else ""

  # ---- critical mode (reliabilities not reported) -------------------------
  if (identical(x$mode, "critical")) {
    cat(sprintf(paste0("  Critical reliability: the disattenuated matrix is valid ",
                       "only for reliability >= %.*f%s (it %s below that).\n"),
                digits, rho_hdl, box_txt, bind_txt(thr$inconsistent_binds)))
    if (cut < 1) {
      cat(sprintf(paste0("  It is also *plausible* (all corrected |r| <= %.2f) ",
                         "only for reliability >= %.*f.\n"),
                  cut, digits, thr$rho_plausible))
    }
    fl <- x$plausible_floor
    if (!is.null(fl)) {
      if (fl < rho_hdl) {
        cat(sprintf(paste0("  Taking %.*f as the lowest reliability these measures ",
          "plausibly attain, the whole range [%.*f, %.*f) yields an INCONSISTENT ",
          "construct matrix: unless every measure was unusually reliable ",
          "(>= %.*f), the disattenuation cannot be valid.\n"),
          digits, fl, digits, fl, digits, rho_hdl, digits, rho_hdl))
      } else {
        cat(sprintf(paste0("  Even the lowest reliability these measures plausibly ",
          "attain (%.*f) clears the %.*f threshold, so no reachable reliability ",
          "makes the disattenuated matrix inconsistent.\n"),
          digits, fl, digits, rho_hdl))
      }
    }
    cat(sprintf("  (lambda_min(observed) = %.*f; max|r| = %.*f.)\n",
                digits, thr$lambda_min, digits, thr$max_r))
    return(invisible(x))
  }

  # ---- reported reliabilities ---------------------------------------------
  rel_txt <- if (x$mode == "common") sprintf("reliability = %.2f", x$reliability)
             else sprintf("per-variable reliabilities (min %.2f)", min(x$reliability))
  cat(sprintf("  Disattenuated (%s): %s.\n", rel_txt, toupper(x$verdict)))
  cat(sprintf("    Largest corrected correlation (at reported values): %.2f.\n",
              x$max_disattenuated))

  fwd <- x$forward; D <- x$disattenuated
  if (identical(x$verdict, "inconsistent")) {
    if (fwd$range_inconsistent) {
      cells <- .fmt_cell(fwd$inconsistent_cells, D)
      cat("    A corrected correlation exceeds 1: ", paste(cells, collapse = "; "),
          " -- an inconsistent correlation.\n", sep = "")
    }
    if (fwd$psd_inconsistent) {
      cat("    The disattenuated matrix is not positive semidefinite ",
          "(no valid construct correlation matrix fits, even within rounding).\n", sep = "")
    }
  } else if (identical(x$verdict, "consistent but implausible")) {
    cells <- .fmt_cell(fwd$implausible_cells, D)
    cat(sprintf("    Valid, but a corrected correlation exceeds the plausibility cutoff %.2f: %s.\n",
                cut, paste(cells, collapse = "; ")))
  }

  # thresholds + headroom (common-reliability model)
  cat(sprintf(paste0("    Valid only if reliability >= %.*f%s (it %s below that)"),
              digits, rho_hdl, box_txt, bind_txt(thr$inconsistent_binds)))
  if (cut < 1) {
    cat(sprintf("; valid AND plausible only if >= %.*f", digits, thr$rho_plausible))
  }
  cat(".\n")
  if (x$mode == "common" && is.finite(x$headroom)) {
    if (x$headroom < 0) {
      cat(sprintf("    The reported %.2f falls %.2f short of the inconsistency boundary.\n",
                  x$reliability, -x$headroom))
    } else {
      cat(sprintf("    The reported %.2f clears the inconsistency boundary by %.2f.\n",
                  x$reliability, x$headroom))
    }
  }
  if (x$mode == "per_variable") {
    cat("    (Thresholds shown are for the common-reliability model; ",
        "per-variable critical reliabilities are a separate step.)\n", sep = "")
  }
  cat(sprintf("    (lambda_min(observed) = %.*f; max|r| = %.*f.)\n",
              digits, thr$lambda_min, digits, thr$max_r))
  invisible(x)
}
