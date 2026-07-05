# -----------------------------------------------------------------------------
# S3 methods and accessors for psd_fault.
# -----------------------------------------------------------------------------

.cells_str <- function(cells) {
  if (is.null(cells) || nrow(cells) == 0L) return("none")
  paste(apply(cells, 1L, function(rw) sprintf("(%d,%d)", rw[1], rw[2])),
        collapse = ", ")
}

#' @export
print.psd_fault <- function(x, digits = 3, ...) {
  cat("<psd_fault>\n")
  v <- x$localization_verdict
  if (identical(v, "none")) {
    cat("  localization: NONE -", x$notes[1], "\n")
    return(invisible(x))
  }

  d <- x$delta
  cat(sprintf("  Impossible given rounding (delta = %g).\n", d))

  # Headline attribution.
  sole <- Filter(function(c) isTRUE(c$sole), x$evidence$sole_culprit_cells)
  find_rec <- function(i, j) {
    for (c in sole) if (c$i == i && c$j == j) return(c)
    NULL
  }
  hl <- switch(v,
    cell = "Attributable to a single cell",
    cell_tentative = "Tentatively attributable to a single cell",
    triad = "Attributable to an inconsistent triangle (three cells)",
    variable = "Attributable to a single variable",
    joint = "Attributable to a joint set of cells (no single cell suffices)",
    diffuse = "Diffuse: no small explanation",
    v)
  cat("  Verdict:", toupper(v), "-", hl, "\n")

  if (v %in% c("cell", "cell_tentative")) {
    rc <- x$implicated$cells[1, ]
    rec <- find_rec(rc[1], rc[2])
    if (!is.null(rec)) {
      cat(sprintf(paste0("    cell (%d,%d): reported %.*g, but given the other cells ",
                         "(within rounding) it must lie in [%.*g, %.*g];\n"),
                  rec$i, rec$j, digits, rec$reported,
                  digits, rec$lo, digits, rec$hi))
      cat(sprintf("    required edit %+.*g, %s the %g rounding tolerance.\n",
                  digits, rec$required_edit,
                  if (abs(rec$required_edit) > d) "far exceeding" else "within",
                  d))
    }
    for (n in x$notes) cat("    -", n, "\n")
  } else if (identical(v, "triad")) {
    cat("    cells:", .cells_str(x$implicated$cells), "\n")
    cat("    ", x$notes[1], "\n", sep = "")
  } else if (identical(v, "variable")) {
    cat(sprintf("    variable %d (its removal restores possibility given rounding).\n",
                x$implicated$variable))
    if (!is.null(x$implicated$cells)) {
      cat("    most likely cell(s) within its column:",
          .cells_str(x$implicated$cells), "\n")
    }
  } else if (identical(v, "joint")) {
    cat("    cells:", .cells_str(x$implicated$cells), "\n")
    cat("    ", x$notes[1], "\n", sep = "")
  } else if (identical(v, "diffuse")) {
    nt <- length(x$evidence$impossible_triples)
    cat(sprintf("    %d impossible triple(s); leave-one-out restoring vars: {%s}.\n",
                nt, paste(x$evidence$lofo_restoring$restoring, collapse = ", ")))
    cat("    See fault_evidence() for the full ranked evidence.\n")
  }

  # Severity (always).
  s <- x$severity
  cat(sprintf("  Severity: largest single edit needed %.*g", digits, s$severity_max))
  if (is.finite(s$severity_max)) {
    cat(sprintf(" (%s; delta = %g)",
                if (s$severity_max > 10 * d) "egregious, rounding cannot excuse it"
                else if (s$severity_max <= 2 * d) "borderline / precision-limited"
                else "beyond rounding", d))
  }
  cat(sprintf(";\n            total mass (Frobenius) %.*g; best achievable lambda_min %.*g.\n",
              digits, s$severity_frob, digits, s$best_lambda_min))
  if (!isTRUE(x$verify)) cat("  (verify = FALSE: impossibility sub-claims are search-based.)\n")
  invisible(x)
}

#' Accessor for the raw fault-localization evidence
#'
#' @param x A `psd_fault` object.
#' @return The `evidence` list: `sole_culprit_cells` (component A),
#'   `impossible_triples` (B), `lofo_restoring` (C), and `sparse_support` (D).
#' @examples
#' R <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
#' fault_evidence(localize_psd_fault(R, decimals = 2))
#' @export
fault_evidence <- function(x) {
  if (!inherits(x, "psd_fault")) {
    stop("`fault_evidence()` expects a <psd_fault> object.", call. = FALSE)
  }
  x$evidence
}
