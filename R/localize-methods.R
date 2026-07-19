# -----------------------------------------------------------------------------
# S3 methods and accessors for psd_fault.
# -----------------------------------------------------------------------------

.cells_str <- function(cells) {
  if (is.null(cells) || nrow(cells) == 0L) {
    return("none")
  }
  paste(
    apply(cells, 1L, function(rw) sprintf("(%d,%d)", rw[1], rw[2])),
    collapse = ", "
  )
}

#' @export
print.psd_fault <- function(x, digits = 3, ...) {
  cat("<psd_fault>\n")
  v <- x$localization_verdict
  d <- x$delta

  if (identical(v, "none")) {
    cat("  localization: NONE -", x$notes[1], "\n")
    for (n in x$notes[-1]) {
      cat("  ", n, "\n", sep = "")
    }
    return(invisible(x))
  }

  cat(sprintf("  Inconsistent given rounding (delta = %g).\n", d))

  sole <- Filter(function(c) isTRUE(c$sole), x$evidence$sole_culprit_cells)
  find_rec <- function(i, j) {
    for (c in sole) {
      if (c$i == i && c$j == j) return(c)
    }
    NULL
  }
  hl <- switch(
    v,
    cell = "Attributable to a single cell",
    cell_tentative = "Tentatively attributable to a single cell",
    triad = "Attributable to an inconsistent triangle (three cells)",
    variable = "Attributable to a single variable",
    joint = "Attributable to a joint set of cells (no single cell suffices)",
    diffuse = "Diffuse: no small explanation",
    v
  )
  cat("  Verdict:", toupper(v), "-", hl, "\n")

  # Verdict-specific detail.
  if (v %in% c("cell", "cell_tentative")) {
    rc <- x$implicated$cells[1, ]
    rec <- find_rec(rc[1], rc[2])
    if (!is.null(rec)) {
      cat(sprintf(
        paste0(
          "    cell (%d,%d): reported %.*g, but given the other cells ",
          "(within rounding) it must lie in [%.*g, %.*g];\n"
        ),
        rec$i,
        rec$j,
        digits,
        rec$reported,
        digits,
        rec$lo,
        digits,
        rec$hi
      ))
      cat(sprintf(
        "    required edit %+.*g, %s the %g rounding tolerance.\n",
        digits,
        rec$required_edit,
        if (abs(rec$required_edit) > d) "far exceeding" else "within",
        d
      ))
    }
  } else if (v %in% c("triad", "joint")) {
    cat("    cells:", .cells_str(x$implicated$cells), "\n")
  } else if (identical(v, "variable")) {
    cat(sprintf(
      "    variable %d (its removal restores consistency given rounding).\n",
      x$implicated$variable
    ))
    if (!is.null(x$implicated$cells)) {
      cat(
        "    most likely cell(s) within its column:",
        .cells_str(x$implicated$cells),
        "\n"
      )
    }
  } else if (identical(v, "diffuse")) {
    nt <- length(x$evidence$inconsistent_triples)
    cat(sprintf(
      "    %d inconsistent triple(s); leave-one-out restoring vars: {%s}.\n",
      nt,
      paste(x$evidence$lofo_restoring$restoring, collapse = ", ")
    ))
  }

  # Structural cause (benign explanations first).
  if (!is.null(x$structural)) {
    cat("  Structural:", .structural_note(x$structural), "\n")
  }

  # Severity (always).
  s <- x$severity
  cat(sprintf(
    "  Severity: largest single edit needed %.*g",
    digits,
    s$severity_max
  ))
  if (is.finite(s$severity_max)) {
    cat(sprintf(
      " (%s; delta = %g)",
      if (s$severity_max > 10 * d) {
        "egregious, rounding cannot excuse it"
      } else if (s$severity_max <= 2 * d) {
        "borderline / precision-limited"
      } else {
        "beyond rounding"
      },
      d
    ))
  }
  cat(sprintf(
    ";\n            total mass (Frobenius) %.*g; best achievable lambda_min %.*g.\n",
    digits,
    s$severity_frob,
    digits,
    s$best_lambda_min
  ))

  # Excusable-imprecision headline: verdict + severity in precision units.
  w <- s$excusable_delta
  if (is.finite(w) && w > d) {
    d_ok <- .excusable_decimals(w)
    cat(sprintf(
      "  Excusable only by mis-reporting beyond +-%.*g per entry",
      digits,
      w
    ))
    if (is.finite(d_ok) && d_ok >= 1) {
      cat(sprintf(
        " -- i.e. only if the values were really rounded to %d decimal place%s or fewer.\n",
        d_ok,
        if (d_ok == 1) "" else "s"
      ))
    } else if (is.finite(d_ok) && d_ok < 1) {
      cat(" -- no conventional rounding precision could excuse it.\n")
    } else {
      cat(".\n")
    }
  }

  # Benign-generator caution (before any inference of a reporting error).
  cat(
    "  Caution: rule out legitimate generators of non-PSD reported matrices\n",
    "  first -- pairwise deletion (per-cell n), polychoric/tetrachoric\n",
    "  estimation, meta-analytically assembled cells, upstream disattenuation.\n",
    sep = ""
  )
  if (
    !is.null(x$structural) &&
      identical(x$structural$severity_class, "substantive")
  ) {
    cat(sprintf(
      "  (The required correction here is ~%sx the rounding step, which those\n  rarely produce.)\n",
      signif(x$structural$severity_ratio, 2)
    ))
  }

  # Notes (verdict corroboration + R^2 attribution).
  drop <- if (!is.null(x$structural)) {
    .structural_note(x$structural)
  } else {
    character(0)
  }
  extra <- setdiff(x$notes, drop)
  for (n in extra) {
    cat("  -", n, "\n")
  }
  if (!isTRUE(x$verify)) {
    cat("  (verify = FALSE: inconsistency sub-claims are search-based.)\n")
  }
  invisible(x)
}

#' Accessor for the raw fault-localization evidence
#'
#' @param x A `psd_fault` object.
#' @return The `evidence` list: `sole_culprit_cells` (component A),
#'   `inconsistent_triples` (B), `lofo_restoring` (C), `sparse_support` (D), and
#'   `rsquared` (per-variable R^2 localizer).
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
