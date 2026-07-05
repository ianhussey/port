# -----------------------------------------------------------------------------
# Fault-localization integration: deterministic ruleset -> localization_verdict.
# -----------------------------------------------------------------------------

.cell_key <- function(i, j) paste(min(i, j), max(i, j), sep = "-")

.new_psd_fault <- function(localization_verdict, implicated, convergence,
                           evidence, severity, delta, p, verify,
                           notes = character(0), base_check = NULL,
                           structural = NULL, plausibility = NULL) {
  structure(
    list(
      localization_verdict = localization_verdict,
      implicated = implicated,          # list(cells = matrix|NULL, variable = int|NULL)
      convergence = convergence,        # "all" / "partial" / NA
      evidence = evidence,              # A-D + R^2 raw outputs
      severity = severity,              # severity_max/frob/best_lambda_min/witness_margin
      structural = structural,          # benign-vs-substantive cause diagnosis
      plausibility = plausibility,      # per-variable R^2 gradient (possible matrices)
      delta = delta,
      p = p,
      verify = verify,
      notes = notes,
      base_check = base_check
    ),
    class = "psd_fault"
  )
}

#' Localize the fault in an impossible-given-rounding correlation matrix
#'
#' When [check_corr_psd()] reports `impossible`, this runs a box-aware
#' fault-localization layer to support semi-automated inference about WHERE the
#' non-PSDness comes from and HOW SEVERE it is. Non-PSD is a global property, so
#' a unique culprit is generally under-identified; the tool reports an honest
#' localization class and never manufactures a unique culprit when the evidence
#' only supports a set.
#'
#' All feasibility computations respect the rounding BOX, never point values:
#' `X_ii = 1`; for cells not being interrogated, `X_ij` in
#' `[R_ij - delta, R_ij + delta]` intersected with `[-1, 1]`; `X` symmetric.
#'
#' @section Evidence (all box-aware):
#' \describe{
#'   \item{A. per-cell interval / sole culprit}{for each cell, the interval it
#'     could take given the others within rounding; a nonempty interval marks a
#'     sole-culprit candidate with a `required_edit` beyond rounding.}
#'   \item{B. impossible triples}{triples whose 3x3 box-max determinant is
#'     provably negative (sound interval bound; no solver).}
#'   \item{C. leave-one-out}{variables whose removal restores feasibility.}
#'   \item{D. sparse correction / severity}{smallest set of cells whose joint
#'     correction restores feasibility, plus severity measures.}
#'   \item{R^2 localizer}{per-variable squared multiple correlation (van Tilburg
#'     & van Tilburg 2023): variables the regression-residual witness direction
#'     rigorously flags as over-determined (`R^2 > 1`), conditional on the
#'     complementary submatrix being PD. A variable-level voter, not an arbiter.}
#' }
#'
#' @section Structural cause and plausibility:
#' Every impossible result also carries a `structural` diagnosis that separates a
#' benign near-boundary rank-deficiency (composite+subscores, a full set of
#' category dummies) from a substantive over-the-boundary violation, using the
#' severity relative to the rounding step and the pattern of the near-dependency
#' direction (van Tilburg & van Tilburg 2023; Lorenzo-Seva & Ferrando 2021). For
#' a `possible` matrix the result reports a `plausibility` gradient: the
#' per-variable reported-point `R^2` and its closest approach to the 100% ceiling.
#'
#' @section Localization verdicts:
#' `cell` (a single cell, methods converge), `cell_tentative` (single cell,
#' partial corroboration), `triad` (three cells of one impossible triangle; any
#' one edit resolves it -- not separable), `variable` (one variable's removal
#' restores feasibility), `joint` (>= 2 cells needed, no single cell suffices),
#' `diffuse` (no small explanation), or `none` (matrix was possible given
#' rounding).
#'
#' @section Rigor:
#' Impossibility-type sub-claims are sound: impossible-triple flags use the
#' interval bound; "freeing a cell alone cannot restore PSD" and "removal still
#' impossible" are, when `verify = TRUE`, confirmed by the base package's
#' witness + Rump machinery, not a solver status. Feasibility-type claims exhibit
#' a Rump-verified witnessing matrix. With `verify = FALSE`, affected claims are
#' labelled search-based.
#'
#' @param x A `corr_psd_check` object, or a numeric correlation matrix.
#' @param decimals,delta Rounding precision when `x` is a matrix (as in
#'   [check_corr_psd()]).
#' @param verify If `TRUE` (default), certify impossibility sub-claims with the
#'   witness + Rump machinery rather than trusting the search.
#' @param sparse_k Maximum cardinality for the sparse-support / joint search.
#' @param tol Magnitude threshold for treating an edit / excess as non-zero.
#'
#' @return An S3 `psd_fault` object; see the print method and [fault_evidence()].
#' @examples
#' R <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
#' localize_psd_fault(R, decimals = 2)
#' @seealso [check_corr_psd()], [impossible_triples()]
#' @export
localize_psd_fault <- function(x, decimals = 2, delta = NULL,
                               verify = TRUE, sparse_k = 3L, tol = 1e-6) {
  # Resolve input to a base check. Localization currently requires a UNIFORM
  # nearest-rounding box (no asymmetric rules, per-cell precision, or NA cells).
  if (inherits(x, "corr_psd_check")) {
    chk <- x
    R <- attr(chk, "R")
    if (is.null(R)) {
      stop("This `corr_psd_check` does not carry its matrix; pass the matrix ",
           "to localize_psd_fault() instead.", call. = FALSE)
    }
    if (!isTRUE(attr(chk, "uniform_box") %||% TRUE)) {
      stop("localize_psd_fault() supports only uniform nearest-rounding boxes ",
           "(no asymmetric rounding, per-cell precision, or NA cells yet).",
           call. = FALSE)
    }
    delta <- chk$delta
  } else {
    R <- .validate_corr(x)
    if (is.null(delta)) delta <- 0.5 * 10^(-decimals)
    chk <- check_corr_psd(R, delta = delta)
  }
  p <- nrow(R)

  # Rule 7: not impossible -> localization not applicable, but still report the
  # plausibility gradient (per-variable R^2 distance to the 100% ceiling, V4).
  if (!identical(chk$verdict, "impossible")) {
    plaus <- .plausibility_gradient(R, delta)
    note <- sprintf("Matrix is '%s' given rounding; localization not applicable.",
                    chk$verdict)
    if (isTRUE(plaus$near_ceiling)) {
      tail <- if (plaus$max_r2 >= 1)
        "exceeds 100% at the exact reported values, but possible within rounding."
      else "possible, but close to the 100% ceiling."
      note <- c(note, sprintf("Plausibility: variable %d has reported-point R^2 = %.1f%% -- %s",
        plaus$max_var, 100 * plaus$max_r2, tail))
    }
    return(.new_psd_fault(
      localization_verdict = "none",
      implicated = list(cells = NULL, variable = NULL),
      convergence = NA_character_,
      evidence = list(sole_culprit_cells = NULL, impossible_triples = NULL,
                      lofo_restoring = NULL, sparse_support = NULL, rsquared = plaus$rsquared),
      severity = list(severity_max = NA_real_, severity_frob = NA_real_,
                      best_lambda_min = NA_real_, witness_margin = chk$margin,
                      excusable_delta = NA_real_),
      delta = delta, p = p, verify = verify, plausibility = plaus,
      notes = note, base_check = chk))
  }

  # ---- Evidence A-D -------------------------------------------------------
  triples <- .scan_impossible_triples(R, delta)             # B
  lofo <- .lofo_restoring(R, delta)                         # C
  cell_int <- .cell_intervals(R, delta, verify = verify)    # A

  sole_records <- Filter(function(c) isTRUE(c$sole), cell_int)
  if (length(sole_records) > 0L) {
    sole_cells <- do.call(rbind, lapply(sole_records, function(c) c(c$i, c$j)))
    sole_edits <- vapply(sole_records, function(c) c$required_edit, numeric(1))
  } else {
    sole_cells <- matrix(integer(0), ncol = 2)
    sole_edits <- numeric(0)
  }

  sparse <- .sparse_support(R, delta, sole_cells, sole_edits,           # D1
                            verify = verify, sparse_k = sparse_k)
  sev_max <- .severity_max(R, delta, verify = verify)
  sev <- list(                                                          # severity
    severity_max = sev_max,
    severity_frob = .severity_frob(R, delta),
    best_lambda_min = .best_lambda_min(R, delta),
    witness_margin = chk$margin,
    # smallest uniform reporting half-width that could excuse the matrix:
    # impossible for every precision finer than this (see ?excusable_delta)
    excusable_delta = delta + sev_max)

  rsq <- .rsquared_evidence(R, delta)                                   # R^2 (V1+V2)
  S_r2 <- which(.rsquared_blamed(rsq))                                  # cleanly-blamed vars
  structural <- .structural_diagnosis(R, delta, sev$severity_max)       # causes (V5)

  # ---- Precompute set primitives -----------------------------------------
  S_sole_keys <- if (nrow(sole_cells) > 0L) {
    apply(sole_cells, 1L, function(rw) .cell_key(rw[1], rw[2]))
  } else character(0)
  S_sparse_keys <- if (!is.null(sparse)) {
    apply(sparse$cells, 1L, function(rw) .cell_key(rw[1], rw[2]))
  } else character(0)

  triple_keysets <- lapply(triples, function(t)
    apply(t$cells, 1L, function(rw) .cell_key(rw[1], rw[2])))
  triples_union <- unique(unlist(triple_keysets))
  common_triple_cells <- if (length(triple_keysets) > 0L) {
    Reduce(intersect, triple_keysets)
  } else character(0)

  evidence <- list(
    sole_culprit_cells = sole_records,
    impossible_triples = triples,
    lofo_restoring = lofo,
    sparse_support = sparse,
    rsquared = rsq)

  key_to_cell <- function(key) as.integer(strsplit(key, "-", fixed = TRUE)[[1]])

  # ---- Deterministic ruleset (first match) --------------------------------
  # The verdict CLASS is decided by A-D (below); the R^2 localizer and structural
  # diagnosis are voters/annotations layered on afterwards, never arbiters.
  decide <- function() {
  # Rule 1 / 2: single sole culprit.
  if (length(S_sole_keys) == 1L) {
    c_key <- S_sole_keys[1]
    in_sparse <- c_key %in% S_sparse_keys
    in_triple <- c_key %in% triples_union
    cell_rc <- matrix(key_to_cell(c_key), nrow = 1)
    if (in_sparse && in_triple) {
      return(.new_psd_fault("cell",
        implicated = list(cells = cell_rc, variable = NULL),
        convergence = "all", evidence = evidence, severity = sev,
        delta = delta, p = p, verify = verify,
        notes = "Single sole culprit corroborated by the sparse correction and an impossible triple.",
        base_check = chk))
    }
    agreed <- c(if (in_sparse) "sparse correction" else NULL,
                if (in_triple) "impossible triple" else NULL)
    disagreed <- c(if (!in_sparse) "sparse correction" else NULL,
                   if (!in_triple) "impossible triple" else NULL)
    return(.new_psd_fault("cell_tentative",
      implicated = list(cells = cell_rc, variable = NULL),
      convergence = "partial", evidence = evidence, severity = sev,
      delta = delta, p = p, verify = verify,
      notes = c(sprintf("Single sole culprit; corroborated by: %s.",
                        if (length(agreed)) paste(agreed, collapse = ", ") else "none"),
                sprintf("Not corroborated by: %s.",
                        if (length(disagreed)) paste(disagreed, collapse = ", ") else "none")),
      base_check = chk))
  }

  # Rule 3: triad -- S_sole is exactly the 3 cells of a single impossible triple.
  if (length(S_sole_keys) == 3L) {
    match_triple <- any(vapply(triple_keysets, function(ks)
      setequal(ks, S_sole_keys), logical(1)))
    if (match_triple) {
      cells_rc <- do.call(rbind, lapply(S_sole_keys, key_to_cell))
      return(.new_psd_fault("triad",
        implicated = list(cells = cells_rc, variable = NULL),
        convergence = "all", evidence = evidence, severity = sev,
        delta = delta, p = p, verify = verify,
        notes = "Not separable: any one of the three edits would resolve the violation.",
        base_check = chk))
    }
  }

  # Rule 4: variable -- exactly one variable's removal restores feasibility.
  if (length(lofo$restoring) == 1L) {
    k <- lofo$restoring[1]
    # cells concentrated in k's column, if any
    in_col <- function(keys) Filter(function(key) k %in% key_to_cell(key), keys)
    col_cells <- unique(c(in_col(S_sole_keys), in_col(S_sparse_keys)))
    likely <- if (length(col_cells) > 0L) do.call(rbind, lapply(col_cells, key_to_cell)) else NULL
    r2_agrees <- k %in% S_r2
    notes <- sprintf("Removing variable %d restores possibility given rounding.", k)
    if (r2_agrees) notes <- c(notes,
      sprintf("Corroborated by the R^2 localizer (variable %d is over-determined).", k))
    return(.new_psd_fault("variable",
      implicated = list(cells = likely, variable = k),
      convergence = if (r2_agrees) "all" else if (!is.null(likely)) "partial" else NA_character_,
      evidence = evidence, severity = sev, delta = delta, p = p, verify = verify,
      notes = notes, base_check = chk))
  }

  # Rule 5: joint -- no single cell suffices but a small set does.
  if (length(S_sole_keys) == 0L && !is.null(sparse) &&
      sparse$cardinality <= sparse_k && sparse$cardinality >= 2L) {
    return(.new_psd_fault("joint",
      implicated = list(cells = sparse$cells, variable = NULL),
      convergence = "partial", evidence = evidence, severity = sev,
      delta = delta, p = p, verify = verify,
      notes = "Requires simultaneous correction of >= 2 cells; no single cell suffices.",
      base_check = chk))
  }

  # Rule 6: diffuse.
  .new_psd_fault("diffuse",
    implicated = list(cells = NULL, variable = NULL),
    convergence = NA_character_, evidence = evidence, severity = sev,
    delta = delta, p = p, verify = verify,
    notes = "No small explanation: ranked evidence only (see fault_evidence()).",
    base_check = chk)
  }

  res <- decide()

  # Layer the causes diagnosis (V5) and the R^2 attribution (V1+V2, V6) onto the
  # verdict as annotations. These never change the verdict CLASS.
  res$structural <- structural
  r2_note <- if (length(S_r2) > 0L) {
    sprintf(paste0("R^2 localizer: variable(s) {%s} are over-determined by the ",
      "others (R^2 > 100%%) -- conditional on the other cells being correct."),
      paste(S_r2, collapse = ", "))
  } else NULL
  res$notes <- c(res$notes, .structural_note(structural), r2_note)
  res
}
