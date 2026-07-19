# -----------------------------------------------------------------------------
# Batch fault-localization over a corpus.
# -----------------------------------------------------------------------------

.implicated_str <- function(lf) {
  cells <- lf$implicated$cells
  var <- lf$implicated$variable
  parts <- character(0)
  if (!is.null(var)) {
    parts <- c(parts, sprintf("var %d", var))
  }
  if (!is.null(cells) && nrow(cells) > 0L) {
    parts <- c(
      parts,
      paste(
        apply(cells, 1L, function(rw) {
          sprintf("(%d,%d)", rw[1], rw[2])
        }),
        collapse = " "
      )
    )
  }
  if (length(parts) == 0L) NA_character_ else paste(parts, collapse = "; ")
}

#' Batch fault localization over a list of matrices
#'
#' Run [localize_psd_fault()] over a corpus and collect the localization
#' verdicts in a tidy table. A per-class summary of how often each verdict
#' occurs is logged as a message.
#'
#' @param mats A (optionally named) list of correlation matrices.
#' @param decimals,delta,verify,sparse_k Passed to [localize_psd_fault()].
#' @param quiet If `TRUE`, suppress the verdict-class summary message.
#'
#' @return A [tibble][tibble::tibble] with columns `id`, `localization_verdict`,
#'   `implicated`, `severity_max`, `severity_frob`, and `message`.
#' @examples
#' triad <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
#' localize_psd_fault_batch(list(triad = triad, ok = diag(3)))
#' @seealso [localize_psd_fault()]
#' @export
localize_psd_fault_batch <- function(
  mats,
  decimals = 2,
  delta = NULL,
  verify = TRUE,
  sparse_k = 3L,
  quiet = FALSE
) {
  if (!is.list(mats)) {
    stop("`mats` must be a list of matrices.", call. = FALSE)
  }
  n <- length(mats)
  ids <- names(mats)
  if (is.null(ids)) {
    ids <- as.character(seq_len(n))
  }
  ids[ids == "" | is.na(ids)] <- as.character(which(ids == "" | is.na(ids)))

  verdict <- character(n)
  implicated <- character(n)
  smax <- rep(NA_real_, n)
  sfrob <- rep(NA_real_, n)
  msg <- rep(NA_character_, n)

  for (i in seq_len(n)) {
    lf <- tryCatch(
      localize_psd_fault(
        mats[[i]],
        decimals = decimals,
        delta = delta,
        verify = verify,
        sparse_k = sparse_k
      ),
      error = function(e) e
    )
    if (inherits(lf, "error")) {
      verdict[i] <- "error"
      msg[i] <- conditionMessage(lf)
      next
    }
    verdict[i] <- lf$localization_verdict
    implicated[i] <- .implicated_str(lf)
    smax[i] <- lf$severity$severity_max
    sfrob[i] <- lf$severity$severity_frob
    msg[i] <- if (length(lf$notes)) lf$notes[1] else NA_character_
  }

  if (!quiet && n > 0L) {
    tab <- table(verdict)
    message(
      "Localization verdict classes: ",
      paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = ", ")
    )
  }

  tibble::tibble(
    id = ids,
    localization_verdict = verdict,
    implicated = implicated,
    severity_max = smax,
    severity_frob = sfrob,
    message = msg
  )
}
