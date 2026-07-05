# -----------------------------------------------------------------------------
# Batch / vectorized screening helper.
# -----------------------------------------------------------------------------

#' Screen a corpus of correlation matrices
#'
#' Run [check_corr_psd()] over a list of matrices and collect the verdicts in a
#' tidy table, for screening a corpus of reported correlation matrices. The
#' `pocs` escalation tier is expected to fire rarely; its frequency is logged
#' (as a message) so you can monitor how often the cheap rigorous tiers were
#' insufficient.
#'
#' @param mats A list of matrices (each a candidate correlation matrix). May be
#'   named; names become the `id` column (otherwise a positional index is used).
#' @param decimals,delta,tau Passed through to [check_corr_psd()].
#' @param on_error How to handle a matrix that fails validation: `"row"`
#'   (default) records a verdict of `"error"` with the message; `"stop"`
#'   re-raises the error.
#' @param quiet If `TRUE`, suppress the escalation-frequency message.
#'
#' @return A [tibble][tibble::tibble] with one row per input matrix and columns
#'   `id`, `verdict`, `tier`, `margin`, `b_upper`, `delta`, `p`, and `message`.
#'
#' @examples
#' good <- diag(3)
#' bad  <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
#' check_corr_psd_batch(list(good = good, bad = bad), decimals = 2)
#'
#' @seealso [check_corr_psd()]
#' @export
check_corr_psd_batch <- function(mats, decimals = 2, delta = NULL, tau = NULL,
                                 on_error = c("row", "stop"), quiet = FALSE) {
  on_error <- match.arg(on_error)
  if (!is.list(mats)) {
    stop("`mats` must be a list of matrices.", call. = FALSE)
  }
  n <- length(mats)
  ids <- names(mats)
  if (is.null(ids)) ids <- as.character(seq_len(n))
  ids[ids == "" | is.na(ids)] <- as.character(which(ids == "" | is.na(ids)))

  verdict <- character(n)
  tier    <- character(n)
  margin  <- rep(NA_real_, n)
  b_upper <- rep(NA_real_, n)
  del     <- rep(NA_real_, n)
  pp      <- rep(NA_integer_, n)
  msg     <- rep(NA_character_, n)

  for (i in seq_len(n)) {
    res <- tryCatch(
      check_corr_psd(mats[[i]], decimals = decimals, delta = delta, tau = tau),
      error = function(e) e)
    if (inherits(res, "error")) {
      if (on_error == "stop") stop(res)
      verdict[i] <- "error"
      tier[i]    <- NA_character_
      msg[i]     <- conditionMessage(res)
      next
    }
    verdict[i] <- res$verdict
    tier[i]    <- res$tier
    margin[i]  <- res$margin
    b_upper[i] <- res$b_upper
    del[i]     <- res$delta
    pp[i]      <- res$p
    msg[i]     <- res$note %||% NA_character_
  }

  n_pocs <- sum(tier == "pocs", na.rm = TRUE)
  if (!quiet && n > 0L) {
    message(sprintf("POCS escalation tier fired on %d of %d matrices (%.1f%%).",
                    n_pocs, n, 100 * n_pocs / n))
  }

  tibble::tibble(
    id = ids,
    verdict = verdict,
    tier = tier,
    margin = margin,
    b_upper = b_upper,
    delta = del,
    p = pp,
    message = msg
  )
}
