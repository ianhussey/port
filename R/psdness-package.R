#' psdness: Certify Whether a Rounded Correlation Matrix Can Be PSD
#'
#' Given a reported correlation matrix whose off-diagonal entries are rounded
#' to a fixed number of decimals, [check_corr_psd()] decides whether *any*
#' positive semidefinite (PSD) matrix is consistent with the induced rounding
#' box. The tool is designed for research-integrity / forensic metascience: its
#' payload is certifying **impossibility** (no PSD matrix fits the rounding
#' box), which is evidence that a reported matrix cannot be a genuine
#' correlation matrix.
#'
#' See [check_corr_psd()] for the main entry point and the package README for
#' the mathematics and the rigor guarantee.
#'
#' @keywords internal
"_PACKAGE"

## Null-coalescing helper (kept internal; avoids an rlang dependency).
`%||%` <- function(a, b) if (is.null(a)) b else a
