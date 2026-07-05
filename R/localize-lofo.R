# -----------------------------------------------------------------------------
# Component C: leave-one-variable-out (variable-level localizer).
#
# For each variable k, drop its row/column and run the base box-inconsistency
# check on the (p-1)x(p-1) submatrix at the same delta. A variable "restores"
# feasibility if its removal makes the submatrix consistent-given-rounding.
# Reuses the base package's rigorous machinery directly: "still inconsistent" is a
# witness-certified claim; "restored" is a Rump-verified consistency.
# -----------------------------------------------------------------------------

# Returns a list with:
#   restoring : integer vector of variables whose removal restores consistency
#   detail    : per-variable verdict on the leave-one-out submatrix
.lofo_restoring <- function(R, delta) {
  p <- nrow(R)
  verdicts <- character(p)
  restoring <- integer(0)
  for (k in seq_len(p)) {
    sub <- R[-k, -k, drop = FALSE]
    res <- check_corr_psd(sub, delta = delta)
    verdicts[k] <- res$verdict
    if (identical(res$verdict, "consistent")) restoring <- c(restoring, k)
  }
  list(restoring = restoring, verdicts = verdicts)
}
