# Monte-Carlo replication counts: full locally and in CI (NOT_CRAN = "true"),
# scaled down to a fast smoke on CRAN so the check stays well within time limits
# while still exercising each soundness/property loop.
n_reps <- function(full, cran = 40L) {
  if (identical(Sys.getenv("NOT_CRAN"), "true")) full else min(full, cran)
}
