# -----------------------------------------------------------------------------
# Rump-style verified positive semidefiniteness (for the POSSIBLE side).
#
# A symmetric matrix A is *provably* PSD in floating point if a standard
# Cholesky factorization of (A - c*I) completes, where c rigorously bounds the
# backward error of the Cholesky (Rump, "Verification of positive
# definiteness", BIT 46:433-452, 2006; Higham 2002, Ch. 10).
#
# Derivation of c: if Cholesky of B = A - c*I completes, producing computed
# factor Rhat, then Rhat' Rhat = B + dB with |dB| <= gamma_{n+1} |Rhat'||Rhat|
# (Higham Thm 10.3). Rhat' Rhat is exactly PSD, so
#     lambda_min(B) >= -||dB||_2 >= -gamma_{n+1} * ||Rhat||_F^2
#                   >= -gamma_{n+1} * trace(B) >= -gamma_{n+1} * trace(A).
# Since lambda_min(A) = lambda_min(B) + c, choosing c >= gamma_{n+1}*trace(A)
# gives lambda_min(A) >= 0, i.e. A is PSD. The test is one-sided and cheap: it
# may fail to certify a genuinely-but-marginally PSD matrix, which is safe (we
# then fall through to "undecided"), and it never certifies a non-PSD matrix.
# -----------------------------------------------------------------------------

# Return TRUE iff A is verified PSD by the one-sided Cholesky test.
.verify_psd <- function(A) {
  n <- nrow(A)
  u <- .unit_roundoff()
  trA <- sum(diag(A))
  if (!is.finite(trA) || trA <= 0) {
    # A PSD matrix with unit diagonal has trace n > 0; a non-positive trace
    # cannot be PSD.
    return(FALSE)
  }
  c <- .gamma(n + 1L, u) * trA
  # Inflate c to cover rounding while forming c and the PSD-vs-PD gap.
  c <- c * 1.0001 + 2 * u * trA
  B <- A - diag(c, n)
  isTRUE(tryCatch({
    chol(B)  # errors on a non-(numerically-)PD matrix
    TRUE
  }, error = function(e) FALSE))
}
