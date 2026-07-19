# -----------------------------------------------------------------------------
# Rigorous a-priori floating-point error model.
#
# The inconsistency verdict must be floating-point-safe: it must not depend on
# a solver tolerance. We achieve this with standard backward/forward error
# bounds (Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed.,
# 2002, Ch. 3). All quantities below are deliberate *over*-estimates of the
# true rounding error, so adding them as slack keeps our computed bound on the
# right (conservative) side of the exact real-arithmetic value.
# -----------------------------------------------------------------------------

# Unit roundoff u = eps/2 for IEEE double precision (round-to-nearest).
.unit_roundoff <- function() .Machine$double.eps / 2

# Higham's gamma_n = n*u / (1 - n*u), valid for n*u < 1.
#
# This is the standard constant that bounds the accumulated relative error of a
# sequence of n floating-point operations. We return Inf if the guard n*u >= 1
# is violated (astronomically large n) so that any downstream comparison fails
# safe (never spuriously declaring inconsistency).
.gamma <- function(n, u = .unit_roundoff()) {
  nu <- n * u
  if (nu >= 1) {
    return(Inf)
  }
  nu / (1 - nu)
}
