# -----------------------------------------------------------------------------
# R^2 / VIF variable-level localizer (van Tilburg & van Tilburg 2023).
#
# For variable i with correlation vector c_i to the others and their submatrix
# M_{-i}, the squared multiple correlation is R^2_i = c_i' M_{-i}^{-1} c_i, and
# the matrix is inconsistent if any R^2_i > 1 (V1). The key identity used here:
# putting variable i last and taking the regression-residual direction
#     v = (M_{-i}^{-1} c_i ; -1),
# one gets  v' M v = 1 - R^2_i. So "R^2_i > 1" is exactly a witness vector, the
# one concentrated on variable i. We therefore obtain the RIGOROUS, box-aware,
# variable-localized test for free by feeding v through the existing box witness
# bound (.witness_box_bound) -- no interval determinant, no verified inversion.
# The inversion only *seeds a direction* and is never trusted; the verdict is the
# sound box-witness sign. Clean attribution to variable i additionally requires
# the complementary submatrix M_{-i} to be PD (V6); otherwise the direction is
# still a sound witness but the "blame variable i" reading is withheld.
# -----------------------------------------------------------------------------

# Regression-residual witness direction for variable i, seeded from the reported
# point matrix. Returns the (unit) direction, R^2_i, and whether M_{-i} is
# verified PD. NULL if the direction cannot be formed.
.rsquared_direction <- function(R, i) {
  p <- nrow(R)
  others <- setdiff(seq_len(p), i)
  M <- R[others, others, drop = FALSE]
  ci <- R[others, i]
  w <- tryCatch(solve(M, ci), error = function(e)
    tryCatch(solve(M + 1e-8 * diag(nrow(M)), ci), error = function(e2) NULL))
  if (is.null(w) || any(!is.finite(w))) return(NULL)
  r2 <- sum(ci * w)                         # R^2_i = c_i' M_{-i}^{-1} c_i
  v <- numeric(p); v[i] <- -1; v[others] <- w
  nrm <- sqrt(sum(v * v))
  list(v = v / nrm, r2 = r2, complement_pd = .verify_psd(M))
}

# Per-variable R^2 evidence. box_inconsistent is the rigorous box-witness verdict
# along the residual direction; complement_pd gates clean attribution (V6).
.rsquared_evidence <- function(R, delta) {
  p <- nrow(R)
  bnd <- .box_bounds(R, delta)
  recs <- vector("list", p)
  for (i in seq_len(p)) {
    d <- .rsquared_direction(R, i)
    if (is.null(d)) {
      recs[[i]] <- list(variable = i, r2 = NA_real_, vif = NA_real_,
                        complement_pd = FALSE, box_inconsistent = NA, b_upper = NA_real_)
      next
    }
    wb <- .witness_box_bound(bnd$lo, bnd$hi, bnd$off, d$v)
    vif <- if (is.finite(d$r2) && d$r2 != 1) 1 / (1 - d$r2) else Inf
    recs[[i]] <- list(variable = i, r2 = d$r2, vif = vif,
                      complement_pd = isTRUE(d$complement_pd),
                      box_inconsistent = isTRUE(wb$B_upper < 0), b_upper = wb$B_upper)
  }
  recs
}

# Variables that the R^2 localizer cleanly blames: the residual direction is a
# rigorous box witness AND the complementary submatrix is PD (V6).
.rsquared_blamed <- function(rsq) {
  vapply(rsq, function(z) isTRUE(z$box_inconsistent) && isTRUE(z$complement_pd),
         logical(1))
}

# Plausibility gradient (V4): for a CONSISTENT matrix, per-variable reported-point
# R^2 and the closest approach to the 100% ceiling. Only variables with a PD
# complement give a meaningful R^2.
.plausibility_gradient <- function(R, delta, ceiling = 0.95) {
  rsq <- .rsquared_evidence(R, delta)
  ok <- vapply(rsq, function(z) isTRUE(z$complement_pd) && is.finite(z$r2), logical(1))
  if (any(ok)) {
    idx <- which(ok)
    r2s <- vapply(rsq[idx], function(z) z$r2, numeric(1))
    m <- which.max(r2s)
    max_var <- idx[m]; max_r2 <- r2s[m]
  } else {
    max_var <- NA_integer_; max_r2 <- NA_real_
  }
  list(rsquared = rsq, max_var = max_var, max_r2 = max_r2,
       near_ceiling = is.finite(max_r2) && max_r2 >= ceiling)
}
