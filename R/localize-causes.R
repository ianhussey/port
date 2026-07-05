# -----------------------------------------------------------------------------
# Causes taxonomy: benign structural rank-deficiency vs substantive violation.
#
# van Tilburg & van Tilburg (2023) stress that impossibility (really rank
# deficiency) is OFTEN BENIGN and does not imply questionable research practices:
# common structural generators include a dummy per level of a categorical
# variable, including both a composite and its subscores, and splicing
# correlations from different populations (see also Lorenzo-Seva & Ferrando
# 2021 on non-positive-definite matrices in factor analysis). For a forensic
# tool the report must surface these BEFORE anything that reads as a
# data-integrity inference.
#
# Two axes, both from quantities already computed:
#   * benign vs substantive: severity_max relative to delta. A structurally
#     singular matrix rounded to d decimals sits within rounding of a valid
#     (boundary) correlation matrix, so its severity is on the rounding scale;
#     a fabricated / grossly inconsistent matrix needs edits far beyond rounding.
#   * generator: the pattern of the smallest-eigenvalue eigenvector (the near-
#     dependency direction). A same-sign combination on a subset ~ variables sum
#     to a constant (ipsative / full dummy set); one dominant entry opposite a
#     cluster ~ a composite built from its subscores; support on 2-3 variables ~
#     a localized (cellular) violation, not a structural dependency.
#
# These are conservative, clearly-labelled heuristics, never a verdict on intent.
# -----------------------------------------------------------------------------

.structural_diagnosis <- function(R, delta, severity_max,
                                  active_frac = 0.25) {
  p <- nrow(R)
  eg <- eigen(R, symmetric = TRUE)
  lam_min <- eg$values[p]
  v <- eg$vectors[, p]
  if (v[which.max(abs(v))] < 0) v <- -v      # sign convention
  amax <- max(abs(v))
  active <- which(abs(v) >= active_frac * amax)
  n_active <- length(active)

  ratio <- if (is.finite(severity_max) && delta > 0) severity_max / delta else NA_real_
  severity_class <- if (is.na(ratio)) "unknown"
                    else if (ratio > 10) "substantive"
                    else if (ratio <= 3) "near-boundary"
                    else "moderate"

  va <- v[active]
  same_sign <- all(sign(va) == sign(va[1]))
  pattern <- if (n_active <= 2L) {
    "localized"                              # a single cell / pair
  } else if (same_sign) {
    "ipsative"                               # same-sign combo ~ sums to a constant
  } else if (n_active == 3L) {
    "localized"                              # mixed-sign triple ~ a triad
  } else {
    ord <- order(abs(va), decreasing = TRUE)
    top_sign <- sign(va[ord[1]])
    rest_sign <- sign(va[ord[-1]])
    if (top_sign != 0 && mean(rest_sign == -top_sign) >= 0.6) "composite" else "diffuse"
  }

  list(lambda_min = lam_min, direction = v, active = active, n_active = n_active,
       severity_ratio = ratio, severity_class = severity_class, pattern = pattern)
}

# One-sentence, honest interpretation of the diagnosis, benign explanations first.
.structural_note <- function(diag) {
  sc <- diag$severity_class
  pat <- diag$pattern
  vars <- paste(diag$active, collapse = ", ")

  if (identical(sc, "substantive")) {
    return(paste0(
      "Substantive: no correlation matrix lies within rounding of the reported ",
      "values (smallest correction is ~", signif(diag$severity_ratio, 2),
      "x the rounding step), so it is not explained by rounding of a structural ",
      "dependency."))
  }

  benign_tail <- paste0(
    " Structural rank-deficiency of this kind is often a benign modelling ",
    "artifact, not evidence of a data problem (van Tilburg & van Tilburg 2023; ",
    "Lorenzo-Seva & Ferrando 2021).")
  scale <- if (identical(sc, "near-boundary"))
    "The matrix sits within rounding of a singular (boundary) correlation matrix. "
  else
    "The matrix is modestly beyond the rounding boundary. "

  gen <- switch(pat,
    ipsative = paste0(
      "The near-dependency is a same-sign combination of variables {", vars,
      "}, consistent with an ipsative / full set-of-dummies structure (they ~sum ",
      "to a constant)."),
    composite = paste0(
      "The near-dependency has one variable opposite a cluster of {", vars,
      "}, consistent with a composite built from its subscores."),
    localized = paste0(
      "The near-dependency is localized to variables {", vars,
      "}, i.e. a small-set (cellular) violation rather than a recognizable ",
      "structural dependency."),
    diffuse = paste0(
      "The near-dependency is spread across variables {", vars,
      "} with no recognizable structural pattern."))

  benign <- pat %in% c("ipsative", "composite")
  paste0(scale, gen, if (benign) benign_tail else "")
}
