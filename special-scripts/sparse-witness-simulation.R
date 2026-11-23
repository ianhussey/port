# Simulation behind sparse-witness.md.
#
# Question: on rounded correlation tables that port's original witness search
# leaves undecided, which extra families of test directions (or constructions)
# turn the verdict into a certificate?
#
# Generator: a low-rank-plus-small-noise correlation matrix (near the PSD
# boundary), one cell perturbed by up to +-0.08, rounded to 2 decimals
# (delta = 0.005). "Undecided" = the original witness search finds no
# certificate (B_upper >= 0) and neither the constructions nor POCS certify
# consistency. The original search is .witness_search(sparse = FALSE).
#
# Part 1 compares the families on these tables. The searches were tuned on
# them, so Part 2 asks how much of the result carries over to fresh tables and
# larger p. Part 3 backs Section 4.2 of forensic-value.qmd: how large a change
# to one cell of a realistic table must be before it is proven impossible, and
# how often a random error makes a table impossible at all.
#
# Run from the package root:  Rscript special-scripts/sparse-witness-simulation.R
# Takes about 10 minutes.

suppressMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(tidyr)
})

delta <- 0.005

# --- Generator ----------------------------------------------------------------

# One table from the generator, or NULL if the perturbed cell left [-1, 1].
# `unique_range` bounds the unique variances. Small ones make the variables
# nearly redundant, which puts the table near the PSD boundary.
draw <- function(p, unique_range = c(0.001, 0.05)) {
  rank <- sample(2:max(2, p - 2), 1)
  loadings <- matrix(rnorm(p * rank), p)
  unique_variances <- runif(p, unique_range[1], unique_range[2])
  population <- cov2cor(tcrossprod(loadings) + diag(unique_variances))
  cell <- sample(p, 2)
  population[cell[1], cell[2]] <- population[cell[2], cell[1]] <-
    population[cell[1], cell[2]] + runif(1, -0.08, 0.08)
  if (any(abs(population[upper.tri(population)]) > 1)) {
    return(NULL)
  }
  round(population, 2)
}

# n draws from the generator, without the ones that left [-1, 1].
draw_tables <- function(p, n, unique_range = c(0.001, 0.05)) {
  map(seq_len(n), \(i) draw(p, unique_range)) |>
    compact()
}

# --- Candidate families ---------------------------------------------------------

proves <- function(cor_table, direction) {
  !is.null(direction) &&
    .witness_bound(cor_table, direction, delta)$B_upper < 0
}

# Variable subsets as index vectors: all bipartitions (each counted once), and
# all subsets of size 3 to `max_size`.
bipartitions <- function(p) {
  map(seq_len(2^(p - 1) - 1), \(mask) which(bitwAnd(mask, 2^(0:(p - 1))) > 0))
}

subsets <- function(p, max_size = p) {
  map(3:max_size, \(size) utils::combn(p, size, simplify = FALSE)) |>
    list_flatten()
}

# Best single test direction, to numerical precision. Within one cell (a
# support and a sign pattern) the box bound is the quadratic form of the vertex
# matrix R + delta * signs signs' (unit diagonal) on the support, so a
# minimiser is an eigenvector of it whose signs match. Enumerating all
# (3^p - 1) / 2 cells therefore finds the global minimum; the result still goes
# through the rigorous bound like any direction.
best_direction <- function(cor_table) {
  p <- nrow(cor_table)
  best_value <- Inf
  best <- NULL
  for (support_mask in seq_len(2^p - 1)) {
    support <- which(bitwAnd(support_mask, 2^(0:(p - 1))) > 0)
    size <- length(support)
    if (size < 2) {
      next
    }
    for (sign_mask in 0:(2^(size - 1) - 1)) {
      flips <- bitwAnd(sign_mask, 2^(seq_len(size - 1) - 1)) > 0
      signs <- c(1, ifelse(flips, -1, 1))
      vertex <- cor_table[support, support] + delta * outer(signs, signs)
      diag(vertex) <- 1
      decomposition <- eigen(vertex, symmetric = TRUE)
      for (j in size:1) {
        if (decomposition$values[j] >= best_value) {
          break
        }
        direction <- decomposition$vectors[, j]
        if (direction[1] < 0) {
          direction <- -direction
        }
        if (all(direction * signs > -1e-12)) {
          best_value <- decomposition$values[j]
          best <- replace(numeric(p), support, direction)
          break
        }
      }
    }
  }
  best
}

inverse_sqrt <- function(matrix) {
  decomposition <- eigen(matrix, symmetric = TRUE)
  decomposition$vectors %*%
    (1 / sqrt(decomposition$values) * t(decomposition$vectors))
}

# Block (canonical-correlation) witness for the partition `part` | rest
# (sparse-witness.md, Section 5.1); NULL unless both diagonal blocks are
# positive definite.
block_direction <- function(cor_table, part) {
  rest <- setdiff(seq_len(nrow(cor_table)), part)
  within_part <- cor_table[part, part, drop = FALSE]
  within_rest <- cor_table[rest, rest, drop = FALSE]
  is_positive_definite <- function(matrix) {
    min(eigen(matrix, symmetric = TRUE, only.values = TRUE)$values) > 1e-10
  }
  if (
    !is_positive_definite(within_part) || !is_positive_definite(within_rest)
  ) {
    return(NULL)
  }
  scale_part <- inverse_sqrt(within_part)
  scale_rest <- inverse_sqrt(within_rest)
  canonical <- svd(
    scale_part %*% cor_table[part, rest, drop = FALSE] %*% scale_rest
  )
  direction <- numeric(nrow(cor_table))
  direction[part] <- scale_part %*% canonical$u[, 1]
  direction[rest] <- -scale_rest %*% canonical$v[, 1]
  direction
}

# Does some multiple s in the regression-residual direction (beta_j; -s) of
# variable j certify the table (sparse-witness.md, Section 5.1)?
residual_s_proves <- function(cor_table, j) {
  residual <- .rsquared_direction(cor_table, j)
  if (is.null(residual)) {
    return(FALSE)
  }
  beta <- residual$v / -residual$v[j]
  bound_at <- function(s) {
    direction <- beta
    direction[j] <- -s
    .witness_bound(cor_table, direction, delta)$B_upper
  }
  stats::optimize(bound_at, c(0.01, 3))$objective < 0
}

# Cross-block shrinkage (sparse-witness.md, Section 5.3): every cell between
# `part` and the rest moved towards zero by up to delta.
cross_shrink <- function(cor_table, part) {
  rest <- setdiff(seq_len(nrow(cor_table)), part)
  shrunk <- cor_table
  shrunk[part, rest] <- sign(cor_table[part, rest]) *
    pmax(0, abs(cor_table[part, rest]) - delta)
  shrunk[rest, part] <- t(shrunk[part, rest])
  shrunk
}

# --- Part 1: families on the tables the original search leaves undecided -------

# Which families certify an undecided table. `shrink` counts consistency
# certificates; all others count inconsistency certificates. The greedy
# elimination is run without polishing, as implemented (each round's winner
# polished on its support), and with every trial polished; the margin-free POCS
# direction (.pocs_witness(), the last resort of check_corr_psd()) with a cap of
# 1000 iterations and as implemented.
try_families <- function(cor_table) {
  p <- nrow(cor_table)
  bound <- function(direction) .witness_bound(cor_table, direction, delta)
  polish <- function(direction) {
    .polish_witness(cor_table, delta, direction, on_support = TRUE)
  }
  box <- .box_bounds(cor_table, delta)
  greedy_seconds <- system.time(
    greedy <- .greedy_subset_witness(cor_table, bound, polish)
  )[["elapsed"]]
  greedy_all_seconds <- system.time(
    greedy_all <- .greedy_subset_witness(cor_table, polish, polish)
  )[["elapsed"]]
  greedy_proves <- isTRUE(greedy$B_upper < 0)
  pocs_proves <- .pocs_witness(box$lo, box$hi, box$off, bound)$B_upper < 0

  tibble(
    residual_s = any(map_lgl(seq_len(p), \(j) residual_s_proves(cor_table, j))),
    block = any(map_lgl(bipartitions(p), \(part) {
      proves(cor_table, block_direction(cor_table, part))
    })),
    shrink = any(map_lgl(bipartitions(p), \(part) {
      .verify_psd(cross_shrink(cor_table, part))
    })),
    subsets_3_4 = any(map_lgl(subsets(p, 4), \(support) {
      proves(cor_table, .subset_bottom_eigvec(cor_table, support))
    })),
    subsets_all = any(map_lgl(subsets(p), \(support) {
      proves(cor_table, .subset_bottom_eigvec(cor_table, support))
    })),
    greedy_unpolished = isTRUE(
      .greedy_subset_witness(cor_table, bound, bound)$B_upper < 0
    ),
    greedy = greedy_proves,
    greedy_all = isTRUE(greedy_all$B_upper < 0),
    pocs_1000 = .pocs_witness(box$lo, box$hi, box$off, bound, 1000L)$B_upper <
      0,
    pocs = pocs_proves,
    pocs_only = pocs_proves && !greedy_proves, # where greedy elimination fails
    best_single = proves(cor_table, best_direction(cor_table)),
    check_corr_psd = check_corr_psd(cor_table, delta = delta)$verdict ==
      "inconsistent",
    greedy_seconds = greedy_seconds,
    greedy_all_seconds = greedy_all_seconds
  )
}

# Verdict of the original search for one table and, if undecided, which
# families certify it. For tables the original search proves inconsistent,
# `lost_without_triples` says whether the certificate depends on the triples
# family once greedy elimination is added.
classify <- function(cor_table) {
  if (.witness_search(cor_table, delta, sparse = FALSE)$B_upper < 0) {
    without_triples <- .witness_search(cor_table, delta, triples_max_p = 0L)
    return(tibble(
      verdict = "inconsistent",
      lost_without_triples = without_triples$B_upper >= 0
    ))
  }
  if (!is.null(.construct_consistent(cor_table, delta))) {
    return(tibble(verdict = "consistent (construction)"))
  }
  if (!is.null(.pocs_consistent(cor_table, delta))) {
    return(tibble(verdict = "consistent (POCS)"))
  }
  tibble(verdict = "undecided", table = list(cor_table)) |>
    bind_cols(try_families(cor_table))
}

set.seed(2026)
for (p in c(5, 6, 8, 10)) {
  results <- draw_tables(p, n = if (p == 10) 1500 else 2500) |>
    map(classify) |>
    list_rbind()
  undecided <- filter(results, verdict == "undecided")

  cat("\n=== p =", p, "===\n")
  print(count(results, verdict))
  cat("Undecided cases certified by each family (shrink: consistency):\n")
  print(
    summarise(undecided, across(residual_s:check_corr_psd, sum)),
    width = Inf
  )
  cat("Seconds in greedy elimination, winner / every trial polished:\n")
  print(summarise(undecided, across(ends_with("_seconds"), sum)))
  cat(
    "Original certificates lost if triples are dropped (greedy kept):",
    sum(results$lost_without_triples, na.rm = TRUE),
    "\n"
  )
  if (p == 6) {
    cat("Example table undecided before, inconsistent now:\n")
    undecided |>
      filter(check_corr_psd) |>
      pull(table) |>
      pluck(1) |>
      dput()
  }
}

# --- Part 2: held-out tables ---------------------------------------------------

# Outcome for one fresh table, or NA if it is not open: the original search
# proves it inconsistent, or it is consistent. A presumed consistent (margin
# above tau) skips POCS in check_corr_psd(), so POCS is run here as in Part 1,
# so that "open" means the same thing in both parts. For POCS certificates,
# "over 1000" marks those a cap of 1000 iterations would have missed.
holdout_outcome <- function(cor_table) {
  if (.witness_search(cor_table, delta, sparse = FALSE)$B_upper < 0) {
    return(NA_character_)
  }
  result <- check_corr_psd(cor_table, delta = delta)
  if (
    result$verdict == "consistent" &&
      (result$certified || !is.null(.pocs_consistent(cor_table, delta)))
  ) {
    return(NA_character_)
  }
  if (result$verdict == "consistent") {
    return("presumed consistent")
  }
  if (result$verdict == "undecided") {
    return("undecided")
  }
  if (result$tier == "witness") {
    return("sparse search")
  }
  box <- .box_bounds(cor_table, delta)
  bound <- function(direction) .witness_bound(cor_table, direction, delta)
  pocs_1000 <- .pocs_witness(box$lo, box$hi, box$off, bound, 1000L)
  if (pocs_1000$B_upper < 0) {
    "POCS, within 1000 iterations"
  } else {
    "POCS, over 1000 iterations"
  }
}

set.seed(2027)
holdout <- tibble(p = c(8, 10, 15, 20)) |>
  mutate(
    outcome = map(p, \(p) map_chr(draw_tables(p, n = 2000), holdout_outcome))
  ) |>
  unnest_longer(outcome) |>
  filter(!is.na(outcome)) |>
  count(p, outcome)

cat(
  "\n=== Held-out tables (open cases by outcome; empty outcomes omitted) ===\n"
)
holdout |>
  pivot_wider(names_from = p, names_prefix = "p = ", values_from = n) |>
  mutate(across(starts_with("p = "), \(n) replace_na(n, 0L))) |>
  print()
cat("Open cases per p:\n")
print(count(holdout, p, wt = n, name = "open"))

# --- Part 3: realistic tables (forensic-value.qmd, Section 4.2) ----------------

# A table in which some variables are tied by an exact linear dependency, plus
# `n_other` variables that share a common factor with them.
dependency_table <- function(kind, n_other, seed) {
  set.seed(seed)
  n <- 300
  common <- rnorm(n)
  core <- if (kind == "total") {
    subscales <- 0.7 * common + matrix(rnorm(n * 4), nrow = n)
    colnames(subscales) <- paste0("sub", 1:4)
    cbind(subscales, total = rowSums(subscales))
  } else {
    pre <- common + rnorm(n)
    post <- 0.6 * pre + rnorm(n)
    cbind(pre, post, change = post - pre)
  }
  others <- 0.3 * common + matrix(rnorm(n * n_other), nrow = n)
  colnames(others) <- paste0("x", seq_len(n_other))
  round(cor(cbind(core, others)), 2)
}

# Smallest change to one cell, in either direction and in steps of 0.005, that
# the previous directions and check_corr_psd() prove impossible (NA: none up to
# 0.4).
detection_threshold <- function(
  cor_table,
  cell,
  steps = seq(0.005, 0.4, by = 0.005)
) {
  found <- c(previous = NA, new = NA)
  for (step in steps) {
    for (value in cor_table[cell[1], cell[2]] + c(-step, step)) {
      if (abs(value) > 1) {
        next
      }
      changed <- cor_table
      changed[cell[1], cell[2]] <- changed[cell[2], cell[1]] <- value
      previous_proves <-
        .witness_search(changed, delta, sparse = FALSE)$B_upper < 0
      if (is.na(found[["previous"]]) && previous_proves) {
        found[["previous"]] <- step
      }
      new_proves <- check_corr_psd(changed, decimals = 2)$verdict ==
        "inconsistent"
      if (is.na(found[["new"]]) && new_proves) {
        found[["new"]] <- step
      }
    }
    if (!anyNA(found)) {
      break
    }
  }
  found
}

# fmt: skip
scenarios <- tribble(
  ~kind,    ~cell,
  "total",  c("sub1", "sub2"),
  "total",  c("sub1", "x1"),
  "change", c("pre", "post"),
  "change", c("pre", "x1")
)

thresholds <- scenarios |>
  expand_grid(n_other = c(3, 15)) |>
  mutate(
    found = pmap(list(kind, n_other, cell), \(kind, n_other, cell) {
      detection_threshold(dependency_table(kind, n_other, seed = 2), cell)
    }),
    cell = map_chr(cell, \(cell) paste(cell, collapse = " with "))
  ) |>
  unnest_wider(found)

cat("\n=== Detection thresholds (forensic-value.qmd, Section 4.2) ===\n")
print(thresholds)

# Verdicts on random tables with nearly redundant variables (the generator
# above) and with realistic unique variances.
random_verdict <- function(cor_table) {
  result <- check_corr_psd(cor_table, decimals = 2)
  if (result$verdict == "consistent" && !result$certified) {
    "presumed consistent"
  } else {
    result$verdict
  }
}

set.seed(2028)
verdicts <- expand_grid(
  p = c(8, 20),
  unique_range = list(c(0.001, 0.05), c(0.3, 0.7))
) |>
  mutate(
    verdict = map2(p, unique_range, \(p, unique_range) {
      map_chr(draw_tables(p, n = 200, unique_range), random_verdict)
    }),
    unique_range = map_chr(unique_range, toString)
  ) |>
  unnest_longer(verdict) |>
  count(p, unique_range, verdict)

cat("\n=== Verdicts on random tables (forensic-value.qmd, Section 4.2) ===\n")
print(verdicts, n = Inf)
