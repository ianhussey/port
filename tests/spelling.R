if (requireNamespace("spelling", quietly = TRUE)) {
  spelling::spell_check_test(vignettes = FALSE, error = TRUE,
                             skip_on_cran = TRUE)
}
