# port Shiny app

A browser front end to the [port](https://github.com/ianhussey/port) R package.
All computation is done by port; this app only handles input, layout, and
download.

## Panels

- **Correlation matrix** — upload a reported matrix, get the verdict
  (`check_corr_psd()`), the excusable rounding half-width (`excusable_delta()`),
  the fault localization (`localize_psd_fault()`), per-cell feasible intervals
  (`implied_interval()`), and rounding-robust inconsistent triples
  (`inconsistent_triples()`). Results download as a multi-sheet `.xlsx`.
- **Three correlations** — the same machinery on a single triple of variables.
- **About** — the method, what an inconsistent verdict does and does not prove,
  and the citation.

## Input format

A `.csv` or `.xlsx` file with one row and one column per variable, in the same
order. Accepted shapes:

- the full symmetric matrix;
- either triangle, with the other half left blank;
- an optional leading column of variable names (column headers are used for the
  variable names when they are informative).

The diagonal is set to 1 whatever is supplied. A cell left blank in *both*
triangles is treated as unreported: port frees it to `[-1, 1]`, so a verdict
then holds whatever the unreported value was. `example_correlation_matrix.csv`
is a worked example (it is inconsistent at 2 decimal places).

## Running it

```r
shiny::runApp("app/shiny_app")
```

Requires `shiny`, `bslib`, `port`, plus `readxl` (Excel upload) and `writexl`
(results download).

## Deploying

```r
rsconnect::deployApp("app/shiny_app")
```

`port` must be installable on the server (from CRAN, or via
`remotes::install_github("ianhussey/port")`).
