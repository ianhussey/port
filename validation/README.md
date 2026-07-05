# psdness validation simulations

Three ADEMP-structured Monte Carlo simulations that validate the operating
characteristics of the `psdness` package (`check_corr_psd()` /
`localize_psd_fault()`). They follow the template in
`understanding-statistics-through-monte-carlo-simulations` (Morris et al. 2019;
Siepe et al. 2024): each `.qmd` has the full ADEMP introduction, a
`tidyr` + `furrr` + `simhelpers` pipeline, and Monte Carlo uncertainty reported
alongside every estimate.

## What each simulation validates

| File | Property | Pass condition |
|---|---|---|
| `validation_sim1_soundness_exact.qmd` | **Soundness on exact valid matrices.** The tool never certifies a genuinely valid (PSD) matrix as impossible. | False-impossibility rate **= 0** (Wilson-bounded). A hard diagnostics chunk halts the render on any violation. |
| `validation_sim2_soundness_rounded.qmd` | **Soundness under rounding**, and the delta-matching boundary condition. Rounding a valid matrix to nearest with the matching symmetric `delta` never yields impossible; an asymmetric rule with a mismatched (too-small) `delta` can — correctly for the box it was given — and a matched `delta` restores soundness. | False-impossibility rate **= 0** under any *matched*-delta condition (hard assertion). Unmatched asymmetric rules may fire; the oracle confirms those boxes are genuinely infeasible. |
| `validation_sim3_power.qmd` | **Power + localization + severity.** Detection power on genuinely infeasible boxes, culprit localization, severity calibration, and a false-positive guard. | Power → 1 as perturbation magnitude grows on oracle-infeasible cases; impossibility rate ≈ 0 on oracle-feasible cases; localization recovers the planted culprit; `severity_max` increases monotonically with magnitude. |

Sims 1 and 2 validate the **soundness / zero-false-impossibility** payload; Sim 3
validates **power, localization, and severity**.

## The independent oracle (ground truth)

Ground truth is the answer to the **box** question — *does a PSD matrix exist
inside the rounding box?* — and is computed independently of the tool under test
(`box_is_feasible()` in `R/sim_helpers.R`). It is **not** the point
`lambda_min >= -1e-10` tolerance used elsewhere, which answers a different
(point, not box) question.

- **Infeasible (rigorous):** the Frobenius distance from the reported point to
  the PSD cone exceeds the box radius `delta * sqrt(p(p-1))` by a margin that
  dominates eigenvalue rounding error, so *no* in-box matrix is PSD (triangle
  inequality).
- **Feasible:** an exhibited PSD matrix (the known generating matrix, or an
  independent nearest-correlation matrix) is verified PSD (eigenvalues) and shown
  to lie inside the box (arithmetic).
- **Uncertain:** a thin near-boundary band where neither certificate fires;
  excluded from rate denominators and reported.

## Generators (`R/sim_helpers.R`)

- `gen_onion(p, eta)` — `clusterGeneration` onion method; `eta = 1` spreads over
  valid matrices, `eta = 0.1` piles on the PSD boundary (near-singular). Coverage
  generator and boundary stress source.
- `gen_factor_model(p, k, loading_sd)` — **primary plausibility generator**: a
  positive general factor plus group factors, giving a positive-manifold,
  heterogeneous, moderate-magnitude structure resembling real psychology
  matrices. The off-diagonal distribution is recorded so plausibility is
  auditable.
- `gen_equicorrelation(p, rho)` — closed-form edge-case oracle; unit tests only.
- `gen_independent(p, n)` — trivial sanity row only.

## Running

The package must be installed first (the `furrr` multisession workers load it):

```r
# from the repository root
devtools::install()      # or: R CMD INSTALL .
```

Then render (from this `validation/` directory):

```sh
quarto render validation_sim1_soundness_exact.qmd
quarto render validation_sim2_soundness_rounded.qmd
quarto render validation_sim3_power.qmd
```

- `K = 1000` per condition is the baseline written into each file; a publishable
  zero-false-impossibility claim needs `K = 10000` (stated in each "Number of
  repetitions" section).
- Set the environment variable `SIM_QUICK=1` to render with a small `K` for a
  fast development pass, e.g. `SIM_QUICK=1 quarto render validation_sim1_soundness_exact.qmd`.
- No CVXR / no MATLAB: the tool and the oracle are POCS/eigenvalue-based only.

## Files

- `R/sim_helpers.R` — generators, contamination (`round_offdiag`, `perturb`,
  `contaminate_data`), the independent oracle (`box_is_feasible`), analysis
  wrappers (`analyse_check`, `analyse_localize`), and utilities
  (`offdiag_summary`, `seed_from`, `wilson_ci`).
- `validation_sim1_soundness_exact.qmd`, `validation_sim2_soundness_rounded.qmd`,
  `validation_sim3_power.qmd`.
