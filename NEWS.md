# port 0.1.0

First release.

`port` decides whether a reported, rounded Pearson correlation matrix is consistent with any positive semidefinite (PSD) matrix inside its rounding box, and grades and localizes the violation when it is not.

## Core check

* `check_corr_psd()` returns a rounding-box verdict of `inconsistent`, `consistent`, or `undecided`, backed by a floating-point-safe witness bound (Higham 2002) on the inconsistency side and a Rump-verified in-box PSD matrix (Rump 2006) on the consistency side, with no dependency on an external SDP solver.
* Every result carries a `certified` flag, and `certificate()` extracts the witness vector.
* Supports asymmetric rounding rules (`rounding = "truncate"/"floor"/"ceiling"`), per-cell precision via a `decimals`/`delta` matrix, and missing (`NA`) cells freed to `[-1, 1]`.
* `check_corr_psd_batch()` screens a list of matrices into a tidy tibble, including the off-diagonal correlation summary (`r_min`, `r_max`, `r_mean`, `r_sd`).

## Severity and localization

* `excusable_delta()` folds the verdict and its severity into one number in precision units.
* `localize_psd_fault()` classifies where the fault sits (`cell`, `cell_tentative`, `triad`, `variable`, `joint`, `diffuse`, or `none`) and separates benign structural causes from substantive violations; `localize_psd_fault_batch()` runs it over a list.
* `fault_evidence()` exposes the raw box-aware evidence, `inconsistent_triples()` the standalone triple scan, and `implied_interval()` the interval a cell (including a missing one) could occupy.

## Disattenuation

* `check_disattenuated_psd()` checks reliability-corrected construct matrices, adding a `consistent but implausible` level, and reports the critical reliability when reliabilities are unreported; `disattenuate()` applies the Spearman correction.
