# port (development version)

* New `check_witness()` evaluates the bound behind an `inconsistent` verdict for a given witness vector, without any search, so a witness from `check_corr_psd()`, from another source or built by hand can be checked directly.

* `check_corr_psd()` certifies tables in which some correlation is pinned at exactly 1 or -1, such as a variable and its duplicate or reverse-coded copy under `rounding = "truncate"`, `"floor"` or `"ceiling"` (or `delta = 0`). Such a box holds no positive definite matrix, so these tables used to be `undecided`. Variables correlated +-1 are equal up to sign in every PSD matrix in the box, so they are now merged, the merged table is certified, and its verified matrix expands exactly to an in-box PSD matrix. The localizers, `implied_interval()` and `check_disattenuated_psd()` use the same step. A new section of `?check_corr_psd` explains why, under nearest rounding, a genuine table is never `undecided` for structural reasons.

* `localize_psd_fault()` now finds joint faults for p >= 7. The search for a small set of cells to correct jointly used to commit cells in index order, so it reported `diffuse` (or a set containing cell (1,2) for no reason) where a two-cell correction existed. It now commits the cell the current witness depends on most.

* `check_disattenuated_psd()` no longer reports `inconsistent` when the smallest attainable corrected correlation is exactly 1. Its intervals are now widened by a floating-point error allowance; before, rounding noise in the "|D| > 1" test flagged some of these valid tables.

* `inconsistent_triples()` (and the triple evidence in `localize_psd_fault()`) flags a triple only when its determinant bound is negative by more than a floating-point error allowance.

* Consistency certificates are checked against the box exactly rather than with a 1e-9 tolerance, and the verified Cholesky test adds a shift covering the few-ulp difference between computed and real box edges.

* In `localize_psd_fault()`, leave-one-out counts a variable as restoring only when its removal gives a certified `consistent`, not a presumed one.

* Tables with missing cells are validated like complete ones: infinite entries are an error, and small asymmetries are symmetrized with a warning instead of being an error.

* `inconsistent_triples()`, `implied_interval()`, `localize_psd_fault()` and `check_disattenuated_psd()` validate `decimals` (and `reliability_decimals`); `implied_interval()` validates `cells`.

* `check_corr_psd_batch()` gains `rounding`, passed to `check_corr_psd()`.

* The witness returned by `check_corr_psd()` is named by the variables when `R` has dimnames.

* The print method shows the general box bound for non-uniform boxes instead of the uniform formula, and the severity labels use the same thresholds as the structural diagnosis.

* The documentation of `excusable_delta()`, the severities and `required_edit` now says that they are certified upper bounds rather than exact values.

* port now requires R >= 4.2.0. Its dependencies already required R >= 4.1.

* `check_corr_psd()` certifies more near-boundary tables as `inconsistent`. When the existing witness directions certify nothing and no in-box PSD matrix is constructed, a greedy backward elimination over principal submatrices searches for a sparse witness and polishes it on its own support. If that finds nothing either and no in-box PSD matrix is found, the direction that separates the box from the PSD cone under alternating projections is tried before the verdict is `undecided`; a certificate found this way has `tier = "pocs"`. Some tables that were `undecided` are now `inconsistent`, and `b_upper` can be smaller than before. Because a smaller `b_upper` can fall below `tau`, some tables that were presumed `consistent` (`certified = FALSE`) now reach the alternating-projection search, which in the cases examined certified them `consistent`. The localizers and `check_disattenuated_psd()` use the sparse search.

* The witness bound for general boxes (used by the localizers, `implied_interval()`, `check_disattenuated_psd()` and `check_corr_psd()` with missing cells, per-cell precision or rounding other than `"nearest"`) no longer builds p x p matrices for every test direction, which took most of their running time. `check_disattenuated_psd()` at p = 70 takes 14 s instead of 33 s. The bound is the same; only floating-point rounding in `b_upper` can differ.

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
