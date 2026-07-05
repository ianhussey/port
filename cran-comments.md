## Submission

This is a new submission.

## Test environments

- local macOS, R release
- (add: win-builder devel/release, R-hub, GitHub Actions before submitting)

## R CMD check results

0 errors | 0 warnings | 1 note (new submission).

## Notes for the reviewer

* The package uses domain terms and author names that a spell-checker does
  not recognise but are spelled correctly and intentional: PSD (positive
  semidefinite), semidefinite, disattenuation, metascience, and the cited
  authors Higham and Rump. These are recorded in `inst/WORDLIST`, so
  `spelling::spell_check_package()` and the incoming spell-check pass with no
  spelling NOTE. The two references are given as `Authors (year) <doi:...>`.
* Long-running Monte-Carlo property tests are scaled down to a fast smoke on
  CRAN (via the `NOT_CRAN` environment variable) and run in full locally and
  in CI, so the check stays well within time limits.
* The package has no compiled code and writes no files outside `tempdir()`.
