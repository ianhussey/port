## Submission

This is a new submission.

## Test environments

- local macOS (aarch64), R 4.5.2 (release)
- win-builder (x86_64-w64-mingw32), R-devel, via `devtools::check_win_devel()`

## R CMD check results

0 errors | 0 warnings | 1 note.

The one NOTE is the standard "New submission" flag; there are no other notes on
CRAN's machines. (A local `R CMD check --as-cran` additionally reports an
HTML-Tidy version note, which is specific to the local `tidy` install and does
not occur on the CRAN check farm.)

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
