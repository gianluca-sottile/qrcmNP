## Submission

This is a new major version of 'qrcmNP' (0.2.2 -> 0.3.0):

* the package no longer accesses non-exported functions of 'qrcm'
  (the required internal routines are now part of 'qrcmNP', with attribution);
* 'qrcm' and 'survival' moved from Depends to Imports ('slp()' and 'plf()' are re-exported);
* new penalized estimator (lasso/group lasso, proximal Newton algorithm), new
  interface and inference for nonlinear models; see NEWS.md.

## Test environments

* local macOS, R release
* GitHub Actions: macOS, Windows, Ubuntu (devel, release, oldrel-1)
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies on CRAN.
