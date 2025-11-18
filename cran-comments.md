## Resubmission

This is a resubmission of SSP 1.1.0, following incoming checks.

We have:

- Updated the `Date:` field in DESCRIPTION.
- Kept the `.github/` folder intentionally to support automated R-hub checks via GitHub Actions.

### Test environments

* Windows 10 x64, R 4.4.0 (ucrt)
* Windows 10 x64, R 4.5.0 (ucrt), R CMD check --as-cran
* Ubuntu 22.04 (via rhub), R-devel (2025-11-12 r89009)
* R CMD check results: 0 ERRORs | 0 WARNINGs | 1 NOTEs (see below)

Found the following hidden files and directories:
    .github
These were most likely included in error. See section 'Package structure' in the 'Writing R Extensions' manual.
  
### Reverse dependencies

Only **ecocbo** depends on **SSP**; it is maintained by the same team and remains compatible.

### Additional comments

- The `.github/` folder is included deliberately for CI workflows using [r-hub/actions](https://github.com/r-hub/actions).
- The `examples()` NOTE reflects simulation-based functions (`simdata()`, `plot_ssp()`, etc.) that require >10s to run, even with minimal input. These are essential to demonstrate usage and cannot be shortened meaningfully.
