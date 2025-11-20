## Resubmission

This is a resubmission of SSP 1.1.0, addressing CRAN feedback.

We have:

- Removed the hidden `.github/` directory from the source tarball.
- Reduced example runtimes to remain under 5 seconds each.
- Preserved informative usage via the vignette, referenced in example help pages.

## Test environments

- Windows 10 x64, R 4.4.0 (ucrt), local, `--as-cran`
- Ubuntu 24.04 (via GitHub Actions), R 4.5.2, `--as-cran`
- Ubuntu 22.04 (via rhub), R-devel (2025-11-12 r89009)
- Windows (via rhub), R-devel (2025-11-12 r89009)

R CMD check results:
- 0 ERRORs ✓
- 0 WARNINGs ✓
- 0 NOTEs ✓

## Reverse dependencies

Only one package (**ecocbo**) depends on **SSP**, and it is maintained by the same team. Compatibility is confirmed.

## Additional comments

- All examples now execute in <5 seconds.
- The vignette provides a complete demonstration of multi-site simulations that were previously shown in examples.
- GitHub CI workflows were tested but excluded from the final package per CRAN guidance.
