# SSP 1.1.0

## Function adjustments

- `simdata()` can now introduce mild randomness via the new `jitter.base` 
  parameter to increase variability among simulated communities. This feature is 
  optional (default `jitter.base = 0`).


# SSP 1.0.2 (2025-04-23)

## Changes requested by CRAN

- Updated cross-references in `.Rd` files to use explicit package
  anchors (e.g., `\link[ggplot2:ggplot2-package]{ggplot2}`) as requested
  by the CRAN team.
- Removed all `\link{}` usages without package anchors in documentation
  files.
- Ensured compatibility with HTML reference manuals.

## Documentation improvements

- All .Rd documentation files were migrated to the **roxygen2** format.
- Improved consistency and accuracy in the documentation of function
  arguments and outputs.
- Added missing imports in the NAMESPACE to fix visibility issues
  (`ggplot`, `geom_point`, `specnumber`, etc.).
- Removed unused `@importFrom` tags and documented only necessary
  dependencies.

## Code adjustments

- Replaced use of `rel()` in `plot_ssp()` to avoid dependency on
  `grid::rel`, which is not exported.
- Reformatted long lines in examples to comply with CRAN’s 100-character
  limit.

## Testing and validation

- `devtools::check()` returns **0 errors \| 0 warnings \| 1 notes**
  (only minor ones related to timestamp verification).
