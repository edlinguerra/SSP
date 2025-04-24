## Resubmission

This is a resubmission of the SSP package (version 1.0.2) to CRAN, following the request to fix cross-references in Rd files.

### Changes in this version:

- Updated Rd cross-references to use explicit package anchors (e.g., `\link[ggplot2:ggplot2-package]{ggplot2}`) as required.
- Migrated all documentation to roxygen2 and ensured consistent formatting.
- Fixed undocumented arguments and cleaned up duplicate documentation tags.
- Removed deprecated use of `rel()` from `grid` package and replaced it with equivalent numeric sizes.
- Ensured all examples respect line length limits (max 100 characters).
- Moved legacy `.Rd` backups to a non-installed directory and excluded them via `.Rbuildignore`.
- Fixed missing `@importFrom` statements for all used functions (e.g., `ggplot`, `geom_text`, etc.).

## Test environments

* Local Windows 10 x64, R 4.4.0 (ucrt)
* Local R 4.5.0 (ucrt) (with `--as-cran`)
* Checked using `devtools::check()` and `R CMD check`

## R CMD check results
Copy/Paste of the R CMD check results

❯ checking for future file timestamps ... NOTE
  unable to verify current time

0 errors ✔ | 0 warnings ✔ | 1 note ✖

## Downstream dependencies

There are currently no known downstream packages depending on SSP.
