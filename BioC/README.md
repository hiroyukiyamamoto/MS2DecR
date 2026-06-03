# MS2DecR Bioconductor preparation notes

This folder keeps a separate Bioconductor-oriented copy of MS2DecR and notes
for a possible future Bioconductor submission.

The project is not being submitted immediately. The goal is to keep the current
status and remaining tasks clear so that the package can be uploaded to GitHub
and improved gradually.

## Folder layout

The original package directory is intentionally kept separate from the
Bioconductor preparation copy.

```text
MS2DecR-main/
  DESCRIPTION          # original package metadata, currently Version 0.1.0
  R/
  man/
  demo/
  BioC/
    README.md          # this note
    MS2DecR-bioc/      # Bioconductor-oriented working copy
```

Use `BioC/MS2DecR-bioc` for future Bioconductor-specific edits. When the
Bioconductor version is ready, selected changes can be moved back into the main
package.

## Current status

As of the current preparation pass, the Bioconductor-oriented copy could be
built and checked with standard R package checks:

```text
R CMD build MS2DecR-bioc
R CMD check --no-manual --no-build-vignettes MS2DecR_0.99.0.tar.gz
Status: OK
```

The checked source tarball was generated as:

```text
C:/Users/hyama/Documents/R/MS2DecR/release/MS2DecR_0.99.0.tar.gz
```

Note: this tarball was generated during local checking. It does not need to be
committed to GitHub.

## Changes already made

- Updated `DESCRIPTION` for a Bioconductor-style development version:
  - `Version: 0.99.0`
  - added `biocViews: MassSpectrometry, Metabolomics, Proteomics, Software`
  - cleaned up imports and suggests
- Replaced broad `exportPattern()` in `NAMESPACE` with explicit exports.
- Added a missing R implementation for `match_spectra()`.
- Reworked demos so they use the included `msdial` dataset instead of local
  absolute file paths.
- Added `demo/00Index`.
- Rewrote demo files as valid UTF-8 without BOM.
- Fixed documentation mismatches in Rd files.
- Added documentation for `matchedfilter()`.
- Moved examples requiring external mzML/MSP files into `\dontrun{}` blocks.
- Removed the broken dependency on `loadings::ica`, which was not available in
  the installed `loadings` namespace.
- Temporarily changed `performICA()` initialization to use `stats::prcomp()`
  so the package can pass standard checks in the current environment.

## Important open decisions

### License

The current package license is:

```text
License: CC BY-NC-ND 4.0
```

This is likely not suitable for Bioconductor because it restricts commercial
use and derivative works. Before submission, choose an open-source license that
allows redistribution and modification.

Likely candidates:

- `Artistic-2.0` - common and conservative for Bioconductor packages
- `GPL-3`
- `MIT + file LICENSE`

Do not change this casually unless the author is comfortable with the legal
meaning of the new license.

### ICA implementation

The original code called `ica()`, apparently expecting an ICA implementation
from another package. In the current local environment, `loadings::ica` did not
exist, so the package could not pass examples.

Current temporary behavior:

```text
performICA() -> internal PCA-based initialization using stats::prcomp()
```

Before submission, decide whether MS2DecR should:

- keep the PCA-based initialization and adjust the wording away from strict ICA,
- import a stable ICA implementation such as `fastICA`, or
- implement a small internal ICA routine if scientifically justified.

If the package continues to describe the workflow as ICA-based, using an actual
ICA backend is preferable.

## Bioconductor tasks for later

1. Install and run `BiocCheck`.

   Suggested future check:

   ```r
   BiocManager::install("BiocCheck")
   BiocCheck::BiocCheck("MS2DecR_0.99.0.tar.gz", `new-package` = TRUE)
   ```

2. Add a formal vignette under `vignettes/`.

   Suggested file:

   ```text
   vignettes/MS2DecR.Rmd
   ```

   It should use the included `msdial` dataset and demonstrate:

   - generating or loading an MS2 matrix
   - ICA/ALS or final chosen deconvolution workflow
   - MS-DIAL-like model construction
   - MS-DIAL-like deconvolution
   - spectral matching
   - interpreting the returned objects

3. Add tests.

   Suggested structure:

   ```text
   tests/testthat.R
   tests/testthat/test-deconvICA.R
   tests/testthat/test-msdial-model.R
   tests/testthat/test-msdial-deconv.R
   tests/testthat/test-match-spectra.R
   ```

4. Add `NEWS.md`.

   Start with:

   ```text
   # MS2DecR 0.99.0

   - Initial Bioconductor submission candidate.
   ```

5. Review package contents.

   Check whether these should be included in the final submission:

   - top-level `MS2DecR.png`
   - `dev/`
   - any large example data
   - old source tarballs inside the package directory

   Development-only files should generally stay outside the package tarball or
   be excluded with `.Rbuildignore`.

6. Add or update metadata before GitHub submission.

   Recommended `DESCRIPTION` fields:

   ```text
   URL: https://github.com/<user>/MS2DecR
   BugReports: https://github.com/<user>/MS2DecR/issues
   ```

7. Make sure the package is not already submitted to CRAN.

   Bioconductor packages should not also be on CRAN.

8. Submit through the Bioconductor Contributions repository when ready.

   General route:

   - push the package to GitHub
   - confirm `R CMD check` is clean
   - run and address `BiocCheck`
   - open a new issue at Bioconductor Contributions
   - respond to reviewer comments

## Useful commands

From:

```text
C:/Users/hyama/Documents/R/MS2DecR/release/MS2DecR-main/BioC
```

Build:

```text
"C:/Program Files/R/R-4.3.1/bin/R.exe" CMD build MS2DecR-bioc
```

Check:

```text
"C:/Program Files/R/R-4.3.1/bin/R.exe" CMD check --no-manual --no-build-vignettes MS2DecR_0.99.0.tar.gz
```

## Notes before uploading to GitHub

- Confirm whether generated tarballs such as `MS2DecR_0.99.0.tar.gz` should be
  committed. Usually source tarballs are not committed to the package repository.
- Add a `.gitignore` / `.Rbuildignore` if needed.
- Keep this `BioC/README.md` as planning notes, or move the important parts to
  GitHub issues after upload.
