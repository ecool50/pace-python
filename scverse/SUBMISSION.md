# scverse ecosystem submission

Open a pull request against
[scverse/ecosystem-packages](https://github.com/scverse/ecosystem-packages)
adding `packages/PACE/meta.yaml` from this directory, with the checklist below
in the PR description.

## Inclusion criteria

- [x] **The code is publicly available under an OSI-approved license.**
      GPL-3.0-or-later, `LICENSE` in the repository root.
- [ ] **The package provides versioned releases.**
      Pending the first tagged release.
- [ ] **The package can be installed from a standard registry.**
      Pending publication of `pace-spatial` to PyPI. Wheels are built by
      `.github/workflows/release.yml` for Linux and macOS on x86-64 and arm64,
      so installing does not require a compiler.
- [x] **Automated tests cover essential functions and a reasonable range of
      inputs.** 166 tests. They are not smoke tests: the fit, the shrinkage,
      both variance decompositions, the driver tables and the pair attribution
      are each held to the R package's own output on a packaged breast cancer
      Xenium crop, and the fit is asserted bit-identical.
- [x] **CI automatically executes these tests on each push or pull request.**
      `.github/workflows/test.yml`: Linux x86-64, Linux arm64 and macOS arm64,
      Python 3.10 and 3.12, six jobs, all green.
- [x] **The package provides API documentation via a website or README.**
      <https://ecool50.github.io/pace-python/> — an overview, a worked
      example, and a reference page per public function generated from the
      docstrings. Built and deployed by `.github/workflows/docs.yml`.
- [x] **The package uses scverse datastructures where appropriate.**
      `pace.fit_anndata(adata, "cell_type")` reads the counts from `X` or a
      layer, coordinates from `obsm["spatial"]`, labels from `obs` and gene
      names from `var_names`, and never modifies the object. `anndata` is an
      optional dependency so that fitting from plain arrays needs only numpy
      and scipy.
- [x] **I am the author or maintainer, and agree to the package being listed.**
- [x] **I agree to abide by the scverse code of conduct.**

## Notes for the reviewers

The compiled core is shared with the Bioconductor package
[PACE](https://github.com/ecool50/PACE) and vendored here under `src/core/`,
synced by `tools/sync_core.py`, which records a sha256 manifest and the
upstream commit. The build pins Eigen 3.4 and refuses any other version:
against Eigen 3.5 the fitted random effects move by about 3e-4, which would
silently break agreement with the R package.
