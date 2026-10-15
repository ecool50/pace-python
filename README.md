# PACE for Python

**Proximity-associated changes in expression**, for imaging spatial
transcriptomics such as Xenium and CosMx.

PACE asks a directed question of spatial data: *when a macrophage sits next to
tumour cells, which of its genes change, and by how much?* It fits a
hierarchical negative binomial mixed model in which every focal-neighbour
relationship borrows strength from the others, so noisy pairs are pulled
towards the consensus rather than reported as discoveries. A per-cell
contamination term absorbs the short-range ambient signal that leaks between
adjacent cells through segmentation error, so what comes back is a change in
the cell's own expression rather than its neighbour's transcripts landing in
it.

This is the Python implementation. It is not a reimplementation: it compiles
the same C++ core as the [Bioconductor package](https://github.com/ecool50/PACE)
and reproduces its fits bit for bit.

```sh
pip install pace-spatial          # fitting and every readout
pip install "pace-spatial[all]"   # plus AnnData, pandas tables and plotting
```

> Not released yet. Until the first version is on PyPI, install from a clone:
> `pip install -e ".[all]"`, and read **Building from source** below first —
> the Eigen version is not a free choice.

`numpy` and `scipy` are all a fit needs. `anndata`, `pandas` and `matplotlib`
are optional, imported only when you reach for the thing that uses them.

## Quickstart

```python
import pace

fit    = pace.fit_anndata(adata, "cell_type")
slopes = pace.pace_shrink(fit)
stats  = pace.fit_statistics(fit, counts)
blocks = pace.pace_decompose(fit, counts, stats=stats)
frame  = pace.single_frame_decomposition(fit, counts, blocks)
drivers = pace.pace_drivers(fit, slopes, blocks, stats)
pairs  = pace.pair_variance(fit, counts, frame)

drivers["Macrophage", "Tumour"].to_dataframe().head()
```

On the breast cancer Xenium crop the R package ships, that last line gives:

```
 rank   gene  driver_score  b_clean   spec   lfsr
    1   MRC1        0.0197  -0.1495 0.9358 0.0000
    2  APOC1        0.0057   0.1302 0.8122 0.0000
    3  CD163        0.0020  -0.0391 0.9474 0.0365
```

MRC1 down and APOC1 up in macrophages near tumour: the tissue-resident marker
falling and the lipid-associated one rising, which is the finding the method
was built to recover. `examples/demo.py` runs the whole thing end to end.

Fitting from plain arrays instead of an `AnnData`:

```python
fit = pace.pace_fit(counts, x, y, cell_type)   # counts is cells by genes
```

## What you get

**`pace_fit` / `fit_anndata`** — the model. `fit_anndata` reads coordinates
from `obsm["spatial"]`, labels from `obs`, gene names from `var_names` and
counts from `X` or a `layer`, and never modifies the object.

**`pace_shrink`** — multivariate adaptive shrinkage of the slopes, one mash fit
per neighbour with the focal types as its conditions. Returns every slope with
a local false sign rate; `slopes.calls()` summarises how many effects are
called and how many of those should have the wrong sign.

**`pace_decompose`** — each gene's fitted variance per focal type, split into
cell-type identity, baseline spatial state, responder spatial state,
contamination and dispersion. **`single_frame_decomposition`** gives the
observed version, where the magnitudes come from the data and only the split
within a cell type comes from the model; that is the table R's
`varianceDecomposition()` returns.

**`pace_drivers`** — per-pair tables ranking the genes that mediate a
relationship by `b² · spec² · focal_mean`. Squaring the slope makes the score
direction-blind, so a gene pushed down ranks alongside one pushed up;
specificity and expression level keep a barely-expressed gene off the top on a
noisy slope. Filtered on lfsr < 0.05 alone.

**`pair_variance`** — how much of each focal type's variance each neighbour
contributes. The focal's spatial share comes from the observed single frame and
is split across neighbours by a Pratt attribution renormalised over the
off-diagonal pairs, so a row sums to that share and the own-type pair is
excluded.

**`plot_*`** — the five figures, reproducing the R package's panels down to
their conventions and hex colours.

## Figures

![PACE figures](examples/pace_demo.png)

```python
pace.plot_decomposition(frame)                 # full bar + the rescaled zoom
pace.plot_pair_heatmap(pairs, fit.cell_types)
pace.plot_drivers(drivers, "Macrophage", "Tumour", single_frame=frame)
pace.plot_proximity(fit, counts, ["MRC1", "APOC1"], "Macrophage", "Tumour")
pace.plot_shrinkage(slopes)
```

Every panel returns matplotlib axes, so saving is ordinary matplotlib. For
publication use PDF with Type 42 fonts, which keeps text editable and avoids
the Type 3 fonts journals reject:

```python
import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
axis, zoom = pace.plot_decomposition(frame)
axis.figure.savefig("decomposition.pdf", bbox_inches="tight")
```

Four of the conventions are R's and are not obvious from the data: the
decomposition leaves the **residual block out**, so bars stop near 40% rather
than 100, and orders focal types by descending spatial share; the pair heatmap
puts **focal on x**; `plot_drivers` takes `single_frame=` to add R's companion
panel, which catches a top gene whose variance is really cell-type identity;
and `plot_proximity` bins by the **number of neighbours** within a radius, at
R's breaks of 0, 2, 4, 6, 8 and up, because the fitted kernel is centred within
image and would give a negative axis.

## How closely it matches R

Bit-identity is a property of a build, not of the package. CI measures all
three cases against the R package's own fits:

| build | agreement |
|---|---|
| arm64 macOS, Apple clang (the reference) | **bit-identical** |
| arm64 Linux, g++ | every tolerance test passes; exact equality lost |
| x86-64 Linux, g++ | about 2.7e-4 |

Eigen vectorises and reduces in a different order on x86-64, and thirteen IRLS
iterations amplify it. One further difference is deliberate: R takes `trigamma`
from its own maths library where the core computes its own. They agree to
1.8e-15, which those iterations amplify to 8e-07 in the fitted effects; against
an R build using the core's trigamma, nothing differs at all.

The readouts are held to the same standard. Shrinkage calls the **same
gene-focal pairs** as R, not merely the same number. Driver tables give the
same genes in the same rank order, to 2e-7. The observed decomposition's
sums of squares match to 1e-11, since they come from the counts rather than
the fit.

## Before comparing to an R fit row by row

**Counts must be counts.** `adata.X` is whatever the last preprocessing step
left there, and a negative binomial fit of log-normalised values converges and
reports slopes that mean nothing. Non-integer values are refused; point `layer`
at the raw counts.

**Cell type order.** `types=None` sorts labels in Python's byte order; R sorts
under `LC_COLLATE`, so a mixed-case label set orders differently. The fit is
equivariant under a permutation, so the biology is unchanged but the row order
is not. Pass `types=` to match a given R fit.

**The condition's reference arm is not a relabelling.** Its first level is what
every other arm is contrasted against, so it decides what the baseline slopes
mean and which way the interaction points. R picks it by collation and Python
would pick it by byte order; on labels as ordinary as `{"Responder",
"non-responder"}` the fitted slopes then differ by 2. Leaving
`condition_levels` unset warns for that reason. A numeric condition is refused,
since R would model it continuously.

**Neighbour scope.** `kernel_per_image=False`, as R's is: the search is global
even across images, which is what fields of view from one section want. Set it
`True` when the images are separate samples.

**Shrinkage settings.** The core has mashr's canonical covariances but not the
data-driven ones, so this is R at `data_driven = FALSE` — worth 11 of 1,637
breast cancer calls. `null_correlation` defaults to `True` as R's does; the
published counts used `False`. On a condition fit what is shrunk is
`<resp_term>:<neighbour>`, so the result is the **responder layer, not the
baseline effect**.

**Large cohorts.** `ambient_mode="stream"` avoids materialising the
cells-by-genes ambient product, which runs to several GB at a million cells. It
is a memory knob only: the streamed fit is bit-identical to the cached one.

## Tests

```sh
pip install -e ".[test]"
python -m pytest tests
```

From a non-editable install, run them from outside the checkout
(`cd /tmp && python -m pytest /path/to/pace-python/tests`): the repository root
holds the source package `pace/`, which shadows the installed one but has no
compiled `_core`, and the resulting ImportError blames a circular import rather
than the real cause.

The exact-equality tests are marked `bit_identical` and skip off the reference
build; `PACE_REFERENCE_BUILD=1` forces them anywhere. The mash tests hold the
bindings to every neighbour slice of both canonical cohorts, compared gene by
gene against mashr: 1,634 breast cancer and 426 melanoma calls.

## Building from source

Two things are not free choices.

**Eigen must be 3.4.** Against 3.5 the random effects move by about 3e-4 and
the two languages no longer agree, so the build refuses anything else:

```sh
EIGEN_INCLUDE_DIR=$(Rscript -e 'cat(system.file("include", package="RcppEigen"))') \
  pip install -e .
```

or fetch it directly with `bash tools/fetch_eigen.sh`. Set
`PACE_ALLOW_ANY_EIGEN=1` to build anyway, accepting that the fits will not
match. Installing a released wheel avoids the question entirely.

**On macOS, use Apple's toolchain.** A conda or Homebrew `clang` brings an old
`ld64` that links a shared object which loads nowhere — the build succeeds and
the failure appears at import as `mis-aligned LINKEDIT string pool`:

```sh
CC=/usr/bin/clang CXX=/usr/bin/clang++ \
  LDSHARED="/usr/bin/clang -bundle -undefined dynamic_lookup" pip install -e .
```

`ld -v` should report `PROJECT:ld-27037` or newer. The R package has the same
requirement for the same reason.

## The shared core

The C++ core is owned by the [R package](https://github.com/ecool50/PACE) and
vendored here under `src/core/`. Keeping the two in step is the one risk of
separate repositories, so it is made visible rather than trusted:

```sh
python tools/sync_core.py ../PACE    # copy it in, record what came
python -m pytest tests               # and prove it still reproduces R
```

`sync_core.py` writes a sha256 manifest with the upstream commit;
`tests/test_core_sync.py` fails if a vendored file is edited here. The guard
that matters is the rest of the suite, whose fixtures are R's own fits. Change
the core in the R repository and sync — never edit `src/core/` here.

## Licence

GPL-3.0-or-later, as the R package is.
