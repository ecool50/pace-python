# PACE

**Proximity-associated changes in expression** for imaging spatial
transcriptomics — Xenium, CosMx, MERFISH and the like.

PACE asks a directed question of spatial data: *when a macrophage sits next to
tumour cells, which of its genes change, and by how much?*

That question is harder than it looks. Two cells of the same type sitting in
different neighbourhoods differ for reasons that have nothing to do with their
neighbours, and in imaging data a cell's expression profile is partly its
neighbours' transcripts, misassigned across the segmentation boundary. PACE
handles both. It fits a hierarchical negative binomial mixed model in which
every focal–neighbour relationship borrows strength from the others, so an
isolated noisy slope is pulled towards the consensus instead of being reported
as a discovery, and a per-cell contamination term absorbs the short-range
ambient signal, so what comes back is a change in the cell's own expression.

```sh
pip install pace-spatial          # fitting and every readout
pip install "pace-spatial[all]"   # plus AnnData, pandas tables and plotting
```

> Not released yet. Until the first version is on PyPI, install from a clone
> with `pip install -e ".[all]"` and read **Building from source** first.

`numpy` and `scipy` are all a fit needs. `anndata`, `pandas` and `matplotlib`
are optional, imported only when you reach for the thing that uses them.

## Quickstart

One call fits the model and every reporting layer:

```python
import pace

fit = pace.pace_fit(adata, cell_type="cell_type")
```

The readouts then come off the fit, as attributes or as functions:

```python
pace.neighbour_slopes(fit)        # every slope, with a local false sign rate
pace.variance_decomposition(fit)  # where each gene's variance comes from
pace.top_drivers(fit, "Macrophage", "Tumour")
pace.pair_variance(fit)           # which neighbours matter to which types
pace.cell_contamination(fit)      # per cell, how much is its neighbours'
```

`pace.top_drivers(fit, "Macrophage", "Tumour").to_dataframe().head()` on a
breast cancer Xenium section:

```
 rank   gene  driver_score  b_clean   spec   lfsr
    1   MRC1        0.0197  -0.1495 0.9358 0.0000
    2  APOC1        0.0057   0.1302 0.8122 0.0000
    3  CD163        0.0020  -0.0391 0.9474 0.0365
```

Macrophages near tumour lose MRC1 and gain APOC1: the tissue-resident marker
falling and the lipid-associated, immunosuppressive one rising. Because the
contamination term has already absorbed the tumour transcripts bleeding into
those macrophages, these are macrophage genes changing, not tumour genes
leaking. `examples/demo.py` runs the whole pipeline end to end.

From plain arrays instead of an `AnnData`:

```python
fit = pace.pace_fit(counts, x, y, cell_type)   # counts is cells by genes
```

### The four stages

`pace_fit()` runs these in turn. Call them yourself to stop partway, or to
re-run a later stage with different settings:

```python
fit = pace.pace_model(adata, cell_type="cell_type")   # the model alone
fit = pace.pace_shrink(fit, null_correlation=False)
fit = pace.pace_decompose(fit, counts)
fit = pace.pace_drivers(fit)
```

Each takes a fit and returns it with one more layer attached.

## The readouts

**`pace_model` — the model.** One slope per
(focal type, neighbour type, gene), with a per-cell contamination loading and
gene-wise overdispersions. Given an `AnnData` it reads coordinates from
`obsm["spatial"]`, labels and condition from `obs`, gene names from
`var_names` and counts from `X` or a `layer`, and never modifies it.

**`pace_shrink` — which slopes to believe.** Raw slopes are mostly noise.
Multivariate adaptive shrinkage pools them, one mixture per neighbour type
with the focal types as its conditions, and returns every slope with a local
false sign rate. `slopes.calls()` reports how many effects are called and how
many of those are expected to have the wrong sign.

**`pace_decompose` — where a gene's variance comes from.** Per focal type and
gene, split into cell-type identity, the baseline spatial effect, the
condition-by-space interaction, contamination, and dispersion. Useful for the
question "is this gene's variation about *where* the cell is at all?"
**`single_frame_decomposition`** gives the observed version, where the
magnitudes are measured from the data and only the split within a cell type
comes from the model.

**`pace_drivers` — which genes carry a relationship.** Per pair, genes ranked
by `b² · spec² · focal_mean`. Squaring the slope makes the score
direction-blind, so a gene pushed down ranks alongside one pushed up;
specificity and expression level keep a barely-expressed gene off the top on a
noisy slope. Filtered at lfsr < 0.05.

**`pair_variance` — which neighbours matter to which cell types.** Each focal
type's spatial variance split across its neighbours, so a row sums to that
type's spatial share and the own-type pair is excluded.

**`plot_*`** — five figures over the above.

Three more: `build_neighbourhood()` and `ambient_field()` return the kernels
and the contamination field without fitting, for checking what the model is
being given; `anchor_genes()` lists the genes that mark each type and so
anchor its contamination; `stash_pace()` copies the readouts into
`adata.uns`, and `compact_fit()` drops the rebuildable bulk from a fit you
want to keep.

## Figures

![PACE figures](examples/pace_demo.png)

```python
pace.plot_decomposition(frame)                 # full bar + the rescaled zoom
pace.plot_pair_heatmap(pairs, fit.cell_types)
pace.plot_drivers(drivers, "Macrophage", "Tumour", single_frame=frame)
pace.plot_proximity(fit, counts, ["MRC1", "APOC1"], "Macrophage", "Tumour")
pace.plot_shrinkage(slopes)
```

Two are worth knowing about. `plot_decomposition` leaves the residual block
out, so its bars stop well short of 100%: the comparison is between the parts
the model attributes, and a residual bar towering over the rest would tell you
nothing. `plot_drivers` takes `single_frame=` to add a companion panel showing
what each top gene's variance is made of, which catches a gene that ranks
highly because the focal type simply expresses a lot of it.

Every panel returns matplotlib axes, so saving is ordinary matplotlib. For
publication, write PDF with Type 42 fonts, which keeps text editable and avoids
the Type 3 fonts journals reject:

```python
import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
axis, zoom = pace.plot_decomposition(frame)
axis.figure.savefig("decomposition.pdf", bbox_inches="tight")
```

## Things worth setting

**Give it counts.** `adata.X` holds whatever the last preprocessing step left
there, and a negative binomial fit of log-normalised values will converge and
report slopes that mean nothing. Non-integer values are refused; point `layer`
at the raw counts.

**Name the reference arm** of a condition. Its first level is what every other
arm is contrasted against, so it decides what the baseline slopes mean and
which way the interaction points — this is not a relabelling. Leaving
`condition_levels` unset picks one for you and warns.

```python
fit = pace.fit_anndata(adata, "cell_type", condition_key="Responder",
                       condition_levels=["non-responder", "Responder"])
```

**`kernel_per_image`** is `False`: the neighbour search runs across the whole
object even when it has several images, which is what fields of view from one
section want. Set it `True` when the images are genuinely separate samples.

**`types=`** fixes the cell type order, and so the row order of everything
downstream. The default sorts the observed labels.

**`ambient_mode="stream"`** avoids materialising the cells-by-genes
contamination field, which runs to several GB at a million cells. It is a
memory setting and nothing else — the fit is identical either way.

**On a condition fit, `pace_shrink` returns the responder layer**, the
condition-by-neighbour interaction, not the baseline proximity effect.

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

## Building from source

Two things are not free choices.

**Eigen must be 3.4.** Later majors vectorise and reduce in a different order,
which moves the fitted effects by about 3e-4, so the build refuses them:

```sh
bash tools/fetch_eigen.sh && EIGEN_INCLUDE_DIR=$PWD/eigen-3.4.0 pip install -e .
```

`PACE_ALLOW_ANY_EIGEN=1` builds anyway. Installing a released wheel avoids the
question entirely.

**On macOS, use Apple's toolchain.** A conda or Homebrew `clang` brings an old
`ld64` that links a shared object which loads nowhere — the build succeeds and
it fails at import with `mis-aligned LINKEDIT string pool`:

```sh
CC=/usr/bin/clang CXX=/usr/bin/clang++ \
  LDSHARED="/usr/bin/clang -bundle -undefined dynamic_lookup" pip install -e .
```

`ld -v` should report `PROJECT:ld-27037` or newer.

## Also available in R

PACE is also a [Bioconductor package](https://github.com/ecool50/PACE). Both
compile the same C++ core, so a fit is the same fit in either language — on a
matched build, identical to the last bit. If you work in both, note that the
cell type order follows the platform's sort, so pass `types=` explicitly when
you want the two to line up row for row.

The core is owned by the R repository and vendored here under `src/core/`,
synced with `tools/sync_core.py`, which records a checksum manifest and the
upstream commit. Change it there and sync; never edit `src/core/` here.

## Licence

GPL-3.0-or-later.
