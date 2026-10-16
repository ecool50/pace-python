# PACE

A hierarchical empirical Bayes framework for **proximity-associated changes in
expression** — how a cell's gene expression shifts near specific neighbouring
cell types — in imaging-based spatial transcriptomics (Xenium, CosMx).

```{image} figures/pace_overview.png
:alt: Overview of the PACE framework
:width: 100%
```

## Overview

**PACE** quantifies how proximity to specific neighbouring cell types reshapes
a cell's gene expression. It fits hierarchical negative binomial mixed models
with partial pooling across cell types, corrects for transcript contamination
between adjacent cells with a per-cell ambient term, and decomposes expression
variance into cell-type identity, spatial cell state, contamination and
residual components. A per-gene driver score ranks the genes that mediate each
focal–neighbour relationship.

Input and output follow [AnnData](https://anndata.readthedocs.io) conventions:
a single call to {func}`~pace.pace_fit` returns a `PaceFit`, read out with
{func}`~pace.neighbour_slopes`, {func}`~pace.variance_decomposition` and
{func}`~pace.top_drivers`.

```python
import pace

fit = pace.pace_fit(adata, cell_type="cell_type")
pace.top_drivers(fit, "Macrophage", "Tumour").to_dataframe().head()
```

## Installation

PACE needs a C++ compiler and Eigen 3.4 for the contamination-correcting
solver. Installing a released wheel avoids both.

```sh
pip install pace-spatial
```

This is the Python implementation. PACE is also an
[R/Bioconductor package](https://github.com/ecool50/PACE); both compile the
same C++ core, so a fit is the same fit in either language.

```{toctree}
:maxdepth: 2
:hidden:

usage
demo
api
```
