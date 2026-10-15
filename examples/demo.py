"""A full PACE run in Python, on the breast cancer Xenium crop the R package ships.

Builds an AnnData from the packaged crop, fits, and draws the figures. Every
number printed here is produced by the Python package alone.

    python examples/demo.py
"""
import gzip
from pathlib import Path

import anndata
import matplotlib
matplotlib.use("Agg")
# Type 42 embeds the fonts as TrueType, so the text stays selectable and
# editable in Illustrator. The default, Type 3, is what journals bounce.
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix

import pace

FIX = Path(__file__).parent.parent / "tests" / "fixtures"
OUT = Path(__file__).parent


def column(name, dtype=float):
    return np.loadtxt(FIX / (name + ".gz"), dtype=dtype, ndmin=1)


def labels(name):
    with gzip.open(FIX / (name + ".gz"), "rt") as handle:
        return [line.strip() for line in handle if line.strip()]


# --- the data, as a user would hold it -------------------------------------
shape = column("raw_y_dim.txt", int)
counts = csc_matrix(
    (column("raw_y_x.txt"), column("raw_y_i.txt", np.int32), column("raw_y_p.txt", np.int32)),
    shape=(int(shape[0]), int(shape[1])))
genes = labels("raw_genes.txt")
cell_type = labels("raw_celltype.txt")

adata = anndata.AnnData(
    X=counts.tocsr(),
    obs=pd.DataFrame({"cell_type": pd.Categorical(cell_type)},
                     index=[f"cell{i}" for i in range(counts.shape[0])]),
    var=pd.DataFrame(index=genes))
adata.obsm["spatial"] = np.column_stack([column("raw_x.txt"), column("raw_y.txt")])
print(adata)
print("cells per type:")
print(adata.obs["cell_type"].value_counts().sort_index().to_string(), "\n")

# --- the pipeline ----------------------------------------------------------
fit = pace.fit_anndata(adata, "cell_type", n_threads=1)
print(fit, "\n")

slopes = pace.pace_shrink(fit, null_correlation=True)   # R's default
print(slopes)
print("  ", slopes.calls(), "\n")

stats = pace.fit_statistics(fit, counts)
blocks = pace.pace_decompose(fit, counts, stats=stats)
frame = pace.single_frame_decomposition(fit, counts, blocks)
drivers = pace.pace_drivers(fit, slopes, blocks, stats)
pairs = pace.pair_variance(fit, counts, frame)
print(frame, "\n", drivers, "\n")

print("variance by focal type (observed, %):")
print(frame.by_focal().round(2).to_string(index=False), "\n")

strongest = np.argsort(pairs["val"])[::-1][:5]
print("strongest spatial pairs:")
for i in strongest:
    print(f"  {pairs['focal'][i]:<16s} <- {pairs['neighbour'][i]:<16s} {pairs['val'][i]:5.2f}%")

print("\nmacrophage drivers near tumour:")
print(drivers["Macrophage", "Tumour"].to_dataframe()
      .loc[:, ["rank", "gene", "driver_score", "b_clean", "spec", "lfsr"]]
      .round(4).to_string(index=False))

# --- the figures -----------------------------------------------------------
# Each panel on its own, and all four as one multi-page PDF. PdfPages keeps
# every page vector rather than stitching rasters together.
from matplotlib.backends.backend_pdf import PdfPages

figures = []

full, zoom = pace.plot_decomposition(frame)
figures.append((full.figure, "py_decomposition.pdf"))

axis = pace.plot_pair_heatmap(pairs, cell_types=fit.cell_types)
figures.append((axis.figure, "py_pairheatmap.pdf"))

left, right = pace.plot_drivers(drivers, "Macrophage", "Tumour", single_frame=frame)
figures.append((left.figure, "py_drivers.pdf"))

panels = pace.plot_proximity(fit, counts, ["MRC1", "APOC1"], "Macrophage", "Tumour")
figures.append((panels[0].figure, "py_proximity.pdf"))

for figure, name in figures:
    figure.savefig(OUT / name, bbox_inches="tight")

combined = OUT / "pace_python_figures.pdf"
with PdfPages(combined) as pdf:
    pdf.infodict().update(Title="PACE for Python: breast cancer Xenium crop",
                          Subject="Figures from examples/demo.py")
    for figure, _ in figures:
        pdf.savefig(figure, bbox_inches="tight")

print("\nwrote", ", ".join(name for _, name in figures))
print("combined into", combined.name, f"({len(figures)} pages)")
