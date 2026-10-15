"""PACE: proximity-associated changes in expression.

The same shape as the R package. `pace_fit()` runs everything and hands back
one object; `pace_model()`, `pace_shrink()`, `pace_decompose()` and
`pace_drivers()` are the four stages if you want them separately. The readouts
are attributes of the fit and also functions over it.
"""

from .fit import PaceFit, pace_model
from .pipeline import (
    cell_contamination,
    cell_types,
    fit_state,
    neighbour_slopes,
    pace_decompose,
    pace_drivers,
    pace_fit,
    pace_params,
    pace_shrink,
    pair_variance,
    top_drivers,
    variance_decomposition,
)
from .primitives import ambient_field, anchor_genes, build_neighbourhood, compact_fit, stash_pace
from .stats import FitStatistics, fit_statistics, fitted_means
from .plots import (plot_decomposition, plot_drivers, plot_pair_heatmap,
                    plot_proximity, plot_shrinkage)
from ._core import mash, driver_scores, pair_variance_pratt

__all__ = [
    # the pipeline
    "pace_fit",
    "pace_model",
    "pace_shrink",
    "pace_decompose",
    "pace_drivers",
    # the readouts
    "neighbour_slopes",
    "variance_decomposition",
    "top_drivers",
    "pair_variance",
    "cell_contamination",
    "anchor_genes",
    "cell_types",
    "pace_params",
    "fit_state",
    "fit_statistics",
    "fitted_means",
    # the figures
    "plot_decomposition",
    "plot_pair_heatmap",
    "plot_drivers",
    "plot_proximity",
    "plot_shrinkage",
    # constructing and inspecting a fit
    "build_neighbourhood",
    "ambient_field",
    "stash_pace",
    "compact_fit",
    "PaceFit",
    "FitStatistics",
    # the numerical core, for anyone who wants a stage on its own
    "mash",
    "driver_scores",
    "pair_variance_pratt",
]
__version__ = "0.1.0"
