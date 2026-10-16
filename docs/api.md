# Reference

Everything the package exports, built from the docstrings.

## The pipeline

`pace_fit` runs all four stages and returns one object. The stages each take a
fit and return it with one more layer attached, for when you want to stop
partway or re-run a later one with different settings.

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   pace.pace_fit
   pace.pace_model
   pace.pace_shrink
   pace.pace_decompose
   pace.pace_drivers
```

## Readouts

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   pace.neighbour_slopes
   pace.variance_decomposition
   pace.top_drivers
   pace.pair_variance
   pace.cell_contamination
   pace.anchor_genes
   pace.cell_types
   pace.pace_params
   pace.fit_state
   pace.fit_statistics
   pace.fitted_means
```

## Figures

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   pace.plot_decomposition
   pace.plot_pair_heatmap
   pace.plot_drivers
   pace.plot_proximity
   pace.plot_shrinkage
```

## Building and inspecting a fit

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   pace.build_neighbourhood
   pace.ambient_field
   pace.stash_pace
   pace.compact_fit
```

## Objects

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   pace.PaceFit
   pace.FitStatistics
```

## The numerical core

Single stages, for anyone who wants one on its own rather than a whole fit.

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   pace.mash
   pace.driver_scores
   pace.pair_variance_pratt
```
