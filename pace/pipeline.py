"""The staged pipeline and the readouts hung off a fit.

PACE runs in four stages. `pace_model()` fits; `pace_shrink()` stabilises the
slopes; `pace_decompose()` splits each gene's variance; `pace_drivers()` ranks
the genes carrying each pair. `pace_fit()` runs all four and hands back one
object with every layer attached, which is how you will usually want it.

The readouts are then attributes of that object, and also functions, so both
`fit.neighbour_slopes` and `neighbour_slopes(fit)` work.
"""
from __future__ import annotations

import numpy as np

from .decompose import pace_decompose as _decompose
from .decompose import pair_variance as _pair_variance
from .decompose import single_frame_decomposition as _single_frame
from .drivers import pace_drivers as _drivers
from .fit import pace_model
from .shrink import pace_shrink as _shrink
from .stats import fit_statistics


def pace_fit(counts, x=None, y=None, cell_type=None, *, pairs=None,
             resp_term=None, null_correlation=True, n_threads=1, verbose=False,
             **kwargs):
    """Fit PACE and attach every reporting layer.

    Runs `pace_model()`, `pace_shrink()`, `pace_decompose()` and
    `pace_drivers()` in turn, which is what you want unless you mean to stop
    partway or re-run a later stage with different settings.

    Takes the same inputs as `pace_model()`: a cells-by-genes count matrix
    with coordinates and labels, or an `AnnData` with `cell_type` naming the
    `obs` column.

    Parameters
    ----------
    pairs
        Which (focal, neighbour) pairs to rank. None does every ordered pair.
    resp_term
        For a condition fit, which condition term carries the responder layer.
    null_correlation
        Passed to `pace_shrink()`.

    Returns
    -------
    PaceFit, with `neighbour_slopes`, `variance_decomposition`,
    `single_frame`, `top_drivers` and `contam_fraction` filled in.
    """
    fit = pace_model(counts, x, y, cell_type, n_threads=n_threads,
                     verbose=verbose, **kwargs)
    # the model holds no counts, so the later stages take the matrix the model
    # was given, exactly as R hands them the SpatialExperiment again
    matrix = fit._counts if fit._counts is not None else counts
    return _complete(fit, matrix, pairs=pairs, resp_term=resp_term,
                     null_correlation=null_correlation, n_threads=n_threads,
                     verbose=verbose)


def _complete(fit, counts, *, pairs, resp_term, null_correlation, n_threads, verbose):
    fit = pace_shrink(fit, resp_term=resp_term,
                      null_correlation=null_correlation, verbose=verbose)
    fit = pace_decompose(fit, counts, resp_term=resp_term, n_threads=n_threads)
    return pace_drivers(fit, counts, pairs=pairs, resp_term=resp_term,
                        n_threads=n_threads)


# --- the stages, each taking a fit and returning it with one more layer ----

def pace_shrink(fit, *, resp_term=None, null_correlation=True, verbose=False):
    """Shrink the neighbour slopes and attach them to the fit."""
    fit.neighbour_slopes = _shrink(fit, resp_term=resp_term,
                                   null_correlation=null_correlation,
                                   verbose=verbose)
    return fit


def pace_decompose(fit, counts, *, resp_term=None, stats=None, n_threads=1):
    """Decompose the variance and attach both tables to the fit."""
    if stats is None:
        stats = fit_statistics(fit, counts)
    fit._statistics = stats
    fit.contam_fraction = stats.contam_fraction
    fit.variance_decomposition = _decompose(fit, counts, resp_term=resp_term,
                                            stats=stats, n_threads=n_threads)
    fit.single_frame = _single_frame(fit, counts, fit.variance_decomposition,
                                     n_threads=n_threads)
    return fit


def pace_drivers(fit, counts=None, *, pairs=None, resp_term=None, n_threads=1):
    """Rank each pair's drivers and attach the tables to the fit."""
    if fit.neighbour_slopes is None or fit.variance_decomposition is None:
        raise ValueError(
            "pace_drivers() needs the shrunken slopes and the decomposition; "
            "run pace_shrink() and pace_decompose() first, or use pace_fit()")
    stats = fit._statistics
    if stats is None:
        if counts is None:
            raise ValueError("pass the counts, or run pace_decompose() first")
        stats = fit_statistics(fit, counts)
        fit._statistics = stats
    fit.top_drivers = _drivers(fit, fit.neighbour_slopes, fit.variance_decomposition,
                               stats, pairs=pairs, resp_term=resp_term,
                               n_threads=n_threads)
    return fit


# --- the readouts, as functions over a fit --------------------------------

def _require(value, what, stage):
    if value is None:
        raise ValueError(f"this fit carries no {what}; run {stage} first")
    return value


def neighbour_slopes(fit):
    """The shrunken slopes, with a local false sign rate per effect."""
    return _require(fit.neighbour_slopes, "shrunken slopes", "pace_shrink()")


def variance_decomposition(fit, observed=True):
    """The per-gene variance decomposition.

    `observed=True` gives the single-frame table, whose magnitudes are
    measured from the data; False gives the fitted five-block one.
    """
    layer = fit.single_frame if observed else fit.variance_decomposition
    return _require(layer, "variance decomposition", "pace_decompose()")


def top_drivers(fit, focal=None, neighbour=None):
    """The per-pair driver tables, or one pair's."""
    tables = _require(fit.top_drivers, "driver tables", "pace_drivers()")
    if focal is None:
        return tables
    return tables[focal, neighbour]


def pair_variance(fit, block="spatial", resp_term=None, n_threads=1):
    """Each neighbour's share of a focal type's spatial variance."""
    frame = _require(fit.single_frame, "variance decomposition", "pace_decompose()")
    return _pair_variance(fit, None, frame, block=block, resp_term=resp_term,
                          n_threads=n_threads)


def cell_contamination(fit):
    """Per cell: the contamination loading, and the fraction it explains.

    Read the fraction rather than the loading. The ambient field is
    cross-cell-type, so a cell with no differently-typed neighbour inside the
    technical kernel has no ambient signal at all; its loading is unidentified
    and shrinks to the prior mean, giving every such cell the same apparently
    middling value. The fraction reports those cells as zero, correctly.
    """
    fraction = _require(fit.contam_fraction, "contamination fractions",
                        "pace_decompose()")
    try:
        import pandas as pd
    except ImportError:
        return {"cell": np.arange(fraction.size),
                "celltype": np.asarray(fit.cell_types, dtype=object)[fit.type_code],
                "rho": fit.rho, "contam_fraction": fraction}
    return pd.DataFrame({
        "cell": np.arange(fraction.size),
        "celltype": np.asarray(fit.cell_types, dtype=object)[fit.type_code],
        "rho": fit.rho,
        "contam_fraction": fraction})


def cell_types(fit):
    """The cell types, in the order they index the random-effect rows."""
    return list(fit.cell_types)


def pace_params(fit):
    """What the fit was run with."""
    return dict(fit.params)


def fit_state(fit):
    """The solver's own output: coefficients, errors, dispersions, loadings."""
    return {"B": fit.B, "U": fit.U, "se_B": fit.se_B, "se_U": fit.se_U,
            "alpha": fit.alpha, "tau": fit.tau, "rho": fit.rho,
            "row_names": fit.row_names, "genes": fit.genes,
            "n_iter": fit.n_iter, "converged": fit.converged}
