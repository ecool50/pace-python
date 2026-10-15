"""Constructing and inspecting a fit, without running one.

The neighbourhood and the contamination field can be built on their own, which
is useful for checking what the model is being given. `stash_pace` and
`compact_fit` are about carrying a fit around once it exists.
"""
from __future__ import annotations

import numpy as np

from . import _core
from .fit import _anndata_coordinates, _anndata_column, _csc


def build_neighbourhood(data, cell_type=None, *, x=None, y=None, image=None,
                        h_bio=30.0, h_tech=5.0, eps=None, spatial_key="spatial",
                        n_threads=1):
    """The kernel-weighted neighbour abundances, before any fitting.

    Returns the Gaussian biological kernel that carries the proximity
    coefficients and the short-range exponential technical kernel behind the
    contamination field, each cells by cell types.
    """
    if hasattr(data, "obsm"):
        x, y = _anndata_coordinates(data, spatial_key)
        labels = _anndata_column(data, cell_type, "cell types")
    else:
        if x is None or y is None:
            raise ValueError("pass x and y, or an AnnData with coordinates in obsm")
        labels = np.asarray(cell_type if cell_type is not None else data, dtype=object)
    labels = np.asarray([str(v) for v in np.ravel(labels)], dtype=object)
    types = sorted(set(labels))
    code = np.ascontiguousarray([types.index(v) for v in labels], dtype=np.int32)
    group = (np.zeros(code.size, dtype=np.int32) if image is None else
             np.ascontiguousarray(
                 [sorted(set(map(str, image))).index(str(v)) for v in image],
                 dtype=np.int32))
    kernels = _core.neighbour_kernels(
        np.asfortranarray(np.asarray(x, dtype=float).ravel()),
        np.asfortranarray(np.asarray(y, dtype=float).ravel()),
        code, len(types), group, per_group=False, h_bio=h_bio, h_tech=h_tech,
        eps=3.0 * h_bio if eps is None else eps, n_threads=n_threads)
    return {"k_bio": kernels["k_bio"], "k_tech": kernels["k_tech"], "cell_types": types}


def ambient_field(data, cell_type=None, *, x=None, y=None, image=None,
                  h_tech=5.0, spatial_key="spatial", n_threads=1):
    """The sparse cross-cell-type contamination field, in compressed-column form.

    Only cells of a different type contribute, which is what makes the
    contamination separable from the biology at all.
    """
    if hasattr(data, "obsm"):
        x, y = _anndata_coordinates(data, spatial_key)
        labels = _anndata_column(data, cell_type, "cell types")
    else:
        labels = np.asarray(cell_type if cell_type is not None else data, dtype=object)
    labels = np.asarray([str(v) for v in np.ravel(labels)], dtype=object)
    types = sorted(set(labels))
    code = np.ascontiguousarray([types.index(v) for v in labels], dtype=np.int32)
    if image is None:
        group, n_images = np.zeros(code.size, dtype=np.int32), 1
    else:
        levels = sorted({str(v) for v in image})
        group = np.ascontiguousarray([levels.index(str(v)) for v in image], dtype=np.int32)
        n_images = len(levels)
    return _core.ambient_field(
        np.asfortranarray(np.asarray(x, dtype=float).ravel()),
        np.asfortranarray(np.asarray(y, dtype=float).ravel()),
        code, group, n_images=n_images, h_tech=h_tech, edge_correct=True,
        n_threads=n_threads)


def anchor_genes(fit, counts, top_n=10):
    """The genes that most mark each cell type, which anchor its contamination.

    The contamination term is identified by genes a type expresses and its
    neighbours do not, so these are worth looking at when a fit's contamination
    estimates seem wrong: a type with no distinctive genes has nothing to
    anchor it.
    """
    from .stats import fit_statistics
    stats = fit._statistics or fit_statistics(fit, counts)
    means = stats.mu_mean
    genes = np.asarray(fit.genes, dtype=object)
    out = {}
    for index, name in enumerate(fit.cell_types):
        others = np.delete(means, index, axis=0)
        # how much higher this type's mean is than the loudest other type
        margin = means[index, :] - others.max(axis=0)
        order = np.argsort(-margin)[:top_n]
        out[name] = [(genes[g], float(margin[g])) for g in order]
    return out


def stash_pace(adata, fit, key="PACE"):
    """Copy a fit's readouts into `adata.uns`, leaving the solver state behind.

    A fit on a real cohort is large and usually outlives the session that made
    it. This keeps the tables with the data instead.
    """
    if not hasattr(adata, "uns"):
        raise TypeError("stash_pace needs an AnnData")
    stash = {"cell_types": list(fit.cell_types), "params": dict(fit.params)}
    if fit.neighbour_slopes is not None:
        stash["neighbour_slopes"] = fit.neighbour_slopes.to_dataframe()
    if fit.single_frame is not None:
        stash["variance_decomposition"] = fit.single_frame.to_dataframe()
    if fit.top_drivers is not None:
        stash["top_drivers"] = fit.top_drivers.to_dataframe()
    if fit.contam_fraction is not None:
        adata.obs[f"{key}_contam_fraction"] = fit.contam_fraction
    adata.uns[key] = stash
    return adata


def compact_fit(fit):
    """Drop what a saved fit does not need, keeping every readout.

    The designs and the kernel are the bulk of a fit and are rebuildable; the
    tables are not. A compacted fit can still be reported from, but not
    re-decomposed.
    """
    import copy
    small = copy.copy(fit)
    for field in ("design", "ambient_weights", "kernel", "fixed", "_counts",
                  "_statistics", "x", "y"):
        setattr(small, field, None)
    return small
