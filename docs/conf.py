"""Sphinx configuration: the Python counterpart to the R package's pkgdown site.

Same shape as that site -- an introduction, a reference built from the
docstrings, and the worked example -- so the two read alike.
"""
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

project = "PACE"
author = "Elijah Willie"
copyright = f"{date.today().year}, Elijah Willie"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",      # the docstrings are numpydoc-flavoured
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "myst_parser",
]

autosummary_generate = True
autodoc_member_order = "bysource"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
myst_enable_extensions = ["colon_fence", "deflist"]
# the README and usage page are written for GitHub, where headings are not
# required to nest; don't fail the build over that
suppress_warnings = ["myst.header"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
    "anndata": ("https://anndata.readthedocs.io/en/stable", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
}

templates_path = ["_templates"]
exclude_patterns = ["_build"]
html_theme = "pydata_sphinx_theme"
html_static_path = ["figures"]
html_title = "PACE"

html_theme_options = {
    "github_url": "https://github.com/ecool50/pace-python",
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "show_prev_next": False,
    "icon_links": [
        {"name": "The R package", "url": "https://github.com/ecool50/PACE",
         "icon": "fa-brands fa-r-project"},
    ],
}
