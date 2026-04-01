# Configuration file for the Sphinx documentation builder.

project = "PyGDSM"
author = "Danny C. Price"
release = "1.6.4"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "nbsphinx",
    "myst_parser",
]

# Allow .md files in the toctree
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

html_theme = "sphinx_rtd_theme"

# Don't re-execute notebooks at build time (data files may not be available)
nbsphinx_execute = "never"

exclude_patterns = ["_build", "**.ipynb_checkpoints"]
