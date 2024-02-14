"""Configuration file for the Sphinx documentation builder."""

import os
import sys

sys.path.insert(0, os.path.abspath("../../"))

# from rustyms import __version__

# Project information
project = "rustyms"
author = "Douwe Schulte, Ralf Gabriels"
github_project_url = "https://github.com/compomics/rustyms/"
github_doc_root = "https://github.com/compomics/rustyms/tree/release/docs/python/source/"

# Version
# release = __version__

# General configuration
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "myst_parser",
    "sphinx_rtd_theme",
]
source_suffix = [".rst", ".md"]
master_doc = "index"
exclude_patterns = ["_build"]

# Options for HTML output
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]

# Autodoc options
autodoc_default_options = {"members": True, "show-inheritance": True}
autodoc_member_order = "bysource"
autodoc_typehints = "both"
autoclass_content = "init"

# Intersphinx options
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}

def setup(app):
    config = {
        "enable_eval_rst": True,
    }
