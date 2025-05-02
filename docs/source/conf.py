# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
from sphinx_gallery.sorting import ExplicitOrder
# from sphinx_gallery.sorting import ExampleTitleSortKey 
from sphinx_gallery.sorting import FileNameSortKey
pykpproot = os.path.abspath('../../')
sys.path.insert(0, pykpproot)


with open('../../pykpp/__init__.py', 'r') as initf:
    for _l in initf.readlines():
        if _l.startswith('__version__ = '):
            release = _l.split(' = ')[-1][1:-1]
            break
    else:
        release = '0.0.0'

# -- Project information -----------------------------------------------------

project = 'pykpp'
copyright = '2023, Barron H. Henderson'
author = 'Barron H. Henderson'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.githubpages',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'IPython.sphinxext.ipython_directive',
    'IPython.sphinxext.ipython_console_highlighting',
    'matplotlib.sphinxext.plot_directive',
    'sphinx_copybutton',
    'sphinx_design',
    #'sphinx_rtd_theme',
    'myst_nb',
    'sphinx_gallery.gen_gallery',
    'sphinx.ext.napoleon',
]

sphinx_gallery_conf = {
    'examples_dirs': '../../examples',
    'gallery_dirs': 'auto_examples',
    'subsection_order': ExplicitOrder([
        '../../examples',
        '../../examples/knote_ae_2015',
    ]),
    'within_subsection_order': FileNameSortKey,
}

# Generate the API documentation when building
autoclass_content = 'both'
autosummary_generate = True
autosummary_imported_members = True

html_sidebars = {
    'userguide': ['searchbox.html', 'sidebar-nav-bs.html'],
    'API': ['searchbox.html', 'sidebar-nav-bs.html'],
    'examples': ['searchbox.html', 'sidebar-nav-bs.html'],
    'notebook-gallery': ['searchbox.html', 'sidebar-nav-bs.html'],
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', '**.ipynb_checkpoints']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
