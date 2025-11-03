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
sys.path.insert(0, os.path.abspath('_scripts'))
import module_index_generator
import ersem_webpage
import subprocess
sys.path.insert(0, os.path.abspath('../../'))


# -- Project information -----------------------------------------------------

project = 'ERSEM'
copyright = '2021, The PML Modelling Group'
author = 'The PML Modelling Group'

# The full version, including alpha/beta/rc tags
release = '2021'

if os.path.isdir("ersem-setups"):
    os.chdir("ersem-setups")
    subprocess.run(['git', 'pull'])
    os.chdir("..")
else:
    subprocess.run(['git', 'clone', 'https://github.com/pmlmodelling/ersem-setups.git'])
# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['m2r2',
              'sphinx.ext.autosectionlabel',
              'sphinxcontrib.bibtex',
              'sphinx.ext.imgmath',
              'sphinx_panels']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

module_index_generator.create('module_index.rst')

on_rtd = os.environ.get('READTHEDOCS') == 'True'
if not on_rtd:
    import ersem_webpage
    ersem_webpage.generator_web_doc('model_info/ERSEM_model.rst')

autosectionlabel_prefix_document = True

bibtex_bibfiles = ['ref.bib']
