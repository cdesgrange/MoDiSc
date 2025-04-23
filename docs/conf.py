# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'MoDiSc'
copyright = '2025, Célia Desgrange'
author = 'Célia Desgrange'

release = '1.0'
version = '1.1.9'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'

# -- Options for logo

html_logo = "_static/logo.jpeg"

# -- Ensure that the folder _static is taken into account
html_static_path = ["_static"]

master_doc = 'index'

# Color
from docutils import nodes
from docutils.parsers.rst import roles

def red_role(name, rawtext, text, lineno, inliner, options={}, content=[]):
    node = nodes.raw('', f'<span style="color: #D52600;">{text}</span>', format='html')
    return [node], []

roles.register_local_role('red', red_role)
