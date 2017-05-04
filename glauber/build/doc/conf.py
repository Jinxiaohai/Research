# -*- coding: utf-8 -*-

import sphinx_rtd_theme

project = 'trento'
version = release = '1.4.0'
copyright = '2015 Jonah E. Bernhard, J. Scott Moreland, Steffen A. Bass'

extensions = ['breathe', 'sphinx.ext.mathjax']
breathe_projects = {project: '/home/xiaohai/Github/Research/glauber/build/doc/doxyxml'}
breathe_default_project = project
breathe_implementation_filename_extensions = ['.cxx']

source_suffix = '.rst'
master_doc = 'index'
templates_path = ['/home/xiaohai/Github/Research/glauber/doc/_templates']

primary_domain = 'cpp'
highlight_language = 'python'
default_role = 'math'
pygments_style = 'sphinx'

html_domain_indices = False
html_use_index = False
html_static_path = ['/home/xiaohai/Github/Research/glauber/doc/_static']
html_favicon = '/home/xiaohai/Github/Research/glauber/doc/_static/favicon.png'

html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_context = dict(
    display_github=True,
    github_user='Duke-QCD',
    github_repo='trento',
    github_version='master',
    conf_py_path='/doc/',
    source_suffix = source_suffix,
)
