#!/usr/bin/env python
""" A set of startup parameters for running notebooks with.

I am constantly running the same first few things in cells. To simplify things
I thought I would just have a startup script that gets run as the first cell in
each notebook.
"""

# Activate the autoreload extension for easy reloading of external packages
%reload_ext autoreload
autoreload 1

# Trun on the water mark
%reload_ext watermark
%watermark -u -d -v
