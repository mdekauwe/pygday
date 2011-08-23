#!/usr/bin/env python
"""
Build the G'Day model

that's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (09.02.2011)"
__email__ = "mdekauwe@gmail.com"

from distutils.core import setup

setup(name="pygday",
    version="1.0",
    description="A Python implementation of the G'DAY (Generic Decomposition And Yield) model",
    long_description="GDAY model, simulates the cycling of C, N and water in plants and soil",
    author="Martin De Kauwe",
    author_email='mdekauwe@gmail.com',
    platforms = ['any'],
    package_dir = {'': 'src'},
    packages = ['gday']
)