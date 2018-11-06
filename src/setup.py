#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

setup(
    name='ncbi_remap',
    version='0.0.1',
    description="Local library for the NCBI Remap project",
    author="Justin M Fear",
    author_email='justin.m.fear@gmail.com',
    url='https://github.com/jfear/ncbi_remap',
    packages=['ncbi_remap'],
    license="MIT license",
    entry_points={
        'console_scripts':
        [
            'bigWigMerge = ncbi_remap.bigWigMerge:main',
        ],
    },
)
