#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import sys

import setuptools
from setuptools import setup

sys.path.insert(0, "pyRosita")
from version import __version__


long_description = \
    """
A simple Python wrapper for using eRosita data.
Read the documentation at https://adina.feinste.in/pyrosita.

"""

with open('requirements.txt') as f:
    install_reqs = f.read().splitlines()


setup(
    name='pyRosita',
    version=__version__,
    license='MIT',
    author='Adina D. Feinstein',
    author_email='adina.d.feinstein@gmail.com',
    packages=[
        'pyRosita',
        ],
    include_package_data=True,
    url='http://github.com/afeinstein20/pyRosita',
    description='Python tools for eRosita Data',
    long_description=long_description,
    long_description_content_type="text/markdown",
    package_data={'': ['README.md', 'LICENSE']},
    install_requires=install_reqs,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.0',
        ],
    )
