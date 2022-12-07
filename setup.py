#!/usr/bin/env python

import sys
# Test Python's version
major, minor = sys.version_info[0:2]
from setuptools import find_packages, setup
from distutils.command.sdist import sdist
cmdclass={'sdist': sdist}
#try:
#    from numpy.distutils.core import setup
#    from numpy.distutils.core import Extension
#except ImportError:
#    pass

def readme():
    with open('README.rst', encoding='utf-8', mode='r') as f:
        return f.read()

setup(name='lava2d',
    packages=find_packages(),
    scripts=[],
    version="1.0.0",
    description='Lava flow simulator',
    long_description=readme(),
    url = 'https://github.com/davemhyman/lava2d',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved',
        'License :: MIT',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
    ],
    author='David Hyman',
    author_email='',
    install_requires=[
        'numpy',
        'scipy',
        'numba',
        'pandas',
        'netcdf4',
        'sklearn',
        'matplotlib',
        'sphinx',
        'progressbar2',
        'numexpr',
        'geotiff',
        'shapely',
        'fiona',

    ],
)
