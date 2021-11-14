from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import os
import setuptools

__version__ = '0.0.1'

setup(
    name='PyGrid2D',
    long_description='',
    install_requires=['numpy', 'scipy', 'argparse', 'matplotlib', 'pytest'],
    packages = ["pygrid2d"],
    package_dir = {"pygrid2d": "pygrid2d"},
    zip_safe=False,
)
