#!/usr/bin/env python3

from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import glob
import os
import pkg_resources
# Note: the _program variable is set in __init__.py.
# it determines the name of the package/final command line tool.
from mink import __version__, _program


setup(name='mink',
      version=__version__,
      packages=find_packages(),
      scripts=["mink/scripts/mapping.py",
      "mink/scripts/mutation_funcs.py",
      "mink/scripts/report_writer.py"],
      package_data={"mink":["data/mapping_files/gadm36_GBR_2.json",
                             "data/mapping_files/NI_counties.geojson",
                             "data/mapping_files/channel_islands.json"]},
      install_requires=[
            "pytools>=2020.1",
            'pandas>=1.0.1',
            "matplotlib>=3.2.1",
            "scipy>=1.4.1",
            "numpy>=1.13.3",
            "geopandas>=0.7.0",
            "descartes>=1.1.0",
            "adjustText>=0.7.3",
            "grip>=4.5.2",
            "tabulate>=0.8.7",
            "seaborn>=0.11.0",
            "epiweeks>=2.1.1"
        ],
      description='Mutation I N K',
      url='https://github.com/COG-UK/mink',
      author='Verity Hill',
      author_email='verity.hill@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = mink.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
