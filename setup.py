#!/usr/bin/env python
# coding=utf-8

from setuptools import setup, find_packages
import sys
import os

here = os.path.abspath(os.path.dirname(__file__))
sys.path.append(here)


CLASSIFIERS = """
Development Status :: 1 - Beta
Intended Audience :: Science/Research
Programming Language :: Python :: 3.12.12
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

INSTALL_REQUIRES = [
    "veros[jax]>=1.6.2",
    "xarray",
    "h5netcdf",
    "dask",
]

EXTRAS_REQUIRE = {
    "tools": [
        "gcm_filters",
        "xesmf",
        "netCDF4",
        "matplotlib",
        "dask[distributed]",
        "cdsapi>=0.7.7",
        "copernicusmarine",
    ],
}

with open(os.path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()


setup(
    name="veros_global_realistic",
    license="MIT",
    author="Elliot Chevet",
    author_email="elliot.chevet@imt-atlantique.fr",
    keywords="oceanography python geophysics ocean-model veros plugin",
    description="Global realistic setup/plugin for Veros",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/elliotchevet/veros_global_realistic",
    python_requires=">=3.12.12",
    packages=find_packages(),
    version="0.1.0",
    install_requires=INSTALL_REQUIRES,
    extras_require=EXTRAS_REQUIRE,
    entry_points={
        "veros.setup_dirs": [
            "global_realistic = veros_global_realistic.setup"
        ]
    },
    package_data={"veros_global_realistic": ["setup/*/assets.json"]},
    classifiers=[c for c in CLASSIFIERS.split("\n") if c],
)
