#!/usr/bin/env python
# coding=utf-8

from setuptools import setup, find_packages
from codecs import open
import sys
import os

here = os.path.abspath(os.path.dirname(__file__))
sys.path.append(here)
import versioneer  # noqa: E402


CLASSIFIERS = """
Development Status :: 1 - Beta
Intended Audience :: Science/Research
License :: OSI Approved :: IMT-Atlantique License
Programming Language :: Python :: 3.12.12
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

INSTALL_REQUIRES = [
    "veros>=1.6"
]

with open(os.path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()


setup(
    name="veros_global_realistic",
    license="IMT-Atlantique",
    author="Elliot Chevet",
    author_email="elliot.chevet@imt-atlantique.fr",
    keywords="oceanography python geophysics ocean-model veros plugin",
    description="Global realistic setup/plugin for Veros",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/elliotchevet/veros_global_realistic",
    python_requires=">=3.12.12",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    install_requires=INSTALL_REQUIRES,
    entry_points={
        "veros.setup_dirs": [
            "global_realistic = veros_global_realistic.setup"
        ]
    },
    package_data={"veros_global_realistic": ["setup/*/assets.json"]},
    classifiers=[c for c in CLASSIFIERS.split("\n") if c],
)
