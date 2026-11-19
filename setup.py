#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import find_packages, setup

with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

requirements = [
    "numpy>=1.16.0",
    "pint>=0.7.0",
    "unyt",
    "xarray>=0.8.0",
    "sympl>=0.5.0",
    "importlib_resources",
]

test_requirements = [
    "pytest>=2.9.2",
    "mock>=2.0.0",
]

setup(
    name="climt",
    version="0.31.0",
    description="CliMT is a Toolkit for building Earth system models in Python.",
    long_description=readme + "\n\n" + history,
    author="Rodrigo Caballero",
    author_email="rodrigo.caballero@misu.su.se",
    url="https://github.com/CliMT/climt",
    packages=find_packages(include=["climt", "climt.*"]),
    package_dir={"climt": "climt"},
    package_data={
        "climt": [
            "_data/ozone_profile.npy",
            "_data/atmospheric_properties/*.toml",
            "_data/cork/correlated_k/*.nc",
            "_data/cork/correlated_k/*.npz",
            "_data/land_mask/*.nc",
            "_data/data_ocean/*.nc",
            "_data/topography/*.nc",
        ]
    },
    install_requires=requirements,
    license="BSD license",
    zip_safe=False,
    keywords="climt",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
    ],
    test_suite="tests",
    tests_require=test_requirements,
)
