#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, Extension
try:
    import numpy as np
    include_dirs = [np.get_include()]
except ImportError:
    include_dirs = []

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'numpy>=1.10',
    'pint>=0.7.0',
    'xarray>=0.8.0',
    'matplotlib',
]

test_requirements = [
    'pytest>=2.9.2',
    'mock>=2.0.0',
]

ext_modules = [
    Extension(
        'climt._components._berger_solar_insolation',
        ['climt/_components/_berger_solar_insolation.pyx'])
]

setup(
    name='climt',
    version='1.0.0',
    description='CliMT is a Toolkit for building Earth system models in Python.',
    long_description=readme + '\n\n' + history,
    author="Rodrigo Caballero",
    author_email='rodrigo.caballero@misu.su.se',
    url='https://github.com/CliMT/climt',
    packages=[
        'climt',
    ],
    package_dir={'climt':
                 'climt'},
    include_package_data=True,
    install_requires=requirements,
    ext_modules=ext_modules,
    include_dirs=include_dirs,
    license="BSD license",
    zip_safe=False,
    keywords='climt',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
