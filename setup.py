#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'numpy>=1.11',
]

test_requirements = [
    'pytest>=2.9.2',
    'mock>=2.0.0',
]

setup(
    name='climt',
    version='1.0.0',
    description="CliMT is a Toolkit for building Earth system models in Python.",
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
    license="BSD license",
    zip_safe=False,
    keywords='climt',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
