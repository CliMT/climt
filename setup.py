#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, Extension
from Cython.Build.Distutils import build_ext as native_build_ext
from wheel.bdist_wheel import bdist_wheel as native_bdist_wheel

try:
    import numpy as np
    include_dirs = [np.get_include()]
except ImportError:
    include_dirs = []

import os
import subprocess
import platform

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'numpy>=1.10',
    'pint>=0.7.0',
    'xarray>=0.8.0',
    'matplotlib',
    'sympl',
]

test_requirements = [
    'pytest>=2.9.2',
    'mock>=2.0.0',
]


# Platform specific settings
operating_system = platform.system()

libraries = ['m', 'gfortran']

compiled_base_dir = 'climt/_lib/'

if operating_system == 'Linux':
    libraries.append('rt')

dir_path = os.getcwd()
compiled_path = os.path.join(dir_path, compiled_base_dir)
lib_path = os.path.join(compiled_path, operating_system+'/')
inc_path = os.path.join(compiled_path, 'include')
include_dirs.append(inc_path)

# Compile libraries

if 'FC' not in os.environ:
    os.environ['FC'] = 'gfortran'

if 'CC' not in os.environ:
    os.environ['CC'] = 'gcc'

os.environ['FFLAGS'] = '-fPIC -fno-range-check'
os.environ['CFLAGS'] = '-fPIC'


# Create a custom build class to build libraries, and patch cython extensions
def build_libraries():

    curr_dir = os.getcwd()
    os.chdir(compiled_path)
    os.environ['PWD'] = compiled_path
    subprocess.call(['make', 'CLIMT_ARCH='+operating_system])
    os.chdir(curr_dir)
    os.environ['PWD'] = curr_dir


def patch_cython_binary():
    if operating_system == 'Darwin':
        subprocess.call(['python', 'mac_os_patch.py'])
    return


# Custom build class
class climt_build_ext(native_build_ext):

    def run(self):
        build_libraries()
        native_build_ext.run(self)


# Custom bdist_wheel class
class climt_bdist_wheel(native_bdist_wheel):

    def run(self):
        self.run_command('build')
        patch_cython_binary()
        native_bdist_wheel.run(self)


# Define extensions to be built
ext_modules = [
    Extension(
        'climt._components._berger_solar_insolation',
        ['climt/_components/_berger_solar_insolation.pyx']),

    Extension(
        'climt._components.gfs._gfs_dynamics',
        sources=['climt/_components/gfs/_gfs_dynamics.pyx'],
        libraries=libraries,
        include_dirs=include_dirs,
        library_dirs=[lib_path],
        extra_link_args=['-fopenmp', lib_path+'/libgfs_dycore.a',
                         lib_path+'/libshtns_omp.a', lib_path+'/libfftw3_omp.a',
                         lib_path+'/libfftw3.a', lib_path+'/libopenblas.a']),

    Extension(
        'climt._components.simple_physics._simple_physics',
        sources=['climt/_components/simple_physics/_simple_physics.pyx'],
        libraries=libraries,
        include_dirs=include_dirs,
        library_dirs=[lib_path],
        extra_link_args=[lib_path+'/libsimple_physics.a']),

    Extension(
        'climt._components.emanuel._emanuel_convection',
        sources=['climt/_components/emanuel/_emanuel_convection.pyx'],
        libraries=libraries,
        include_dirs=include_dirs,
        library_dirs=[lib_path],
        extra_link_args=[lib_path+'/libemanuel.a']),

    Extension(
        'climt._components.rrtmg.lw._rrtmg_lw',
        sources=['climt/_components/rrtmg/lw/_rrtm_lw.pyx'],
        libraries=libraries,
        include_dirs=include_dirs,
        library_dirs=[lib_path],
        extra_link_args=[lib_path+'/librrtmg_lw.a']),

    Extension(
        'climt._components.rrtmg.sw._rrtmg_sw',
        sources=['climt/_components/rrtmg/sw/_rrtm_sw.pyx'],
        libraries=libraries,
        include_dirs=include_dirs,
        library_dirs=[lib_path],
        extra_link_args=[lib_path+'/librrtmg_sw.a']),

]

'''
fortran_ext = {
    'simple_physics': 'climt/_components/simple_physics',
    'rrtmg_longwave': 'climt/_components/rrtmg/lw',
    'rrtmg_shortwave': 'climt/_components/rrtmg/sw',
    'emanuel_convection': 'climt/_components/emanuel',
}

dir_path = os.getenv('PWD', '')

# os.environ['FC'] = 'gfortran'
# os.environ['FFLAGS'] = '-fPIC -fno-range-check'
# os.environ['CFLAGS'] = '-fPIC'
# os.environ['CLIMT_OPTIMIZE'] = ' -O3'
# os.environ['CLIMT_NO_OPTIMIZE'] = ' -O0'
# os.environ['LDFLAGS'] = '-lgfortran'

for module in fortran_ext.keys():
    mycwd = os.getcwd()
    os.chdir(os.path.join(dir_path, fortran_ext[module]))
    build_process = subprocess.call(['python', 'setup.py', 'build_ext', '--inplace'])
    os.chdir(mycwd)
'''

setup(
    name='climt',
    version='0.9',
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
    cmdclass={
        'build_ext': climt_build_ext,
        'bdist_wheel': climt_bdist_wheel,
    },
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
