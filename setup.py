#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, Extension
from wheel.bdist_wheel import bdist_wheel as native_bdist_wheel
import os
import subprocess
import platform
import re
try:
    from pip import main as pip_main
except:
    from pip._internal import main as pip_main

try:
    from Cython.Build.Distutils import build_ext as native_build_ext
except ImportError:
    print('Suitable Cython unavailable, installing...')
    pip_main(['install', 'cython'])
    from Cython.Build.Distutils import build_ext as native_build_ext

try:
    import numpy as np
except ImportError:
    print('Suitable numpy unavailable, installing...')
    pip_main(['install', 'numpy'])
    import numpy as np


include_dirs = [np.get_include()]

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'numpy>=0.10',
    'pint>=0.7.0',
    'xarray>=0.8.0',
    'sympl==0.3.2',
    'cython>=0.25',
    'scipy>=0.18.1',
]

test_requirements = [
    'pytest>=2.9.2',
    'mock>=2.0.0',
]


# Platform specific settings
def guess_compiler_name(env_name):

    search_string = ''
    if env_name == 'FC':
        search_string = 'gfortran-\d$'
    if env_name == 'CC':
        search_string = 'gcc-\d$'

    for root, dirs, files in os.walk('/usr/local/Cellar/gcc/'):

        for line in files:
            if re.match(search_string, line):
                print('Using ', env_name, '= ', line)
                os.environ[env_name] = line


operating_system = platform.system()

libraries = ['m', 'gfortran']
default_link_args = []
default_compile_args = []

compiled_base_dir = 'climt/_lib'

if operating_system == 'Linux':
    libraries = ['m', 'gfortran', 'rt']
    default_link_args = ['-lgfortran']

if operating_system == 'Windows':
    compiled_base_dir = 'climt\\_lib'

dir_path = os.getcwd()
compiled_path = os.path.join(dir_path, compiled_base_dir)
lib_path = os.path.join(compiled_path, operating_system+'/')
inc_path = os.path.join(compiled_path, 'include')
include_dirs.append(inc_path)
openblas_path = lib_path+'/libopenblas.a'


# Compile libraries

if 'FC' not in os.environ:
    if operating_system == 'Darwin':
        guess_compiler_name('FC')
    else:
        os.environ['FC'] = 'gfortran'
        os.environ['F77'] = 'gfortran'

if 'CC' not in os.environ:
    if operating_system == 'Darwin':
        guess_compiler_name('CC')
    else:
        os.environ['CC'] = 'gcc'

if operating_system == 'Windows' and os.environ.get('APPVEYOR') == 'True':
    # os.environ['CC'] = 'x86_64-w64-mingw32-gcc.exe'
    # os.environ['FC'] = 'x86_64-w64-mingw32-gfortran.exe'
    # os.environ['AR'] = 'x86_64-w64-mingw32-gcc-ar.exe'
    libraries = []
    openblas_path = os.path.join(os.environ['COMPILER_PATH'], '../lib/libopenblas.a')
    default_link_args = ['-l:libgfortran.a', '-l:libquadmath.a', '-l:libm.a']
    default_compile_args = ['-DMS_WIN64']

os.environ['FFLAGS'] = '-fPIC -fno-range-check'
os.environ['CFLAGS'] = '-fPIC'

print('Compilers: ', os.environ['CC'], os.environ['FC'])


# Create a custom build class to build libraries, and patch cython extensions
def build_libraries():

    if os.environ.get('READTHEDOCS') == 'True':
        return

    curr_dir = os.getcwd()
    os.chdir(compiled_path)
    os.environ['PWD'] = compiled_path
    if subprocess.call(['make', 'CLIMT_ARCH='+operating_system]):
        raise RuntimeError('Library build failed, exiting')
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
        # patch_cython_binary()
        native_bdist_wheel.run(self)


# Define extensions to be built
if os.environ.get('READTHEDOCS') == 'True':
    ext_modules = []
else:
    ext_modules = [
        Extension(
            'climt._components._berger_solar_insolation',
            ['climt/_components/_berger_solar_insolation.pyx']),

        Extension(
            'climt._components.simple_physics._simple_physics',
            sources=['climt/_components/simple_physics/_simple_physics.pyx'],
            libraries=libraries,
            include_dirs=include_dirs,
            extra_compile_args=default_compile_args,
            library_dirs=[lib_path],
            extra_link_args=[lib_path+'/libsimple_physics.a'] + default_link_args),

        Extension(
            'climt._components.emanuel._emanuel_convection',
            sources=['climt/_components/emanuel/_emanuel_convection.pyx'],
            libraries=libraries,
            include_dirs=include_dirs,
            extra_compile_args=default_compile_args,
            library_dirs=[lib_path],
            extra_link_args=[lib_path+'/libemanuel.a'] + default_link_args),

        Extension(
            'climt._components.rrtmg.lw._rrtmg_lw',
            sources=['climt/_components/rrtmg/lw/_rrtmg_lw.pyx'],
            libraries=libraries,
            include_dirs=include_dirs,
            extra_compile_args=default_compile_args + ['-fopenmp'],
            library_dirs=[lib_path],
            extra_link_args=[lib_path+'/librrtmg_lw.a', '-fopenmp'] + default_link_args),

        Extension(
            'climt._components.gfs._gfs_dynamics',
            sources=['climt/_components/gfs/_gfs_dynamics.pyx'],
            libraries=libraries,
            include_dirs=include_dirs,
            extra_compile_args=['-fopenmp'] + default_compile_args,
            library_dirs=[lib_path],
            extra_link_args=['-fopenmp', lib_path+'/libgfs_dycore.a',
                             lib_path+'/libshtns_omp.a', lib_path+'/libfftw3_omp.a',
                             lib_path+'/libfftw3.a', openblas_path] + default_link_args),

        # lib_path+'/libshtns_omp.a', openblas_path, os.environ['COMPILER_PATH']+'../lib/libfftw3.a'] + default_link_args),

        Extension(
            'climt._components.rrtmg.sw._rrtmg_sw',
            sources=['climt/_components/rrtmg/sw/_rrtmg_sw.pyx'],
            libraries=libraries,
            include_dirs=include_dirs,
            extra_compile_args=default_compile_args + ['-fopenmp'],
            library_dirs=[lib_path],
            extra_link_args=[lib_path+'/librrtmg_sw.a', '-fopenmp'] + default_link_args),

        Extension(
            'climt._components.dcmip._dcmip',
            sources=['climt/_components/dcmip/_dcmip.pyx'],
            libraries=libraries,
            include_dirs=include_dirs,
            extra_compile_args=default_compile_args,
            library_dirs=[lib_path],
            extra_link_args=[lib_path+'/libdcmip.a'] + default_link_args),
    ]

setup(
    name='climt',
    version='0.14.7',
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
        'Programming Language :: Python :: 3.6',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
