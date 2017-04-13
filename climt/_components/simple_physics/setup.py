from setuptools import setup, Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system
import os

fc = os.getenv('FC', 'gfortran ')
fflags = os.getenv('FFLAGS', ' -fPIC -fno-range-check ')
f_opt_flags = os.getenv('CLIMT_OPTIMIZE_FLAG', '-O3')
f_no_opt_flags = os.getenv('CLIMT_NO_OPTIMIZE_FLAG', ' -O0 ')
ldflags = os.getenv('LDFLAGS', '-lgfortran')
cflags = os.getenv('CFLAGS', '-fPIC')

# compile the fortran modules without linking
fortran_mod_comp = fc+' -Wall simple_physics_custom.f90 -c -o simple_physics_custom.o'+fflags+' '+f_opt_flags
print(fortran_mod_comp)
system(fortran_mod_comp)

ext_modules = [
    Extension(  # module name:
        '_simple_physics',
        # source file:
        ['_simple_physics.pyx'],
        # other compile args for gcc
        extra_compile_args=[cflags, f_opt_flags, ldflags],
        # other files to link to
        extra_link_args=['simple_physics_custom.o', ldflags])]

setup(name='_simple_physics',
      cmdclass={'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs=[get_include()],
      ext_modules=ext_modules)
