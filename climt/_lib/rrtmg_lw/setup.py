from setuptools import setup, Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system
import os

# compile the fortran modules without linking
module_list = [
    'parkind.f90',
    'parrrtm.f90',
    'rrlw_cld.f90',
    'rrlw_con.f90',
    'rrlw_kg01.f90',
    'rrlw_kg02.f90',
    'rrlw_kg03.f90',
    'rrlw_kg04.f90',
    'rrlw_kg05.f90',
    'rrlw_kg06.f90',
    'rrlw_kg07.f90',
    'rrlw_kg08.f90',
    'rrlw_kg09.f90',
    'rrlw_kg10.f90',
    'rrlw_kg11.f90',
    'rrlw_kg12.f90',
    'rrlw_kg13.f90',
    'rrlw_kg14.f90',
    'rrlw_kg15.f90',
    'rrlw_kg16.f90',
    'rrlw_ncpar.f90',
    'rrlw_ref.f90',
    'rrlw_tbl.f90',
    'rrlw_vsn.f90',
    'rrlw_wvn.f90']

sources_list = [
    'rrtmg_lw_cldprop.f90',
    'rrtmg_lw_cldprmc.f90',
    'rrtmg_lw_rtrn.f90',
    'rrtmg_lw_rtrnmr.f90',
    'rrtmg_lw_rtrnmc.f90',
    'rrtmg_lw_setcoef.f90',
    'rrtmg_lw_taumol.f90',
    'rrtmg_lw_rad.nomcica.f90',
    'mcica_random_numbers.f90',
    'rrtmg_lw_init.f90',
    'mcica_subcol_gen_lw.f90',
    'rrtmg_lw_rad.f90',
    'rrtmg_lw_c_binder.f90']

unoptimised_sources_list = [
    'rrtmg_lw_k_g.f90',
]
object_file_list = []

fc = os.getenv('FC', 'gfortran ')
fflags = os.getenv('FFLAGS', ' -fPIC -fno-range-check ')
cflags = os.getenv('CFLAGS', '-fPIC')
f_opt_flags = os.getenv('CLIMT_OPTIMIZE_FLAG', '-O3')
f_no_opt_flags = os.getenv('CLIMT_NO_OPTIMIZE_FLAG', ' -O0 ')
ldflags = os.getenv('LDFLAGS', '-lgfortran')

print('Compiling Modules')
for module in module_list:

    output_file = module[:-3]+'o'
    object_file_list.append(output_file)

    compilation_command = fc+module+' -c -o '+output_file+' '+f_opt_flags+fflags
    print(compilation_command)
    system(compilation_command)

print('Compiling Sources')
for source in sources_list:

    output_file = source[:-3]+'o'
    object_file_list.append(output_file)

    compilation_command = fc+source+' -c -o '+output_file+' '+f_opt_flags+fflags
    print(compilation_command)
    system(compilation_command)

print('Compiling k coefficient tables')
for source in unoptimised_sources_list:

    output_file = source[:-3]+'o'
    object_file_list.append(output_file)

    compilation_command = fc+source+' -c -o '+output_file+f_no_opt_flags+fflags
    print(compilation_command)
    system(compilation_command)

link_args_list = object_file_list + [ldflags]

ext_modules = [
    Extension(  # module name:
        '_rrtm_lw',
        # source file:
        ['_rrtm_lw.pyx'],
        # other compile args for gcc
        extra_compile_args=[cflags, f_opt_flags, ldflags],
        # other files to link to
        extra_link_args=link_args_list)]

setup(name='_rrtm_lw',
      cmdclass={'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs=[get_include()],
      ext_modules=ext_modules)
