from setuptools import setup, Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system
import os

# compile the fortran modules without linking
module_list = [
    'parkind.f90',
    'parrrsw.f90',
    'rrsw_cld.f90',
    'rrsw_con.f90',
    'rrsw_kg16.f90',
    'rrsw_kg17.f90',
    'rrsw_kg18.f90',
    'rrsw_kg19.f90',
    'rrsw_kg20.f90',
    'rrsw_kg21.f90',
    'rrsw_kg22.f90',
    'rrsw_kg23.f90',
    'rrsw_kg24.f90',
    'rrsw_kg25.f90',
    'rrsw_kg26.f90',
    'rrsw_kg27.f90',
    'rrsw_kg28.f90',
    'rrsw_kg29.f90',
    'rrsw_ncpar.f90',
    'rrsw_ref.f90',
    'rrsw_tbl.f90',
    'rrsw_vsn.f90',
    'rrsw_aer.f90',
    'rrsw_wvn.f90']

sources_list = [
    'rrtmg_sw_cldprop.f90',
    'rrtmg_sw_cldprmc.f90',
    'rrtmg_sw_taumol.f90',
    'rrtmg_sw_vrtqdr.f90',
    'rrtmg_sw_reftra.f90',
    'rrtmg_sw_spcvmc.f90',
    'rrtmg_sw_setcoef.f90',
    'rrtmg_sw_spcvrt.f90',
    'rrtmg_sw_rad.nomcica.f90',
    'mcica_random_numbers.f90',
    'rrtmg_sw_init.f90',
    'mcica_subcol_gen_sw.f90',
    'rrtmg_sw_rad.f90',
    'rrtmg_sw_c_binder.f90']

unoptimised_sources_list = [
    'rrtmg_sw_k_g.f90',
]
object_file_list = []

fc = os.getenv('FC', 'gfortran ')
fflags = os.getenv('FFLAGS', ' -fPIC -fno-range-check ')
f_opt_flags = os.getenv('CLIMT_OPTIMIZE_FLAG', '-O3')
f_no_opt_flags = os.getenv('CLIMT_NO_OPTIMIZE_FLAG', ' -O0 ')
ldflags = os.getenv('LDFLAGS', '-lgfortran')
cflags = os.getenv('CFLAGS', '-fPIC')

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

    compilation_command = fc+source+' -c -o '+output_file+' '+f_no_opt_flags+fflags
    print(compilation_command)
    system(compilation_command)

link_args_list = object_file_list + [ldflags]

ext_modules = [
    Extension(  # module name:
        '_rrtm_sw',
        # source file:
        ['_rrtm_sw.pyx'],
        # other compile args for gcc
        extra_compile_args=[cflags, f_opt_flags, ldflags],
        # other files to link to
        extra_link_args=link_args_list)]

setup(name='_rrtm_sw',
      cmdclass={'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs=[get_include()],
      ext_modules=ext_modules)
