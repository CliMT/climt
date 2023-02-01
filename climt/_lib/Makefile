FFTW_DIR = fftw-3.3.6-pl2
SHT_DIR = shtns-2.7-r557
BLAS_DIR = OpenBLAS-0.2.19
SIM_PHYS_DIR = simple_physics
EMAN_DIR = emanuel
RRTMG_LW_DIR = rrtmg_lw
RRTMG_SW_DIR = rrtmg_sw
DCMIP_DIR=dcmip

FFTW_GZ = fftw-3.3.6-pl2.tar.gz
SHT_GZ = shtns-2.7-r557.tar.gz
BLAS_GZ = OpenBLAS-0.2.19.tar.gz

BASE_DIR = $(PWD)

CLIMT_FFTW_FLAGS = --enable-openmp --enable-avx --enable-avx2

LIB_DIR= $(BASE_DIR)/$(CLIMT_ARCH)/
CLIMT_CFLAGS =  "-I$(BASE_DIR) -L$(LIB_DIR) $(CFLAGS)"
LDFLAGS = "-L$(LIB_DIR)"


ifeq ($(CLIMT_ARCH),Windows)
	BASE_DIR = $(subst \,/,$(PWD))
	LIB_DIR= $(BASE_DIR)/$(CLIMT_ARCH)/
	CLIMT_FFTW_FLAGS = --with-our-malloc --build=x86_64-w64-mingw32 --host=x86_64-w64-mingw32 --enable-static --with-g77-wrappers --with-windows-f77-mangling --enable-sse2 --enable-avx
endif


all: simple_physics_lib emanuel_lib rrtmg_lw_lib rrtmg_sw_lib libdcmip

# Simple Physics Configuration

simple_physics_lib: $(CLIMT_ARCH)/libsimple_physics.a

$(CLIMT_ARCH)/libsimple_physics.a: $(SIM_PHYS_DIR)/*.f90
	cd $(SIM_PHYS_DIR); make; cp libsimple_physics.a $(LIB_DIR); touch ../../_components/simple_physics/_simple_physics.pyx

# Emanuel Convection Configuration

emanuel_lib: $(CLIMT_ARCH)/libemanuel.a

$(CLIMT_ARCH)/libemanuel.a: $(EMAN_DIR)/*.f90
	cd $(EMAN_DIR); make; cp libemanuel.a $(LIB_DIR); touch ../../_components/emanuel/_emanuel_convection.pyx

# RRTMG LW Configuration

rrtmg_lw_lib: $(CLIMT_ARCH)/librrtmg_lw.a

$(CLIMT_ARCH)/librrtmg_lw.a: $(RRTMG_LW_DIR)/*.f90
	cd $(RRTMG_LW_DIR); make; cp librrtmg_lw.a $(LIB_DIR); touch ../../_components/rrtmg/lw/_rrtmg_lw.pyx

# RRTMG SW Configuration

rrtmg_sw_lib: $(CLIMT_ARCH)/librrtmg_sw.a

$(CLIMT_ARCH)/librrtmg_sw.a: $(RRTMG_SW_DIR)/*.f90
	cd $(RRTMG_SW_DIR); make; cp librrtmg_sw.a $(LIB_DIR); touch ../../_components/rrtmg/sw/_rrtmg_sw.pyx

# DCMIP Configuration

libdcmip: $(CLIMT_ARCH)/libdcmip.a

$(CLIMT_ARCH)/libdcmip.a: $(DCMIP_DIR)/*.f90
	cd $(DCMIP_DIR); make; cp libdcmip.a $(LIB_DIR); touch ../../_components/dcmip/_dcmip.pyx



clean:
	if [ -d $(FFTW_DIR) ]; then rm -Rf $(FFTW_DIR); fi
	if [ -d $(SHT_DIR) ]; then rm -Rf $(SHT_DIR); fi
	if [ -d $(BLAS_DIR) ]; then rm -Rf $(BLAS_DIR); fi
	cd $(SIM_PHYS_DIR); make clean
	cd $(EMAN_DIR); make clean
	cd $(RRTMG_LW_DIR); make clean
	cd $(RRTMG_SW_DIR); make clean
	cd $(DCMIP_DIR); make clean

clean_libs:
	if [ -d  include ]; then rm -Rf include; fi
	if [ -d  share ]; then rm -Rf share; fi
	if [ -d  bin ]; then rm -Rf bin; fi
	if [ -d  lib ]; then rm -Rf lib; fi
	rm $(LIB_DIR)/*.a
	rm $(LIB_DIR)/.conf*
