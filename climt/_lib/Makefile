FFTW_DIR = fftw-3.3.6-pl2
SHT_DIR = shtns-2.7-r557
BLAS_DIR = OpenBLAS-0.2.19
GFS_DIR = GFS
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

CLIMT_LIB_TARGETS = blas_lib fftw_lib

ifeq ($(CLIMT_ARCH),Windows)
	BASE_DIR = $(subst \,/,$(PWD))
	LIB_DIR= $(BASE_DIR)/$(CLIMT_ARCH)/
	CLIMT_FFTW_FLAGS = --with-our-malloc --build=x86_64-w64-mingw32 --host=x86_64-w64-mingw32 --enable-static --with-g77-wrappers --with-windows-f77-mangling --enable-sse2 --enable-avx
	CLIMT_LIB_TARGETS =
endif


all: $(CLIMT_LIB_TARGETS) sht_lib gfs_lib simple_physics_lib emanuel_lib rrtmg_lw_lib rrtmg_sw_lib libdcmip

# FFTW3 Configuration

$(CLIMT_ARCH)/.configured_fftw:
	if [ ! -d $(FFTW_DIR) ]; then tar -xvzf $(FFTW_GZ) > log 2>&1; fi; cd $(FFTW_DIR); ./configure CFLAGS=$(CLIMT_CFLAGS) LDFLAGS=$(LDFLAGS) --prefix="$(BASE_DIR)" $(CLIMT_FFTW_FLAGS) > log.config.fftw 2>&1 && touch ../$(CLIMT_ARCH)/.configured_fftw

$(CLIMT_ARCH)/libfftw3_omp.a: $(CLIMT_ARCH)/.configured_fftw
	cd $(FFTW_DIR); make install > log.install.fftw 2>&1; cp ../lib/libfftw3_omp.a $(LIB_DIR)

$(CLIMT_ARCH)/libfftw3.a: $(CLIMT_ARCH)/.configured_fftw
	cp lib/libfftw3.a $(LIB_DIR)

fftw_lib: $(CLIMT_ARCH)/libfftw3_omp.a $(CLIMT_ARCH)/libfftw3.a


# SHTNS Configuration

$(CLIMT_ARCH)/.configured_sht:
	if [ ! -d $(SHT_DIR) ]; then tar -xvzf $(SHT_GZ) ; fi; cd $(SHT_DIR); ./configure CFLAGS=$(CLIMT_CFLAGS) LDFLAGS=$(LDFLAGS) --prefix=$(BASE_DIR) --disable-python --enable-openmp && touch ../$(CLIMT_ARCH)/.configured_sht

$(CLIMT_ARCH)/libshtns_omp.a: $(CLIMT_ARCH)/.configured_sht
	cd $(SHT_DIR); make install > log 2>&1; cp ../lib/libshtns_omp.a $(LIB_DIR); ls -ltrh $(LIB_DIR)

sht_lib: $(CLIMT_ARCH)/libshtns_omp.a

# BLAS Configuration

blas_lib: $(CLIMT_ARCH)/libopenblas.a

$(CLIMT_ARCH)/libopenblas.a:
	if [ ! -d $(BLAS_DIR) ]; then tar -xvzf $(BLAS_GZ) > log 2>&1; fi; cd $(BLAS_DIR); CFLAGS=$(CLIMT_CFLAGS) make NO_SHARED=1 NO_LAPACK=0 > log 2>&1 ; make PREFIX=$(BASE_DIR) NO_SHARED=1 NO_LAPACK=0 install > log.install.blas 2>&1 ; cp ../lib/libopenblas.a $(LIB_DIR)

# GFS Configuration

gfs_lib: $(CLIMT_ARCH)/libgfs_dycore.a

$(CLIMT_ARCH)/libgfs_dycore.a: $(GFS_DIR)/*.f90
	cd $(GFS_DIR); make python_lib; cp libgfs_dycore.a $(LIB_DIR); touch ../../_components/gfs/_gfs_dynamics.pyx

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
	cd $(GFS_DIR); make clean
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
