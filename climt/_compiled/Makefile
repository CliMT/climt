FFTW_DIR = fftw-3.3.6-pl2
SHT_DIR = shtns-2.7-r557
BLAS_DIR = OpenBLAS-0.2.19
GFS_DIR = GFS

FFTW_GZ = fftw-3.3.6-pl2.tar.gz
SHT_GZ = shtns-2.7-r557.tar.gz
BLAS_GZ = OpenBLAS-0.2.19.tar.gz

BASE_DIR = $(PWD)

CFLAGS = "-I$(BASE_DIR) -L$(BASE_DIR) -fPIC"
LDFLAGS = "-L$(BASE_DIR)"

all: fftw_lib sht_lib blas_lib gfs_lib

# FFTW3 Configuration

.configured_fftw:
	if [ ! -d $(FFTW_DIR) ]; then tar -xvzf $(FFTW_GZ); fi; cd $(FFTW_DIR); ./configure CFLAGS=$(CFLAGS) LDFLAGS=$(LDFLAGS) --prefix=$(BASE_DIR) --enable-openmp --enable-avx --enable-avx2 && cd -; touch .configured_fftw

lib/libfftw3_omp.a: .configured_fftw
	cd $(FFTW_DIR); make install

fftw_lib: lib/libfftw3_omp.a

# SHTNS Configuration

.configured_sht:
	if [ ! -d $(SHT_DIR) ]; then tar -xvzf $(SHT_GZ); fi; cd $(SHT_DIR); ./configure CFLAGS=$(CFLAGS) LDFLAGS=$(LDFLAGS) --prefix=$(BASE_DIR) --disable-python --enable-openmp && cd -; touch .configured_sht

lib/libshtns_omp.a: .configured_sht
	cd $(SHT_DIR); make install

sht_lib: lib/libshtns_omp.a

# BLAS Configuration

blas_lib: lib/libopenblas.a

lib/libopenblas.a:
	if [ ! -d $(BLAS_DIR) ]; then tar -xvzf $(BLAS_GZ); fi; cd $(BLAS_DIR); make NO_SHARED=1 NO_LAPACK=0; make PREFIX=$(BASE_DIR) NO_SHARED=1 NO_LAPACK=0 install

# GFS Configuration

gfs_lib: lib/libgfsDycore.a

lib/libgfsDycore.a:
	cd $(GFS_DIR); make python_lib; cp libgfsDycore.a ../lib

clean:
	if [ -d $(FFTW_DIR) ]; then rm -Rf $(FFTW_DIR); fi
	if [ -d $(SHT_DIR) ]; then rm -Rf $(SHT_DIR); fi
	if [ -d $(BLAS_DIR) ]; then rm -Rf $(BLAS_DIR); fi
	if [ -d  include ]; then rm -Rf include; fi
	if [ -d  share ]; then rm -Rf share; fi
	if [ -d  bin ]; then rm -Rf bin; fi
	if [ -d  lib ]; then rm -Rf lib; fi
	if [ -e  .configured_fftw ]; then rm .configured_fftw; fi
	if [ -e  .configured_sht ]; then rm .configured_sht; fi
	cd $(GFS_DIR); make clean
