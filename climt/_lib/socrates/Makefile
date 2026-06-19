# Makefile to compile socrates and runes

ifeq ($(origin FC), default)
    FC = gfortran
endif

ifeq ($(origin FFLAGS), default)
    FFLAGS = -O3 -fPIC
endif

BASE_DIR = $(PWD)
FILES_DIR = $(BASE_DIR)/socrates_2511/bin
PLATFORM_LIB_DIR = $(BASE_DIR)/../$(CLIMT_ARCH)

# NetCDF detection
NF_CONFIG := $(shell which nf-config 2>/dev/null)
NC_CONFIG := $(shell which nc-config 2>/dev/null)

ifneq ($(NF_CONFIG),)
    INCCDF_PATH := $(shell nf-config --includedir)
    LIBCDF_PATH := $(shell nf-config --flibs | grep -o '\-L[^ ]*' | head -n1 | sed 's/-L//')
    ifeq ($(LIBCDF_PATH),)
        LIBCDF_PATH := $(shell nf-config --prefix)/lib
    endif
    LIBCDF_NAME := netcdff
    BUILD_TARGET := all
else ifneq ($(NC_CONFIG),)
    INCCDF_PATH := $(shell nc-config --includedir)
    LIBCDF_PATH := $(shell nc-config --libs | grep -o '\-L[^ ]*' | head -n1 | sed 's/-L//')
    ifeq ($(LIBCDF_PATH),)
        LIBCDF_PATH := $(shell nc-config --prefix)/lib
    endif
    LIBCDF_NAME := netcdf
    BUILD_TARGET := all
else
    INCCDF_PATH := /usr/include
    LIBCDF_PATH := /usr/lib
    LIBCDF_NAME := netcdff
    BUILD_TARGET := radlib.a cklib.a scatlib.a cosplib.a flexchemlib.a
endif

.PHONY: all clean compile_librunes

all: compile_librunes

compile_librunes:
	# Clean any existing socrates build directory
	rm -rf socrates_2511
	# Extract the tarball fresh
	tar -xf socrates_2511.tar.xz
	# In-place modify make/Mk_cmd in the source tree to configure the compiler, flags, and NetCDF paths
	sed -i 's|^FORTCOMP\s*=.*|FORTCOMP = $(FC) -fPIC -O -c|' socrates_2511/make/Mk_cmd
	sed -i 's|^LINK\s*=.*|LINK = $(FC) -fPIC -O|' socrates_2511/make/Mk_cmd
	sed -i 's|^INCCDF_PATH\s*=.*|INCCDF_PATH = $(INCCDF_PATH)|' socrates_2511/make/Mk_cmd
	sed -i 's|^LIBCDF_PATH\s*=.*|LIBCDF_PATH = $(LIBCDF_PATH)|' socrates_2511/make/Mk_cmd
	sed -i 's|^LIBCDF_NAME\s*=.*|LIBCDF_NAME = $(LIBCDF_NAME)|' socrates_2511/make/Mk_cmd

	# Run build_code to compile socrates
	cd socrates_2511 && ./build_code $(BUILD_TARGET)
	# Copy our wrapper files into bin/ (where all compile occurs) to be absolutely sure
	cp patchfiles/runes_driver.f90 $(FILES_DIR)/runes_driver_copy.f90
	cp patchfiles/runeswrap.f90 $(FILES_DIR)/runeswrap.f90
	# Compile wrapper files and pack everything into librunes.a
	cd $(FILES_DIR) && \
		$(FC) $(FFLAGS) -I. -c runes_driver_copy.f90 && \
		$(FC) $(FFLAGS) -I. -c runeswrap.f90 && \
		ar rcs librunes.a *.o
	# Copy librunes.a to the platform directory
	mkdir -p $(PLATFORM_LIB_DIR)
	cp $(FILES_DIR)/librunes.a $(PLATFORM_LIB_DIR)/

clean:
	rm -rf socrates_2511
	rm -f $(PLATFORM_LIB_DIR)/librunes.a
