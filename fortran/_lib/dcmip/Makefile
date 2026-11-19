# Build simple physics library


all: libdcmip

libdcmip: test4.o test5.o
	ar rcs libdcmip.a test4.o test5.o

test4.o: dcmip_initial_conditions_test_4_v3.f90 
	$(FC) $(FFLAGS) -O3 dcmip_initial_conditions_test_4_v3.f90  -c -o test4.o

test5.o: dcmip_initial_conditions_test_5_v1.f90 
	$(FC) $(FFLAGS) -O3 dcmip_initial_conditions_test_5_v1.f90  -c -o test5.o
clean:
	rm -f test4.o test5.o *.mod *genmod* libdcmip.a
