# Build emanuel convection library


all: lib_emanuel

lib_emanuel: emanuel.o
	ar rcs libemanuel.a emanuel.o

emanuel.o: convect43c.f90
	$(FC) $(FFLAGS) -O3 convect43c.f90 -c -o emanuel.o

clean:
	rm -f emanuel.o *.mod *genmod* libemanuel.a	
