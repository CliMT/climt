# Build simple physics library


all: lib_simple_physics

lib_simple_physics: simple_physics_custom.o
	ar rcs libsimple_physics.a simple_physics_custom.o

simple_physics_custom.o: simple_physics_custom.f90
	$(FC) $(FFLAGS) -O3 simple_physics_custom.f90 -c -o simple_physics_custom.o

clean:
	rm -f simple_physics_custom.o *.mod *genmod* libsimple_physics.a	
