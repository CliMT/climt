OMPFLAG		 = -fopenmp
DEBUGFLAGS   =-fimplicit-none
GFS_FFLAGS	 = $(FFLAGS) -march=native -ffast-math -pipe -fomit-frame-pointer -fbacktrace
SHAREFLAGS 	 = #-static

.SUFFIXES: .o .f90 .F90

OBJS=	 params.o \
	 kinds.o\
	 dyn_init.o\
	 dyn_run.o\
	 physcons.o\
	 phy_data.o\
	 spectral_data.o\
	 grid_data.o\
	 pressure_data.o\
	 finalize.o\
	 dyn_finalize.o\
	 semimp_data.o\
	 run.o\
         phy_run.o\
         phy_init.o\
	 phy_finalize.o\
	 shtns.o
	 #init.o\

all: $(python_lib)

python_lib : $(OBJS)
	ar rcs libgfs_dycore.a $(OBJS) 

force_look:
	true

clean: force_look
	rm -f $(OBJS) *.mod *genmod* libgfs_dycore.a


.f90.o:
	$(FC) -Wall $(GFS_FFLAGS) $(DEBUGFLAGS) $(SHAREFLAGS) $(OMPFLAG) -c $<

.F90.o:
	$(FC) -Wall $(GFS_FFLAGS) $(DEBUGFLAGS)  $(SHAREFLAGS) $(OMPFLAG) -c $<

.f.o:
	$(FC) -Wall $(GFS_FFLAGS) $(DEBUGFLAGS) $(SHAREFLAGS)  $(OMPFLAG) -free -c $<

%.o : %.mod

kinds.o:	kinds.f90

physcons.o:	physcons.f90 kinds.o

constants.o:	constants.f90 kinds.o

shtns.o:	shtns.f90 kinds.o

params.o:	params.f90 kinds.o 

init.o:  	init.f90 params.o dyn_init.o phy_init.o phy_data.o physcons.o kinds.o 

finalize.o:  	finalize.f90 params.o dyn_finalize.o phy_finalize.o 

spectral_data.o:  	spectral_data.f90 kinds.o

phy_data.o:  	phy_data.f90 kinds.o params.o physcons.o

semimp_data.o: semimp_data.f90 kinds.o shtns.o pressure_data.o physcons.o

grid_data.o:  	grid_data.f90 kinds.o params.o

pressure_data.o:  	pressure_data.f90 kinds.o params.o physcons.o grid_data.o

dyn_init.o:  	dyn_init.f90 params.o shtns.o spectral_data.o pressure_data.o physcons.o grid_data.o semimp_data.o 
		$(FC) -Wall $(GFS_FFLAGS) $(OMPFLAG) -c dyn_init.f90

phy_init.o:     phy_init.f90 phy_data.o

phy_finalize.o:     phy_finalize.f90 

phy_run.o:      phy_run.f90 pressure_data.o shtns.o spectral_data.o physcons.o kinds.o params.o grid_data.o 
		$(FC) -Wall $(GFS_FFLAGS) $(OMPFLAG) -c phy_run.f90
	        
dyn_finalize.o:  dyn_finalize.f90 shtns.o spectral_data.o pressure_data.o grid_data.o 

dyn_run.o:	dyn_run.f90 pressure_data.o shtns.o spectral_data.o physcons.o kinds.o params.o grid_data.o 

run.o:   	run.f90 dyn_run.o phy_run.o kinds.o params.o  

