# Build RRTMG LW library

.SUFFIXES: .f90 .o

OBJS =     parkind.o\
    parrrtm.o\
    rrlw_cld.o\
    rrlw_con.o\
    rrlw_kg01.o\
    rrlw_kg02.o\
    rrlw_kg03.o\
    rrlw_kg04.o\
    rrlw_kg05.o\
    rrlw_kg06.o\
    rrlw_kg07.o\
    rrlw_kg08.o\
    rrlw_kg09.o\
    rrlw_kg10.o\
    rrlw_kg11.o\
    rrlw_kg12.o\
    rrlw_kg13.o\
    rrlw_kg14.o\
    rrlw_kg15.o\
    rrlw_kg16.o\
    rrlw_ncpar.o\
    rrlw_ref.o\
    rrlw_tbl.o\
    rrlw_vsn.o\
    rrlw_wvn.o\
    rrtmg_lw_cldprop.o\
    rrtmg_lw_cldprmc.o\
    rrtmg_lw_rtrn.o\
    rrtmg_lw_rtrnmr.o\
    rrtmg_lw_rtrnmc.o\
    rrtmg_lw_setcoef.o\
    rrtmg_lw_taumol.o\
    rrtmg_lw_rad.nomcica.o\
    mcica_random_numbers.o\
    rrtmg_lw_init.o\
    mcica_subcol_gen_lw.o\
    rrtmg_lw_rad.o\
    rrtmg_lw_c_binder.o

OBJ_K_TABLES = rrtmg_lw_k_g.o

all: lib_rrtmglw

lib_rrtmglw: $(OBJS) $(OBJ_K_TABLES)
	ar rcs librrtmg_lw.a  $(OBJS) $(OBJ_K_TABLES)

rrtmg_lw_k_g.o: rrtmg_lw_k_g.f90
	$(FC) $(FFLAGS) -O0 $< -c -o $@

%.o: %.mod

.f90.o:
	$(FC) $(FFLAGS) -O3 $< -c -o $@

clean:
	rm -f $(OBJS) $(OBJ_K_TABLES) *.mod *genmod* librrtmg_lw.a
