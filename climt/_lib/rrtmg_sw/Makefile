# Build RRTMG LW library

.SUFFIXES: .f90 .o

OBJS =     parkind.o\
    parrrsw.o\
    rrsw_cld.o\
    rrsw_con.o\
    rrsw_kg16.o\
    rrsw_kg17.o\
    rrsw_kg18.o\
    rrsw_kg19.o\
    rrsw_kg20.o\
    rrsw_kg21.o\
    rrsw_kg22.o\
    rrsw_kg23.o\
    rrsw_kg24.o\
    rrsw_kg25.o\
    rrsw_kg26.o\
    rrsw_kg27.o\
    rrsw_kg28.o\
    rrsw_kg29.o\
    rrsw_ncpar.o\
    rrsw_ref.o\
    rrsw_tbl.o\
    rrsw_vsn.o\
    rrsw_aer.o\
    rrsw_wvn.o\
    rrtmg_sw_cldprop.o\
    rrtmg_sw_cldprmc.o\
    rrtmg_sw_taumol.o\
    rrtmg_sw_vrtqdr.o\
    rrtmg_sw_reftra.o\
    rrtmg_sw_spcvmc.o\
    rrtmg_sw_setcoef.o\
    rrtmg_sw_spcvrt.o\
    rrtmg_sw_rad.nomcica.o\
    mcica_random_numbers.o\
    rrtmg_sw_init.o\
    mcica_subcol_gen_sw.o\
    rrtmg_sw_rad.o\
    rrtmg_sw_c_binder.o

OBJ_K_TABLES = rrtmg_sw_k_g.o
 
all: lib_rrtmgsw

lib_rrtmgsw: $(OBJS) $(OBJ_K_TABLES)
	ar rcs librrtmg_sw.a  $(OBJS) $(OBJ_K_TABLES)

rrtmg_sw_k_g.o: rrtmg_sw_k_g.f90
	$(FC) $(FFLAGS) -O0 $< -c -o $@

.f90.o:
	$(FC) $(FFLAGS) -O3 $< -c -o $@

clean:
	rm -f $(OBJS) $(OBJ_K_TABLES) *.mod *genmod* librrtmg_sw.a
