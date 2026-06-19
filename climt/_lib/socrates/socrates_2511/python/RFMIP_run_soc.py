# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
'''
Converts the netCDF files from RFMIP into netCDF input files for the SOCRATES code.
Calculates radiative fluxes using the GA7 or GA9 and reference spectral files.
Creates netCDF files for the output fluxes to submit to RFMIP.
'''

from netCDF4 import Dataset
import nctools as nc
import numpy as np
import sys, math, subprocess, time

def copyVar(nc_in, nc_out, name, newname=None) :   
	if newname is None : 
		newname = name
	nc_out.createVariable(newname, nc_in.variables[name].dtype, nc_in.variables[name].dimensions) 
	nc_out.variables[newname].setncatts(nc_in.variables[name].__dict__)
	nc_out.variables[newname][:] = nc_in.variables[name][:] 

def write_dim(ncdf_file, n_val, vals, name, vtype):
    dimension = ncdf_file.createDimension(name, n_val)
    variable = ncdf_file.createVariable(name, vtype, name)
    variable[:] = vals

def write_var(ncdf_file, vals, name, vtype, dims, units, standard_name, longname, coords):
    variable = ncdf_file.createVariable(name, vtype, dims)
    if (type(units) is str):
        variable.units = units
    if (type(standard_name) is str):
        variable.standard_name = standard_name
    if (type(longname) is str):
        variable.long_name = longname
    if (type(coords) is str):
        variable.coordinates = coords
    variable[:] = vals
    
def main_socrates_run(infile, configuration, equivalence):
    
    rfdat = Dataset(infile)
    basename = 'rfmip'
    nlevs = len(rfdat.dimensions['level'])
    nlays = len(rfdat.dimensions['layer'])
    nprof = len(rfdat.dimensions['site'])
    nexpt = len(rfdat.dimensions['expt'])
    print('Number of levels: ', nlevs)
    print('Number of layers: ', nlays)
    print('Number of profiles: ', nprof)
    print('Number of experiments: ', nexpt)

    expt_label    = rfdat.variables['expt_label'][:]
    p             = rfdat.variables['pres_layer'][:,:].T # layer, profile
    pl            = rfdat.variables['pres_level'][:,:].T # levels, profile
    lwalb         = 1.0 - rfdat.variables['surface_emissivity'][:] # profile
    swalb         = rfdat.variables['surface_albedo'][:] # profile
    szen          = rfdat.variables['solar_zenith_angle'][:] # profile
    stoa          = rfdat.variables['total_solar_irradiance'][:] # profile
    t             = rfdat.variables['temp_layer'][:,:,:].T # layer, profile, expt
    tl            = rfdat.variables['temp_level'][:,:,:].T # levels, profile, expt
    tstar         = rfdat.variables['surface_temperature'][:,:].T # profile, expt
    h2o_mol       = (rfdat.variables['water_vapor'][:,:,:].T # layer, profile, expt
             * float(rfdat.variables['water_vapor'].units))
    co2_mol       = (rfdat.variables['carbon_dioxide_GM'][:] # expt
             * float(rfdat.variables['carbon_dioxide_GM'].units))
    o3_mol        = (rfdat.variables['ozone'][:,:,:].T # layer, profile, expt
             * float(rfdat.variables['ozone'].units))
    n2o_mol       = (rfdat.variables['nitrous_oxide_GM'][:] # expt
             * float(rfdat.variables['nitrous_oxide_GM'].units))
    ch4_mol       = (rfdat.variables['methane_GM'][:] # expt
             * float(rfdat.variables['methane_GM'].units))
    o2_mol        = (rfdat.variables['oxygen_GM'][:] # expt
             * float(rfdat.variables['oxygen_GM'].units))
    cfc11_mol     = (rfdat.variables['cfc11_GM'][:] # expt
             * float(rfdat.variables['cfc11_GM'].units))
    cfc11eq_mol   = (rfdat.variables['cfc11eq_GM'][:] # expt
             * float(rfdat.variables['cfc11eq_GM'].units))
    cfc12_mol     = (rfdat.variables['cfc12_GM'][:] # expt
             * float(rfdat.variables['cfc12_GM'].units))
    cfc12eq_mol   = (rfdat.variables['cfc12eq_GM'][:] # expt
             * float(rfdat.variables['cfc12eq_GM'].units))
    cfc113_mol    = (rfdat.variables['cfc113_GM'][:] # expt
             * float(rfdat.variables['cfc113_GM'].units))
    hcfc22_mol    = (rfdat.variables['hcfc22_GM'][:] # expt
             * float(rfdat.variables['hcfc22_GM'].units))
    hfc134a_mol   = (rfdat.variables['hfc134a_GM'][:] # expt
             * float(rfdat.variables['hfc134a_GM'].units))
    hfc134aeq_mol = (rfdat.variables['hfc134aeq_GM'][:] # expt
             * float(rfdat.variables['hfc134aeq_GM'].units))

    # Convert molar ratios to mass mixing ratios
    mol_weight_air = 28.966
    h2o       = h2o_mol       * 18.0153   / mol_weight_air
    co2       = co2_mol       * 44.0100   / mol_weight_air
    o3        = o3_mol        * 47.9982   / mol_weight_air
    n2o       = n2o_mol       * 44.0128   / mol_weight_air
    ch4       = ch4_mol       * 16.0430   / mol_weight_air
    o2        = o2_mol        * 31.9988   / mol_weight_air
    cfc11     = cfc11_mol     * 137.3686  / mol_weight_air
    cfc11eq   = cfc11eq_mol   * 137.3686  / mol_weight_air
    cfc12     = cfc12_mol     * 120.9140  / mol_weight_air
    cfc12eq   = cfc12eq_mol   * 120.9140  / mol_weight_air
    cfc113    = cfc113_mol    * 187.3765  / mol_weight_air
    hcfc22    = hcfc22_mol    * 86.46892  / mol_weight_air
    hfc134a   = hfc134a_mol   * 102.03184 / mol_weight_air
    hfc134aeq = hfc134aeq_mol * 102.03184 / mol_weight_air

    # Set unlit points to zero irradiance
    unlit = szen > 90.0
    szen[unlit] = 90.0 
    stoa[unlit] = 0.0

    # Output 3D fields using dummy dimensions:
    pstar = np.array(pl[nlevs-1,2])
    pdim = p[:,2]
    pldim = pl[:,2]
    expt = np.arange(nexpt)
    lon = np.arange(nexpt) + 1
    lat = np.arange(nprof) + 1
    basis = 1

    # Create separate files for each experiment:
    label=expt_label.copy()
    label[0]='PD'
    label[1]='PI'
    label[3]='future'
    label[8]='PI_CO2'
    label[9]='PI_CH4'
    label[10]='PI_N2O'
    label[11]='PI_O3'
    label[12]='PI_HCs'
    label[13]='p4K'
    label[14]='p4KconRH'
    label[15]='PI_all'
    label[16]='future_all'
    label[17]='LGM'

    for i in np.arange(nexpt):
        lon = i+1
        basename=label[i]
        subprocess.call('rm '+basename+'.*', shell=True)

        # Pressure:
        nc.ncout3d(basename + '.p', lon, lat, pdim, p, longname= 'Pressure', units= 'Pa')
        nc.ncout3d(basename + '.pl', lon, lat, pldim, pl, longname= 'Pressure', units= 'Pa')

        # Temperature:
        nc.ncout3d(basename + '.tstar', lon, lat, pstar, tstar[:,i], longname= 'Surface temperature', units= 'K')
        nc.ncout3d(basename + '.t', lon, lat, pdim, t[:,:,i], longname= 'Temperature', units= 'K')
        nc.ncout3d(basename + '.tl', lon, lat, pldim, tl[:,:,i], longname= 'Temperature', units= 'K')

        # Solar beam:
        nc.ncout2d(basename + '.stoa', lon, lat, stoa, longname= 'Solar Irradiance', units= 'WM-2')
        nc.ncout2d(basename + '.szen', lon, lat, szen, longname= 'Solar Zenith Angle', units= 'degrees')

        # Surface albedos:
        nc.ncout_surf(basename + '.surflw', lon, lat, basis, lwalb)
        nc.ncout_surf(basename + '.surfsw', lon, lat, basis, swalb)

        # Gases:
        nc.ncout3d(basename + '.q',       lon, lat, pdim, h2o[:,:,i],   longname = 'Specific Humidity')
        nc.ncout3d(basename + '.co2',     lon, lat, pdim, co2[i],       longname = 'Carbon Dioxide MMR')
        nc.ncout3d(basename + '.o3',      lon, lat, pdim, o3[:,:,i],    longname = 'Ozone MMR')
        nc.ncout3d(basename + '.n2o',     lon, lat, pdim, n2o[i],       longname = 'Dinitrogen oxide MMR')
        nc.ncout3d(basename + '.ch4',     lon, lat, pdim, ch4[i],       longname = 'Methane MMR')
        nc.ncout3d(basename + '.o2',      lon, lat, pdim, o2[i],        longname = 'Oxygen MMR')
        if equivalence == 2:
            print("Using Equivalence 2: CFC12 + CFC11eq")
            nc.ncout3d(basename + '.cfc12',   lon, lat, pdim, cfc12[i],     longname = 'CFC-12 MMR')
            nc.ncout3d(basename + '.cfc11', lon, lat, pdim, cfc11eq[i], longname = 'CFC-11 MMR')
        elif equivalence == 3:
            print("Using Equivalence 3: CFC12eq + HFC134aeq")
            nc.ncout3d(basename + '.cfc12',   lon, lat, pdim, cfc12eq[i],     longname = 'CFC-12 MMR')
            nc.ncout3d(basename + '.hfc134a', lon, lat, pdim, hfc134aeq[i], longname = 'HFC-134a MMR')
        else:
            print("Set Equivalence 2 or 3. Equivalence 1 currently not supported by SOCRATES")
            sys.exit(1)
        #    nc.ncout3d(basename + '.cfc11',   lon, lat, pdim, cfc11[i],     longname = 'CFC-11 MMR')
        
        #    nc.ncout3d(basename + '.cfc113',  lon, lat, pdim, cfc113[i],    longname = 'CFC-113 MMR')
        #    nc.ncout3d(basename + '.hcfc22',  lon, lat, pdim, hcfc22[i],    longname = 'HCFC-22 MMR')
        #    nc.ncout3d(basename + '.hfc134a', lon, lat, pdim, hfc134a[i],   longname = 'HFC-134a MMR')
        #    nc.ncout3d(basename + '.cfc12',   lon, lat, pdim, cfc12eq[i],   longname = 'CFC-12 MMR')
     
    
    if configuration == 'ga7':
        print('Running GA7 configuration: LW and SW')
    elif configuration == 'ga7_ref':
        print('Running GA7 reference configuration: LW (300 band) and SW (260 band)')
    elif configuration == 'ga9':
        print('Running GA9 configuration: LW and SW')
    elif configuration == 'ga9_ref':
        print('Running GA9 reference configuration: LW (300 band) and SW (260 band)')
    else:
        print('Configurations currently supported: ga7, ga7_ref, ga9, ga9_ref')
        sys.exit(1)
    subprocess.call('date')
    procs=[]
    for i in np.arange(nexpt):
        basename=str(label[i])
        if configuration == 'ga7':
            procs.append( subprocess.Popen(
                'resrm '+basename+' && '+
                'cp '+basename+'.surflw '+basename+'.surf && '+
                'Cl_run_cdf -B '+basename+' -s $RAD_DATA/spectra/ga7/sp_lw_ga7 -q -R 1 9 -I -g 4 -c -C 5 -z 2 && '+
                'fmove '+basename+' '+basename+'_lw_ga7 && '+
                'cp '+basename+'.surfsw '+basename+'.surf && '+
                'Cl_run_cdf -B '+basename+' -s $RAD_DATA/spectra/ga7/sp_sw_ga7 -R 1 6 -S -r -g 4 -c -C 5 && '+
                'fmove '+basename+' '+basename+'_sw_ga7 && '+
                'echo "'+basename+': '+expt_label[i]+'"', shell=True) )
        elif configuration == 'ga7_ref':
            procs.append( subprocess.Popen(
                'resrm '+basename+' && '+
                'cp '+basename+'.surflw '+basename+'.surf && '+
                'Cl_run_cdf -B '+basename+' -s $RAD_DATA/spectra/ga7_ref/sp_lw_300_jm2 -q -R 1 300 -I -g 2 -c -C 5 -z 2 && ' +
                'fmove '+basename+' '+basename+'_lw_300 && '+
                'cp '+basename+'.surfsw '+basename+'.surf && '+
                'Cl_run_cdf -B '+basename+' -s $RAD_DATA/spectra/ga7_ref/sp_sw_260_jm2 -R 1 260 -S -r -g 2 -c -C 5 && ' +
                'fmove '+basename+' '+basename+'_sw_260 && '+
                'echo "'+basename+': '+expt_label[i]+'"', shell=True) )
        elif configuration == 'ga9':
            procs.append( subprocess.Popen(
                'resrm '+basename+' && '+
                'cp '+basename+'.surflw '+basename+'.surf && '+
                'Cl_run_cdf -B '+basename+' -s $RAD_DATA/spectra/ga9/sp_lw_ga9 -q -R 1 9 -I -g 4 -c -C 5 -z 2 && '+
                'fmove '+basename+' '+basename+'_lw_ga9 && '+
                'cp '+basename+'.surfsw '+basename+'.surf && '+
                'Cl_run_cdf -B '+basename+' -s $RAD_DATA/spectra/ga9/sp_sw_ga9 -R 1 6 -S -r -g 4 -c -C 5 && '+
                'fmove '+basename+' '+basename+'_sw_ga9 && '+
                'echo "'+basename+': '+expt_label[i]+'"', shell=True) )
        elif configuration == 'ga9_ref':
            procs.append( subprocess.Popen(
                'resrm '+basename+' && '+
                'cp '+basename+'.surflw '+basename+'.surf && '+
                'Cl_run_cdf -B '+basename+' -s $RAD_DATA/spectra/ga9_ref/sp_lw_300_jm3 -q -R 1 300 -I -g 2 -c -C 5 -z 2 && ' +
                'fmove '+basename+' '+basename+'_lw_ga9_300 && '+
                'cp '+basename+'.surfsw '+basename+'.surf && '+
                'Cl_run_cdf -B '+basename+' -s $RAD_DATA/spectra/ga9_ref/sp_sw_260_jm3 -R 1 260 -S -r -g 2 -c -C 5 && ' +
                'fmove '+basename+' '+basename+'_sw_ga9_260 && '+
                'echo "'+basename+': '+expt_label[i]+'"', shell=True) )
        else:
            print('Configurations currently supported: ga7, ga7_ref, ga9, ga9_ref')
            sys.exit(1)

        # Simple parallel processing (run 6 experiments at once)
        running = 6
        while running > 5:
            running = 0
            for proc in procs:
                if proc.poll() is None:
                    running +=1
            time.sleep(0.1)

    # Ensure all processes have finished
    for proc in procs:
        proc.wait()
    subprocess.call('date')
    
if __name__ == '__main__':
    if (len(sys.argv) > 2):
        infile = sys.argv[1]
        configuration = sys.argv[2]
        equivalence = 3 # Forcing eqivalence for CFCs; SOCRATES default is 3.
        main_socrates_run(infile, configuration, equivalence)
    elif (len(sys.argv) > 1):
        infile = sys.argv[1]
        configuration = 'ga7'
        equivalence = 3 # Forcing eqivalence for CFCs; SOCRATES default is 3.
        main_socrates_run(infile, configuration, equivalence)
    else:
        raise RuntimeError('\nplease enter an RFMIP netCDF filename \n'
                           'e.g. python RFMIP_run_soc.py multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc\n'
                           'wget http://aims3.llnl.gov/thredds/fileServer/user_pub_work/input4MIPs/CMIP6/RFMIP/UColorado/UColorado-RFMIP-1-2/atmos/fx/multiple/none/v20190401/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc')





