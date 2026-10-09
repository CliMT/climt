from netCDF4 import Dataset
import numpy as np
import time, uuid, argparse
import urllib.request, json
# ---------------------------------------------------------------------------------
# Copy a variable and all its attributes from one netCDF file to another
#
def copyVar(nc_in, nc_out, name, newname=None) :
	if newname is None :
		newname = name
	nc_out.createVariable(newname, nc_in.variables[name].dtype, nc_in.variables[name].dimensions)
	nc_out.variables[newname].setncatts(nc_in.variables[name].__dict__)
	nc_out.variables[newname][:] = nc_in.variables[name][:]
# ---------------------------------------------------------------------------------



#
# Check that source_id is valid 
#   Use source_id to obtain other text
#   https://raw.githubusercontent.com/PCMDI/cmip6-cmor-tables/master/Tables/CMIP6_CV.json
def create_attributes(args):
    with urllib.request.urlopen("file:CMIP6_CV.json") as url:
        cmip6 = json.loads(url.read().decode())


    if(args.source_id in cmip6['CV']['source_id'].keys()):
        source_id      = args.source_id
        source         = cmip6['CV']['source_id'][source_id]['source']
        institution_id = cmip6['CV']['source_id'][source_id]['institution_id'][0]
        institution    = cmip6['CV']['institution_id'][institution_id]
        physics_index = np.int32(args.physics_index)
        forcing_index = np.int32(args.forcing_index)
    else:
        print("source_id {} is not in CMIP6 Controlled Vocabulary".format(args.source_id))
        sys.exit(1)

    if (forcing_index < 1 ) or (forcing_index > 3):
        print('forcing_index must be 1, 2, or 3 (1 = all available greenhouse gases; 2 =  CO2, CH4, N2O, CFC12, CFC11eq;  3 = CO2, CH4, N2O, CFC12eq, HFC-134eq)')
        sys.exit(1)

    if (physics_index < 1 ):
        print('physics_index must be positive')
        sys.exit(1)

    #
    # Model/institution specific attributes
    #
    if physics_index == 1:
        variant_info_lw="Spectral file: sp_lw_ga7"
        variant_info_sw="Spectral file: sp_sw_ga7"
    elif physics_index == 2:
        variant_info_lw="Spectral file: sp_lw_300_jm2"
        variant_info_sw="Spectral file: sp_sw_260_jm2"
    elif physics_index == 3:
        variant_info_lw="Spectral file: sp_lw_ga9"
        variant_info_sw="Spectral file: sp_sw_ga9"
    elif physics_index == 4:
        variant_info_lw="Spectral file: sp_lw_300_jm3"
        variant_info_sw="Spectral file: sp_sw_260_jm3"
    
    variant_label  = "r1i1p{0}f{1}".format(physics_index,forcing_index)
    model_attrs = {
    "institution_id"  :institution_id,
    "institution"     :institution,
    "source_id"       :source_id,
    "source"          :source_id,
    "further_info_url":"https://furtherinfo.es-doc.org/CMIP6." + institution_id + "." + source_id + ".rad-irf.none." + variant_label,
    "forcing_index"   :np.int32(forcing_index),
    "license"         :"CMIP6 model data produced by " + institution_id + " is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License (https://creativecommons.org/licenses). " +
                     "Consult https://pcmdi.llnl.gov/CMIP6/TermsOfUse for terms of use governing CMIP6 output, including citation requirements and proper acknowledgment. " +
                     "Further information about this data, including some limitations, can be found via the further_info_url (recorded as a global attribute in this file) and at https://pcmdi.llnl.gov/." +
                     "The data producers and data providers make no warranty, either express or implied, including, but not limited to, warranties of merchantability and fitness for a particular purpose." +
                     "All liabilities arising from the supply of the information (including any liability arising in negligence) are excluded to the fullest extent permitted by law." }

    #
    # Required attributes, uniform across submissions
    #
    std_attrs = {
    "data_specs_version":"01.00.29",
    "physics_index":np.int32(physics_index)}

    # Submission attrs
    sub_attrs_lw = {
        'creation_date':time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()),
        'tracking_id'  : '/'.join(['hdl:21.14100',str(uuid.uuid4())]),
        "variant_label":variant_label,
            "variant_info":variant_info_lw}
    sub_attrs_sw = {
        'creation_date':time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()),
        'tracking_id'  : '/'.join(['hdl:21.14100',str(uuid.uuid4())]),
        "variant_label":variant_label,
            "variant_info":variant_info_sw}

    # Attributes are taken from https://docs.google.com/document/d/1h0r8RZr_f3-8egBMMh7aqLwy3snpD6_MrDz1q8n5XUk/edit
    # Data reference syntax attributes
    drs_attrs = {
    "activity_id"  :"RFMIP",   # (from CMIP6_activity_id.json)
    "product"      :"model-output",
    "experiment_id":"rad-irf", # (from CMIP6_experiment_id.json)
    "table_id"     :"Efx",     # (per http://clipc-services.ceda.ac.uk/dreq/u/efc0de22-5629-11e6-9079-ac72891c3257.html)
    "frequency"    :"fx",
    "sub_experiment_id":"none"}

    expt_attrs = {
    "Conventions"         :"CF-1.7 CMIP-6.2",
    "mip_era"             :"CMIP6",
    "experiment"          :"offline assessment of radiative transfer parmeterizations in clear skies",
    "sub_experiment"      :"none",
    "product"             :"model-output",
    "realization_index"   :np.int32(1),
    "initialization_index":np.int32(1),
    "source_type"         :"RAD",
    "nominal_resolution"  :"10 km",
    "realm"               :"atmos",
    "grid_label"          :"gn",
    "grid"                :"columns sampled from ERA-Interim, radiative fluxes computed independently"}

    return source_id, physics_index, variant_label, model_attrs, std_attrs, sub_attrs_lw, sub_attrs_sw, drs_attrs, expt_attrs

def write_files(atmos_file, source_id, physics_index, variant_label, drs_attrs, std_attrs, expt_attrs, model_attrs, sub_attrs_lw, sub_attrs_sw):
    nlevs = len(atmos_file.dimensions['level'])
    nlays = len(atmos_file.dimensions['layer'])
    nprof = len(atmos_file.dimensions['site'])
    nexpt = len(atmos_file.dimensions['expt'])
    print('Number of levels: ', nlevs)
    print('Number of layers: ', nlays)
    print('Number of profiles: ', nprof)
    print('Number of experiments: ', nexpt)

    expt_label = atmos_file.variables['expt_label'][:]
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
    # Read in fluxes
    lw_down = np.zeros((nexpt,nprof,nlevs))
    lw_up   = np.zeros((nexpt,nprof,nlevs))
    sw_down = np.zeros((nexpt,nprof,nlevs))
    sw_up   = np.zeros((nexpt,nprof,nlevs))

    if physics_index==1:
        for i in np.arange(nexpt):
            basename=str(label[i])
            lw_down[i,:,:] = Dataset(basename+'_lw_ga7.dflx').variables['dflx'][:,:,0].T
            lw_up[i,:,:]   = Dataset(basename+'_lw_ga7.uflx').variables['uflx'][:,:,0].T
            sw_down[i,:,:] = Dataset(basename+'_sw_ga7.vflx').variables['vflx'][:,:,0].T
            sw_up[i,:,:]   = Dataset(basename+'_sw_ga7.uflx').variables['uflx'][:,:,0].T
    elif physics_index==2:
        for i in np.arange(nexpt):
            basename=str(label[i])
            lw_down[i,:,:] = Dataset(basename+'_lw_300.dflx').variables['dflx'][:,:,0].T
            lw_up[i,:,:]   = Dataset(basename+'_lw_300.uflx').variables['uflx'][:,:,0].T
            sw_down[i,:,:] = Dataset(basename+'_sw_260.vflx').variables['vflx'][:,:,0].T
            sw_up[i,:,:]   = Dataset(basename+'_sw_260.uflx').variables['uflx'][:,:,0].T
    elif physics_index==3:
        for i in np.arange(nexpt):
            basename=str(label[i])
            lw_down[i,:,:] = Dataset(basename+'_lw_ga9.dflx').variables['dflx'][:,:,0].T
            lw_up[i,:,:]   = Dataset(basename+'_lw_ga9.uflx').variables['uflx'][:,:,0].T
            sw_down[i,:,:] = Dataset(basename+'_sw_ga9.vflx').variables['vflx'][:,:,0].T
            sw_up[i,:,:]   = Dataset(basename+'_sw_ga9.uflx').variables['uflx'][:,:,0].T
    elif physics_index==4:
        for i in np.arange(nexpt):
            basename=str(label[i])
            lw_down[i,:,:] = Dataset(basename+'_lw_ga9_300.dflx').variables['dflx'][:,:,0].T
            lw_up[i,:,:]   = Dataset(basename+'_lw_ga9_300.uflx').variables['uflx'][:,:,0].T
            sw_down[i,:,:] = Dataset(basename+'_sw_ga9_260.vflx').variables['vflx'][:,:,0].T
            sw_up[i,:,:]   = Dataset(basename+'_sw_ga9_260.uflx').variables['uflx'][:,:,0].T

    short_names = ['rlu','rsu', 'rld', 'rsd']
    stand_names = ['upwelling_longwave_flux_in_air','upwelling_shortwave_flux_in_air',
                'downwelling_longwave_flux_in_air','downwelling_shortwave_flux_in_air']
    variable_data = [lw_up,sw_up,lw_down,sw_down]
    
    for short, std, data in zip(short_names, stand_names,variable_data) :
        # File name is constructed per https://docs.google.com/document/d/1h0r8RZr_f3-8egBMMh7aqLwy3snpD6_MrDz1q8n5XUk/edit#
        # fixed strings are table_id, experiment_id, grid_label
        out_file_name = short + "_Efx_" + source_id + "_rad-irf_" + variant_label + "_gn" + ".nc"
        print('Creating ' + out_file_name)
        out_file = Dataset(out_file_name, mode='w', FORMAT='NETCDF4_CLASSIC')
        out_file.setncatts(drs_attrs)
        out_file.setncatts(std_attrs)
        out_file.setncatts(expt_attrs)
        out_file.setncatts(model_attrs)
        if "rl" in short:
            out_file.setncatts(sub_attrs_lw)
        elif "rs" in short:
            out_file.setncatts(sub_attrs_sw)
        out_file.setncatts({'variable_id'  :short})
        d = out_file.createDimension('expt',  atmos_file.dimensions['expt'].size)
        d = out_file.createDimension('site',  atmos_file.dimensions['site'].size)
        d = out_file.createDimension('level', atmos_file.dimensions['level'].size)
        copyVar(atmos_file, out_file, 'lat')
        copyVar(atmos_file, out_file, 'lon')
        copyVar(atmos_file, out_file, 'time')
        copyVar(atmos_file, out_file, 'pres_level', 'plev')
        v = out_file.createVariable(short,   'f4', ('expt', 'site', 'level'))
        v[:] = data[:]
        v.setncatts({'variable_id'  :short,
                    'standard_name':std,
                    'units'        :'W m-2',
                    #'_FillValue'   :np.float32(-1.e+03),
                    'missing_value':np.float32(-1.e+03),
                    "cell_methods" :"area: point",
                    'coordinates'  :'lon lat time'})
        copyVar(atmos_file, out_file, 'profile_weight')
        out_file.close()

def main():
    atmos_file = Dataset('multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc', mode='r')
    # Available from https://www.earthsystemcog.org/projects/rfmip/resources/
    # or from https://esgf-node.llnl.gov/search/input4mips/ ; search for "RFMIP"

    parser = argparse.ArgumentParser(description='Create CMIP6/ESGF-compliant output files for RFMIP-IRF.')
    parser.add_argument('--source_id', type=str, \
                        default = "HadGEM3-GC31-LL",
                        help='Source ID, must match CMIP Controlled Vocabulary at https://github.com/WCRP-CMIP/CMIP6_CVs/blob/master/CMIP6_source_id.json')
    parser.add_argument('--forcing_index', type=int, \
                        default = 3,      \
                        help='Forcing index (1 = all available greenhouse gases; 2 =  CO2, CH4, N2O, CFC11eq;  3 = CO2, CH4, N2O, CFC12eq, HFC-134eq)')
    parser.add_argument('--physics_index', type=int, \
                    default = 1,      \
                    help='Physics index: 1 = ga7; 2 = ga7_ref; 3 = ga9; 4 = ga9_ref')
    args = parser.parse_args()

    source_id, physics_index, variant_label, model_attrs, std_attrs, sub_attrs_lw, sub_attrs_sw, drs_attrs, expt_attrs = create_attributes(args)
    write_files(atmos_file, source_id, physics_index, variant_label, drs_attrs, std_attrs, expt_attrs, model_attrs, sub_attrs_lw, sub_attrs_sw)
    
if __name__ == "__main__":
    main()
