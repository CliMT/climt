# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
'''
This file contains python functions that can be used to produce
multicolumn netCDF data compatible with l_run_cdf.
'''

from netCDF4 import Dataset
import numpy as np
import re

def write_dim(ncdf_file, n_val, vals, name, vtype, dims, units, title):
    dimension = ncdf_file.createDimension(name, n_val)
    variable = ncdf_file.createVariable(name, vtype, dims)
    variable.units = units
    variable.title = title
    variable[:] = vals

def write_var(ncdf_file, vals, name, vtype, dims, units, title):
    variable = ncdf_file.createVariable(name, vtype, dims)
    if (type(units) is str):
        variable.units = units
    if (type(title) is str):
        variable.title = title
    variable[:] = vals

def ncout_surf(file, lon, lat, basis, alb):

# ; Program to create netCDF files of surface albedo weights.
# ; Normally 'file' will have the extension .surf, and 'basis' = 1.
# ; 'alb' should then be an array of surface albedo values or a single value.

    levels = np.size(basis)
    n_lon = np.size(lon)
    n_lat = np.size(lat)
    nvals = np.size(alb)

    if (nvals == 1):
        albs = np.zeros(n_lon*n_lat*levels)
        albs += np.sum(alb)
    elif (nvals == n_lon*n_lat*levels):
        albs = alb
    else:
        raise RuntimeError(' Error in ncout_surf: arrays dont match', nvals
                            , n_lon * n_lat * levels)

    print('ncout_surf - file: ',file)

    ncdf_file = create_cdf(file)

    write_dim(ncdf_file, n_lon, lon, 'lon', 'f4', 'lon', 'degree', 'LONGITUDE')
    write_dim(ncdf_file, n_lat, lat, 'lat', 'f4', 'lat', 'degree', 'LATITUDE')
    write_dim(ncdf_file, levels, basis, 'basis', 'i2', 'basis', 'None',
              'BASIS FUNCTION')

    write_var(ncdf_file, albs, 'alb', 'f4', ('basis', 'lat', 'lon'),'None', 
              'ALBEDO WEIGHTS')

    ncdf_file.close()

#-------------------------------------------------------------------------------

def ncout_spectral_surf(file, lon, lat, bands, alb):

# Program to create netCDF files of surface albedo weights.
# Normally 'file' will have the extension .surf. 'bands' is the number bands.
# 'alb' should then be an array of surface albedo values for each band.

    basis = 1
    n_lon = np.size(lon)
    n_lat = np.size(lat)
    nvals = np.size(alb)

    if (nvals == 1):
        albs = np.zeros(n_lon * n_lat * np.sum(bands))
        albs += alb
    elif (nvals == n_lon * n_lat * bands):
        albs = alb
    else:
        raise RuntimeError(' Error in ncout_surf: arrays dont match', nvals
                           , n_lon * n_lat * bands)

    print('ncout_spectral_surf - file: ', file)

    ncdf_file = create_cdf(file)
    write_dim(ncdf_file, n_lon, lon, 'lon', 'f4', 'lon', 'degree', 'LONGITUDE')
    write_dim(ncdf_file, n_lat, lat, 'lat', 'f4', 'lat', 'degree', 'LATITUDE')
    write_dim(ncdf_file, 1, basis, 'basis', 'i2', 'basis', 'None',
              'BASIS FUNCTION')
    write_dim(ncdf_file, np.sum(bands), np.arange(bands) + 1, 'bands', 'i2', 
              'bands', 'None', 'BANDS')

    write_var(ncdf_file, albs, 'alb', 'f4', ('bands', 'basis', 'lat', 'lon'),
              'None', 'ALBEDO WEIGHTS')

    ncdf_file.close()

#-------------------------------------------------------------------------------

def ncout2d(file, lon, lat, val, name = None, longname = None
            , units = None):
# ; Program to create netCDF files of single level fields.
# ; For example (sol=1365.0 and lon, lat are arrays):
# ; ncout2d, 'out.stoa', lon, lat, sol, longname='Solar Irradiance', units='WM-2'
# ; (the optional argument 'name' is missing here and will be set from
# ; the file extension, i.e. name='stoa').

    n_lon = np.size(lon)
    n_lat = np.size(lat)
    nvals = np.size(val)
    if (nvals == 1):
        vals = np.zeros(n_lon*n_lat)
        vals += np.sum(val)
    elif (nvals == n_lon*n_lat):
        vals = val
    else:
        raise RuntimeError(' Error in ncout2d: arrays dont match', nvals
                           , n_lon * n_lat)

    if (type(name) is not str):
        patterns = re.compile(r'\\|\.')
        name = file[patterns.search(file).start() + 1:]
    print('ncout2d - file: ',file)

    ncdf_file = create_cdf(file)

    write_dim(ncdf_file, n_lon, lon, 'lon', 'f4', 'lon', 'degree', 'LONGITUDE')
    write_dim(ncdf_file, n_lat, lat, 'lat', 'f4', 'lat', 'degree', 'LATITUDE')

    write_var(ncdf_file, vals, name, 'f4', ('lat', 'lon'), units, longname)

    ncdf_file.close()   

#-------------------------------------------------------------------------------

def ncout3d(file, lon, lat, p, val, name = None
            , longname = None, units = None):

# ; Program to create netCDF files of 3d fields on pressure levels.
# ; For example (lon, lat, p, and t are arrays):
# ; ncout3d, 'out.t', lon, lat, p, t, longname='Temperature', units='K'
# ; (the optional argument 'name' is missing here and will be set from
# ; the file extension, i.e. name='t').

    levels = np.size(p)
    n_lon = np.size(lon)
    n_lat = np.size(lat)
    nvals = np.size(val)

    if (nvals == 1):
        vals = np.zeros(n_lon * n_lat * levels).reshape(levels, n_lat, n_lon)
        vals += np.sum(val)
    elif (nvals == levels):
        vals = np.zeros(n_lon * n_lat * levels).reshape(levels, n_lat, n_lon)
        for i in np.arange(levels):
            vals[i, :, :] = val[i]
    elif (nvals == n_lon * n_lat * levels):
        vals = val.reshape(levels, n_lat, n_lon)
    else:
        raise RuntimeError(' Error in ncout3d: arrays dont match', nvals
                           , n_lon * n_lat * levels)

    if (p.size > 1):
        order = np.argsort(p)
        if (any (p[order] != p)):
            p = p[order]
            vals = vals[order,:,:]
    
    if (type(name) is not str):
        patterns = re.compile(r'\\|\.')
        name = file[patterns.search(file).start() + 1:]

    print('ncout3d - file: ',file)

    ncdf_file = create_cdf(file)

    write_dim(ncdf_file, n_lon, lon, 'lon', 'f4', 'lon', 'degree', 'LONGITUDE')
    write_dim(ncdf_file, n_lat, lat, 'lat', 'f4', 'lat', 'degree', 'LATITUDE')
    write_dim(ncdf_file, levels, p, 'plev', 'f4', 'plev', 'Pa', 'PRESSURE')

    title = None
    if (type(units) is not str):
        units = 'None'
    if (type(longname) is str):
        title = longname
    write_var(ncdf_file, vals, name, 'f4', ('plev', 'lat', 'lon'), units, title)

    ncdf_file.close()   

#-------------------------------------------------------------------------------

def ncout_opt_prop(file, lon, lat, p, bands, absp, scat, phf):

# ; Program to create netCDF files of prescribed optical properties 
# ; on pressure levels.
# ; For example (lon, lat, p, absp, scat and phf are arrays):
# ; ncout_opt_prop, 'out.op_soot', lon, lat, p, 6, absp, scat, phf

    levels = np.size(p)
    n_lon = np.size(lon)
    n_lat = np.size(lat)

    nvals = np.size(absp)
    if (nvals == levels*bands):
        absp.reshape(bands, levels)
        absp_vals = np.zeros(n_lon * n_lat * levels * bands).reshape(bands, 
                    levels, n_lat, n_lon)
        for i in np.arange(n_lon):
            for j in np.arange(n_lat):
                absp_vals[:, :, j, i] = absp
    elif (nvals == n_lon*n_lat*levels*bands):
        absp_vals = absp.reshape(bands, levels, n_lat, n_lon)
    else:
        raise RuntimeError(' Error in ncout_opt_prop: absp arrays dont match'
                           , nvals, n_lon*n_lat*levels*bands)

    nvals = np.size(scat)
    if (nvals == levels*bands):
        scat.reshape(bands, levels)
        scat_vals = np.zeros(n_lon * n_lat * levels * bands).reshape(bands,
                    levels, n_lat, n_lon)
        for i in np.arange(n_lon):
            for j in np.arange(n_lat):
                scat_vals[:, :, j, i] = scat
    elif (nvals == n_lon * n_lat * levels * bands):
        scat_vals = scat.reshape(bands, levels, n_lat, n_lon)
    else:
        raise RuntimeError(' Error in ncout_opt_prop: scat arrays dont match'
                           , nvals, n_lon * n_lat * levels * bands)

    nvals = np.size(phf)
    if (nvals == 1):
        phf_vals = np.zeros(n_lon * n_lat * levels * bands).reshape(bands, 1, 
                   levels, n_lat, n_lon)
        phf_vals += np.sum(phf)    
    elif (nvals == levels*bands):
        phf.reshape(levels, bands)
        temp = np.zeros(n_lon * n_lat * levels * bands).reshape(bands, levels,
                        n_lat, n_lon)
        for i in np.arange(n_lon):
            for j in np.arange(n_lat):
                temp[:, :, j, i] = phf
    elif (nvals == n_lon * n_lat * levels * bands):
        phf_vals = phf.reshape(bands, 1, levels, n_lat, n_lon)
    else:
        raise RuntimeError(' Error in ncout_opt_prop: phf arrays dont match'
                           , nvals, n_lon * n_lat * levels * bands)

    order=np.argsort(p)
    if (any (p[order] != p)):
        p = p[order]
        absp_vals = absp_vals[:, order, :, :]
        scat_vals = scat_vals[:, order, :, :]
        phf_vals = phf_vals[:, :, order, :, :]

    print('ncout_opt_prop - file: ', file)

    ncdf_file = create_cdf(file)

    write_dim(ncdf_file, n_lon, lon, 'lon', 'f4', 'lon', 'degree', 'LONGITUDE')
    write_dim(ncdf_file, n_lat, lat, 'lat', 'f4', 'lat', 'degree', 'LATITUDE')
    write_dim(ncdf_file, levels, p, 'plev', 'f4', 'plev', 'Pa', 'PRESSURE')
    write_dim(ncdf_file, 1, 1, 'mom', 'i2', 'mom', 'none', 'moment')
    print('bands=', bands)
    print('np.sum(bands)=', np.sum(bands))
    print('np.arange(bands)+1=', np.arange(bands)+1)
    write_dim(ncdf_file, bands, np.arange(bands)+1, 'band', 'i2', 'band',
              'none', 'band')

    write_var(ncdf_file, absp_vals, 'abs', 'f4', ('band', 'plev', 'lat', 'lon')
              , 'M-1', 'absorption')
    write_var(ncdf_file, scat_vals, 'scat', 'f4', ('band', 'plev', 'lat', 'lon')
              , 'M-1', 'scattering')
    write_var(ncdf_file, phf_vals, 'phf', 'f4', ('band', 'mom', 'plev', 'lat', 
              'lon'), 'none', 'phase function')

    ncdf_file.close()   

#-------------------------------------------------------------------------------

def ncout_view(file, lon, lat, direction, level, pol, azim, rlev):

# Program to create netCDF ".view" files.

    n_lon = np.size(lon)
    n_lat = np.size(lat)
    n_dir  = np.size(direction)
    n_lvl  = np.size(level)
    n_pol  = np.size(pol)
    n_azm  = np.size(azim)
    n_rlv  = np.size(rlev)

    if (n_pol == 1):
        pols = np.zeros(n_lon * n_lat * n_dir).reshape(n_dir, n_lat, n_lon)
        pols += np.sum(pol)
    elif (n_pol == n_lon * n_lat * n_dir):
        pols = pol
    else:
        raise RuntimeError(' Error in ncout_view: arrays dont match'
                           , n_pol, n_lon * n_lat * n_dir)

    if (n_azm == 1):
        azims = np.zeros(n_lon*n_lat*n_dir).reshape(n_dir, n_lat, n_lon)
        azims += np.sum(azim)
    elif (n_azm == n_lon * n_lat * n_dir):
        azims = azim
    else:
        raise RuntimeError(' Error in ncout_view: arrays dont match'
                           , n_azm, n_lon * n_lat * n_dir)

    if (n_rlv == 1):
        rlevs = np.zeros(n_lvl)
        rlevs += np.sum(rlev)
    elif (n_rlv == n_lvl):
        rlevs = rlev
    else:
        raise RuntimeError(' Error in ncout_view: arrays dont match'
                           , n_rlv, n_lvl)

    ncdf_file = create_cdf(file)

    write_dim(ncdf_file, n_lon, lon, 'lon', 'f4', 'lon', 'degree', 'LONGITUDE')
    write_dim(ncdf_file, n_lat, lat, 'lat', 'f4', 'lat', 'degree', 'LATITUDE')
    write_dim(ncdf_file, levels, p, 'level', 'i2', 'level', 'None', 'LEVEL')
    write_dim(ncdf_file, n_dir, direction, 'direction', 'i2', 'direction', 'None', 'MOMENT')

    write_var(ncdf_file, pols, 'pol', 'f4', ('direction', 'lat', 'lon')
              , 'degree', 'POLAR VIEWING ANGLE')
    write_var(ncdf_file, azims, 'azim', 'f4', ('direction', 'lat', 'lon')
              , 'degree', 'AZIMUTHAL VIEWING ANGLE')
    write_var(ncdf_file, rlevs, 'rlev', 'f4', ('level')
              , 'None', 'VIEWING LEVEL')

    ncdf_file.close()   

#-------------------------------------------------------------------------------

def create_cdf(path):
    ncdf_file = Dataset(path, 'w', format = 'NETCDF3_CLASSIC', clobber = 'true')
    return ncdf_file
