; *****************************COPYRIGHT*******************************
; (C) Crown copyright Met Office. All rights reserved.
; For further details please refer to the file COPYRIGHT.txt
; which you should have received as part of this distribution.
; *****************************COPYRIGHT*******************************

; nctools
;
; This file contains IDL routines that can be used to produce
; multicolumn netCDF data compatible with l_run_cdf.

;-------------------------------------------------------------------------------
PRO ncout_surf, file, lon, lat, basis, alb

; Program to create netCDF files of surface albedo weights.
; Normally 'file' will have the extension .surf, and 'basis' = 1.
; 'alb' should then be an array of surface albedo values or a single value.

FORWARD_FUNCTION create_cdf

levels = n_elements(basis)
n_lon = n_elements(lon)
n_lat = n_elements(lat)
nvals  = n_elements(alb)
IF nvals EQ 1 THEN BEGIN
   albs = REPLICATE(total(alb), n_lon, n_lat, levels)
ENDIF ELSE IF nvals EQ n_lon*n_lat*levels THEN BEGIN
   albs = reform(alb,[n_lon,n_lat,levels])
ENDIF ELSE BEGIN
   print, ' Error in ncout_surf: arrays dont match', nvals, n_lon*n_lat*levels
   retall
ENDELSE

;print, 'ncout_surf - file: ',file

cdfid = create_cdf(file)
lon_dimid = ncdf_dimdef(cdfid, 'lon', n_lon)
lat_dimid = ncdf_dimdef(cdfid, 'lat', n_lat)
basis_dimid = ncdf_dimdef(cdfid, 'basis', levels)

lon_varid  = ncdf_vardef(cdfid, 'lon',  lon_dimid,/FLOAT)
ncdf_attput, cdfid, lon_varid, 'units', 'degree', LENGTH=6, /CHAR
ncdf_attput, cdfid, lon_varid, 'title', 'LONGITUDE', LENGTH=9, /CHAR

lat_varid  = ncdf_vardef(cdfid, 'lat',  lat_dimid,  /FLOAT)
ncdf_attput, cdfid, lat_varid, 'units', 'degree', LENGTH=6, /CHAR
ncdf_attput, cdfid, lat_varid, 'title', 'LATITUDE', LENGTH=8, /CHAR

basis_varid = ncdf_vardef(cdfid, 'basis', basis_dimid, /SHORT)
ncdf_attput,cdfid, basis_varid, 'units', 'None', LENGTH=4, /CHAR
ncdf_attput,cdfid, basis_varid, 'title', 'BASIS FUNCTION', LENGTH=14, /CHAR

alb_varid = ncdf_vardef(cdfid, 'alb', [lon_dimid,lat_dimid,basis_dimid], /FLOAT)
                            
ncdf_attput,cdfid, alb_varid, 'units', 'None', LENGTH=4, /CHAR
ncdf_attput,cdfid, alb_varid, 'title', 'ALBEDO WEIGHTS', LENGTH=14, /CHAR
ncdf_control,cdfid,/endef

ncdf_varput,cdfid, lon_varid, lon 
ncdf_varput,cdfid, lat_varid, lat
ncdf_varput,cdfid, basis_varid, basis
ncdf_varput,cdfid, alb_varid, albs

ncdf_close,cdfid
end ; ncout_surf

;-------------------------------------------------------------------------------
PRO ncout_spectral_surf, file, lon, lat, bands, alb

; Program to create netCDF files of surface albedo weights.
; Normally 'file' will have the extension .surf. 'bands' is the number bands.
; 'alb' should then be an array of surface albedo values for each band.

FORWARD_FUNCTION create_cdf

basis=1
n_lon = n_elements(lon)
n_lat = n_elements(lat)
nvals  = n_elements(alb)
IF nvals EQ 1 THEN BEGIN
   albs = REPLICATE(total(alb), n_lon, n_lat, 1, total(bands))
ENDIF ELSE IF nvals EQ n_lon*n_lat*bands THEN BEGIN
   albs = reform(alb,[n_lon,n_lat,1,total(bands)])
ENDIF ELSE BEGIN
   print, ' Error in ncout_surf: arrays dont match', nvals, n_lon*n_lat*bands
   retall
ENDELSE

;print, 'ncout_spectral_surf - file: ',file

cdfid = create_cdf(file)
lon_dimid = ncdf_dimdef(cdfid, 'lon', n_lon)
lat_dimid = ncdf_dimdef(cdfid, 'lat', n_lat)
basis_dimid = ncdf_dimdef(cdfid, 'basis', 1)
bands_dimid = ncdf_dimdef(cdfid, 'bands', total(bands))

lon_varid  = ncdf_vardef(cdfid, 'lon',  lon_dimid,/FLOAT)
ncdf_attput, cdfid, lon_varid, 'units', 'degree', LENGTH=6, /CHAR
ncdf_attput, cdfid, lon_varid, 'title', 'LONGITUDE', LENGTH=9, /CHAR

lat_varid  = ncdf_vardef(cdfid, 'lat',  lat_dimid,  /FLOAT)
ncdf_attput, cdfid, lat_varid, 'units', 'degree', LENGTH=6, /CHAR
ncdf_attput, cdfid, lat_varid, 'title', 'LATITUDE', LENGTH=8, /CHAR

basis_varid = ncdf_vardef(cdfid, 'basis', basis_dimid, /SHORT)
ncdf_attput,cdfid, basis_varid, 'units', 'None', LENGTH=4, /CHAR
ncdf_attput,cdfid, basis_varid, 'title', 'BASIS FUNCTION', LENGTH=14, /CHAR

bands_varid = ncdf_vardef(cdfid, 'bands', bands_dimid, /SHORT)
ncdf_attput,cdfid, bands_varid, 'units', 'None', LENGTH=4, /CHAR
ncdf_attput,cdfid, bands_varid, 'title', 'BANDS', LENGTH=5, /CHAR

alb_varid = ncdf_vardef(cdfid, 'alb',$
            [lon_dimid,lat_dimid,basis_dimid,bands_dimid], /FLOAT)
ncdf_attput,cdfid, alb_varid, 'units', 'None', LENGTH=4, /CHAR
ncdf_attput,cdfid, alb_varid, 'title', 'ALBEDO WEIGHTS', LENGTH=14, /CHAR

ncdf_control,cdfid,/endef

ncdf_varput,cdfid, lon_varid, lon 
ncdf_varput,cdfid, lat_varid, lat
ncdf_varput,cdfid, basis_varid, basis
ncdf_varput,cdfid, bands_varid, indgen(bands)+1
ncdf_varput,cdfid, alb_varid, albs

ncdf_close,cdfid
end ; ncout_spectral_surf

;-------------------------------------------------------------------------------
PRO ncout2d, file, lon, lat, val, $
           name=name, longname=longname, units=units

; Program to create netCDF files of single level fields.
; For example (sol=1365.0 and lon, lat are arrays):
; ncout2d, 'out.stoa', lon, lat, sol, longname='Solar Irradiance', units='WM-2'
; (the optional argument 'name' is missing here and will be set from
; the file extension, i.e. name='stoa').


FORWARD_FUNCTION create_cdf

n_lon = n_elements(lon)
n_lat = n_elements(lat)
nvals  = n_elements(val)
IF nvals EQ 1 THEN BEGIN
   vals = REPLICATE(total(val), n_lon, n_lat)
ENDIF ELSE IF nvals EQ n_lon*n_lat THEN BEGIN
   vals = reform(val,[n_lon,n_lat])
ENDIF ELSE BEGIN
   print, ' Error in ncout2d: arrays dont match', nvals, n_lon*n_lat
   retall
ENDELSE

IF size(name, /type) NE 7 THEN BEGIN
  parts  = strsplit(file,'\.',/extract)
  name = parts(1)
ENDIF

;print, 'ncout2d - file: ',file

cdfid = create_cdf(file)
lon_dimid = ncdf_dimdef(cdfid, 'lon', n_lon)
lat_dimid = ncdf_dimdef(cdfid, 'lat', n_lat)

lon_varid  = ncdf_vardef(cdfid, 'lon', lon_dimid,  /FLOAT)
ncdf_attput,cdfid, lon_varid, 'units', 'degree', LENGTH=6, /CHAR
ncdf_attput,cdfid, lon_varid, 'title', 'LONGITUDE', LENGTH=9, /CHAR

lat_varid  = ncdf_vardef(cdfid, 'lat', lat_dimid,  /FLOAT)
ncdf_attput,cdfid, lat_varid, 'units', 'degree', LENGTH=6, /CHAR
ncdf_attput,cdfid, lat_varid, 'title', 'LATITUDE', LENGTH=8, /CHAR

val_varid  = ncdf_vardef(cdfid, name, [lon_dimid,lat_dimid], /FLOAT)
IF size(units, /type) EQ 7 THEN $
  ncdf_attput,cdfid, val_varid, 'units', units, LENGTH=strlen(units), /CHAR
IF size(longname, /type) EQ 7 THEN $
  ncdf_attput,cdfid, val_varid, 'title', longname, LENGTH=strlen(longname), /CHAR
                            
ncdf_control,cdfid,/endef

ncdf_varput,cdfid, lon_varid, lon
ncdf_varput,cdfid, lat_varid, lat
ncdf_varput,cdfid, val_varid, vals

ncdf_close,cdfid

end ; ncout2d

;-------------------------------------------------------------------------------
PRO ncout3d, file, lon, lat, p, val, $
           name=name, longname=longname, units=units

; Program to create netCDF files of 3d fields on pressure levels.
; For example (lon, lat, p, and t are arrays):
; ncout3d, 'out.t', lon, lat, p, t, longname='Temperature', units='K'
; (the optional argument 'name' is missing here and will be set from
; the file extension, i.e. name='t').

FORWARD_FUNCTION create_cdf

levels = n_elements(p)
n_lon = n_elements(lon)
n_lat = n_elements(lat)
nvals  = n_elements(val)
IF nvals EQ 1 THEN BEGIN
   vals = REPLICATE(total(val), n_lon, n_lat, levels)
ENDIF ELSE IF nvals EQ levels THEN BEGIN
   val = reform(val,levels)
   vals=fltarr(n_lon, n_lat, levels)
   FOR i=0,levels-1 DO BEGIN
     FOR j=0,n_lat-1 DO BEGIN
       REPLICATE_INPLACE,vals,val(i),1,[0,j,i]
     ENDFOR
   ENDFOR
ENDIF ELSE IF nvals EQ n_lon*n_lat*levels THEN BEGIN
   vals = reform(val,[n_lon,n_lat,levels])
ENDIF ELSE BEGIN
   print, ' Error in ncout3d: arrays dont match', nvals, n_lon*n_lat*levels
   retall
ENDELSE

plev=[p]
order=SORT(plev)
IF TOTAL(plev(order) NE plev) THEN BEGIN
   plev = plev(order)
   vals = vals(*,*,order)
ENDIF

IF size(name, /type) NE 7 THEN BEGIN
  parts  = strsplit(file,'\.',/extract)
  name = parts(1)
ENDIF

;print, 'ncout3d - file: ',file

cdfid = create_cdf(file)
lon_dimid = ncdf_dimdef(cdfid, 'lon', n_lon)
lat_dimid = ncdf_dimdef(cdfid, 'lat', n_lat)
plev_dimid = ncdf_dimdef(cdfid, 'plev', levels)

lon_varid  = ncdf_vardef(cdfid, 'lon', lon_dimid,  /FLOAT)
ncdf_attput,cdfid, lon_varid, 'units', 'degree', LENGTH=6, /CHAR
ncdf_attput,cdfid, lon_varid, 'title', 'LONGITUDE', LENGTH=9, /CHAR

lat_varid  = ncdf_vardef(cdfid, 'lat', lat_dimid,  /FLOAT)
ncdf_attput,cdfid, lat_varid, 'units', 'degree', LENGTH=6, /CHAR
ncdf_attput,cdfid, lat_varid, 'title', 'LATITUDE', LENGTH=8, /CHAR

plev_varid = ncdf_vardef(cdfid, 'plev', plev_dimid, /FLOAT)
ncdf_attput,cdfid, plev_varid, 'units', 'Pa', LENGTH=2, /CHAR
ncdf_attput,cdfid, plev_varid, 'title', 'PRESSURE', LENGTH=8, /CHAR

val_varid  = ncdf_vardef(cdfid, name, [lon_dimid,lat_dimid,plev_dimid], /FLOAT)

IF size(units, /type) NE 7 THEN units='None'
ncdf_attput,cdfid, val_varid, 'units', units, LENGTH=strlen(units), /CHAR

IF size(longname, /type) EQ 7 THEN $
  ncdf_attput,cdfid, val_varid, 'title', longname, LENGTH=strlen(longname), /CHAR
ncdf_control,cdfid,/endef

ncdf_varput,cdfid, lon_varid, lon
ncdf_varput,cdfid, lat_varid, lat
ncdf_varput,cdfid, plev_varid, plev
ncdf_varput,cdfid, val_varid, vals

ncdf_close,cdfid

end ; ncout3d

;-------------------------------------------------------------------------------
PRO ncout_opt_prop, file, lon, lat, p, bands, absp, scat, phf

; Program to create netCDF files of prescribed optical properties 
; on pressure levels.
; For example (lon, lat, p, absp, scat and phf are arrays):
; ncout_opt_prop, 'out.op_soot', lon, lat, p, 6, absp, scat, phf

FORWARD_FUNCTION create_cdf

levels = n_elements(p)
n_lon = n_elements(lon)
n_lat = n_elements(lat)
nvals  = n_elements(absp)
IF nvals EQ levels*bands THEN BEGIN
   absp = reform(absp,[levels,bands])
   abs_vals = fltarr(n_lon,n_lat,levels,bands)
   FOR j=0,n_lat-1 DO BEGIN
     FOR i=0,n_lon-1 DO BEGIN
       abs_vals(i,j,*,*)=absp(*,*)
     ENDFOR
   ENDFOR
ENDIF ELSE IF nvals EQ n_lon*n_lat*levels*bands THEN BEGIN
   abs_vals = reform(absp,[n_lon,n_lat,levels,bands])
ENDIF ELSE BEGIN
   print, ' Error in ncout_opt_prop: absp arrays dont match', $
     nvals,n_lon*n_lat*levels*bands
   retall
ENDELSE
nvals  = n_elements(scat)
IF nvals EQ levels*bands THEN BEGIN
   scat = reform(scat,[levels,bands])
   scat_vals = fltarr(n_lon,n_lat,levels,bands)
   FOR j=0,n_lat-1 DO BEGIN
     FOR i=0,n_lon-1 DO BEGIN
       scat_vals(i,j,*,*)=scat(*,*)
     ENDFOR
   ENDFOR
ENDIF ELSE IF nvals EQ n_lon*n_lat*levels*bands THEN BEGIN
   scat_vals = reform(scat,[n_lon,n_lat,levels,bands])
ENDIF ELSE BEGIN
   print, ' Error in ncout_opt_prop: scat arrays dont match', $
     nvals,n_lon*n_lat*levels*bands
   retall
ENDELSE
nvals  = n_elements(phf)
IF nvals EQ 1 THEN BEGIN
   phf_vals = REPLICATE(total(phf), n_lon, n_lat, levels, 1, bands)
ENDIF ELSE IF nvals EQ levels*bands THEN BEGIN
   phf = reform(phf,[levels,bands])
   temp = fltarr(n_lon,n_lat,levels,bands)
   FOR j=0,n_lat-1 DO BEGIN
     FOR i=0,n_lon-1 DO BEGIN
       temp(i,j,*,*)=phf(*,*)
     ENDFOR
   ENDFOR
   phf_vals = reform(temp,[n_lon,n_lat,levels,1,bands])
ENDIF ELSE IF nvals EQ n_lon*n_lat*levels*bands THEN BEGIN
   phf_vals = reform(phf,[n_lon,n_lat,levels,1,bands])
ENDIF ELSE BEGIN
   print, ' Error in ncout_opt_prop: phf arrays dont match', $
     nvals,n_lon*n_lat*levels*bands
   retall
ENDELSE

plev=[p]
order=SORT(plev)
IF TOTAL(plev(order) NE plev) THEN BEGIN
   plev = plev(order)
   abs_vals = abs_vals(*,*,order,*)
   scat_vals = scat_vals(*,*,order,*)
   phf_vals = phf_vals(*,*,order,*,*)
ENDIF

;print, 'ncout_opt_prop - file: ',file

cdfid = create_cdf(file)
lon_dimid = ncdf_dimdef(cdfid, 'lon', n_lon)
lat_dimid = ncdf_dimdef(cdfid, 'lat', n_lat)
plev_dimid = ncdf_dimdef(cdfid, 'plev', levels)
mom_dimid = ncdf_dimdef(cdfid, 'mom', 1)
band_dimid = ncdf_dimdef(cdfid, 'band', bands)

lon_varid  = ncdf_vardef(cdfid, 'lon', lon_dimid, /FLOAT)
ncdf_attput,cdfid, lon_varid, 'units', 'degree', LENGTH=6, /CHAR
ncdf_attput,cdfid, lon_varid, 'title', 'longitude', LENGTH=9, /CHAR

lat_varid  = ncdf_vardef(cdfid, 'lat', lat_dimid, /FLOAT)
ncdf_attput,cdfid, lat_varid, 'units', 'degree', LENGTH=6, /CHAR
ncdf_attput,cdfid, lat_varid, 'title', 'latitude', LENGTH=8, /CHAR

plev_varid = ncdf_vardef(cdfid, 'plev', plev_dimid, /FLOAT)
ncdf_attput,cdfid, plev_varid, 'units', 'Pa', LENGTH=2, /CHAR
ncdf_attput,cdfid, plev_varid, 'title', 'pressure', LENGTH=8, /CHAR

mom_varid = ncdf_vardef(cdfid, 'mom', mom_dimid, /SHORT)
ncdf_attput,cdfid, mom_varid, 'units', 'none', LENGTH=4, /CHAR
ncdf_attput,cdfid, mom_varid, 'title', 'moment', LENGTH=6, /CHAR

band_varid = ncdf_vardef(cdfid, 'band', band_dimid, /SHORT)
ncdf_attput,cdfid, band_varid, 'units', 'none', LENGTH=4, /CHAR
ncdf_attput,cdfid, band_varid, 'title', 'band', LENGTH=4, /CHAR

abs_varid  = ncdf_vardef(cdfid, 'abs', $
               [lon_dimid,lat_dimid,plev_dimid,band_dimid], /FLOAT)
ncdf_attput,cdfid, abs_varid, 'units', 'M-1', LENGTH=3, /CHAR
ncdf_attput,cdfid, abs_varid, 'title', 'absorption', LENGTH=10, /CHAR

scat_varid  = ncdf_vardef(cdfid, 'scat', $
               [lon_dimid,lat_dimid,plev_dimid,band_dimid], /FLOAT)
ncdf_attput,cdfid, scat_varid, 'units', 'M-1', LENGTH=3, /CHAR
ncdf_attput,cdfid, scat_varid, 'title', 'scattering', LENGTH=10, /CHAR

phf_varid  = ncdf_vardef(cdfid, 'phf', $
               [lon_dimid,lat_dimid,plev_dimid,mom_dimid,band_dimid], /FLOAT)
ncdf_attput,cdfid, phf_varid, 'units', 'none', LENGTH=4, /CHAR
ncdf_attput,cdfid, phf_varid, 'title', 'phase function', LENGTH=14, /CHAR
ncdf_control,cdfid,/endef

ncdf_varput,cdfid, lon_varid,  lon
ncdf_varput,cdfid, lat_varid,  lat
ncdf_varput,cdfid, plev_varid, plev
ncdf_varput,cdfid, mom_varid,  1
ncdf_varput,cdfid, band_varid, indgen(bands)+1

ncdf_varput,cdfid, abs_varid,  abs_vals
ncdf_varput,cdfid, scat_varid, scat_vals
ncdf_varput,cdfid, phf_varid,  phf_vals

ncdf_close,cdfid

end ; ncout_opt_prop

;-------------------------------------------------------------------------------
PRO ncout_view, file, lon, lat, direction, level, pol, azim, rlev

; Program to create netCDF ".view" files.

FORWARD_FUNCTION create_cdf

n_lon = n_elements(lon)
n_lat = n_elements(lat)
n_dir  = n_elements(direction)
n_lvl  = n_elements(level)
n_pol  = n_elements(pol)
n_azm  = n_elements(azim)
n_rlv  = n_elements(rlev)

IF n_pol EQ 1 THEN BEGIN
   pols = REPLICATE(total(pol),n_lon,n_lat,n_dir)
ENDIF ELSE IF n_pol EQ n_lon*n_lat*n_dir THEN BEGIN
   pols = reform(pol,[n_lon,n_lat,n_dir])
ENDIF ELSE BEGIN
   print, ' Error in ncout_view: arrays dont match', n_pol, n_lon*n_lat*n_dir
   retall
ENDELSE

IF n_azm EQ 1 THEN BEGIN
   azims = REPLICATE(total(azim),n_lon,n_lat,n_dir)
ENDIF ELSE IF n_pol EQ n_lon*n_lat*n_dir THEN BEGIN
   azims = reform(azim,[n_lon,n_lat,n_dir])
ENDIF ELSE BEGIN
   print, ' Error in ncout_view: arrays dont match', n_azm, n_lon*n_lat*n_dir
   retall
ENDELSE

IF n_rlv EQ 1 THEN BEGIN
   rlevs = REPLICATE(total(rlev),n_lvl)
ENDIF ELSE IF n_rlv EQ n_lvl THEN BEGIN
   rlevs = rlev
ENDIF ELSE BEGIN
   print, ' Error in ncout_view: arrays dont match', n_rlv, n_lvl
   retall
ENDELSE

cdfid = create_cdf(file)
lon_dimid       = ncdf_dimdef(cdfid, 'lon', n_lon)
lat_dimid       = ncdf_dimdef(cdfid, 'lat', n_lat)
direction_dimid = ncdf_dimdef(cdfid, 'direction', n_dir)
level_dimid     = ncdf_dimdef(cdfid, 'level', n_lvl)

lon_varid  = ncdf_vardef(cdfid, 'lon',  lon_dimid,/FLOAT)
ncdf_attput, cdfid, lon_varid, 'units', 'degree', LENGTH=6, /CHAR
ncdf_attput, cdfid, lon_varid, 'title', 'LONGITUDE', LENGTH=9, /CHAR

lat_varid  = ncdf_vardef(cdfid, 'lat',  lat_dimid,  /FLOAT)
ncdf_attput, cdfid, lat_varid, 'units', 'degree', LENGTH=6, /CHAR
ncdf_attput, cdfid, lat_varid, 'title', 'LATITUDE', LENGTH=8, /CHAR

direction_varid = ncdf_vardef(cdfid, 'direction', direction_dimid, /SHORT)
ncdf_attput,cdfid, direction_varid, 'units', 'None', LENGTH=4, /CHAR
ncdf_attput,cdfid, direction_varid, 'title', 'MOMENT', LENGTH=6, /CHAR

level_varid = ncdf_vardef(cdfid, 'level', level_dimid, /SHORT)
ncdf_attput,cdfid, level_varid, 'units', 'None', LENGTH=4, /CHAR
ncdf_attput,cdfid, level_varid, 'title', 'LEVEL', LENGTH=5, /CHAR

pol_varid = ncdf_vardef(cdfid, 'pol',[lon_dimid,lat_dimid,direction_dimid], /FLOAT)
ncdf_attput,cdfid, pol_varid, 'units', 'degree', LENGTH=6, /CHAR
ncdf_attput,cdfid, pol_varid, 'title', 'POLAR VIEWING ANGLE', LENGTH=19, /CHAR

azim_varid = ncdf_vardef(cdfid, 'azim',[lon_dimid,lat_dimid,direction_dimid], /FLOAT)
ncdf_attput,cdfid, azim_varid, 'units', 'degree', LENGTH=6, /CHAR
ncdf_attput,cdfid, azim_varid, 'title', 'AZIMUTHAL VIEWING ANGLE', LENGTH=23, /CHAR

rlev_varid = ncdf_vardef(cdfid, 'rlev', level_dimid, /FLOAT)
ncdf_attput,cdfid, rlev_varid, 'units', 'None', LENGTH=4, /CHAR
ncdf_attput,cdfid, rlev_varid, 'title', 'VIEWING LEVEL', LENGTH=13, /CHAR

ncdf_control,cdfid,/endef

ncdf_varput,cdfid, lon_varid, lon 
ncdf_varput,cdfid, lat_varid, lat
ncdf_varput,cdfid, direction_varid, direction
ncdf_varput,cdfid, level_varid, level 
ncdf_varput,cdfid, pol_varid, pols
ncdf_varput,cdfid, azim_varid, azims
ncdf_varput,cdfid, rlev_varid, rlevs

ncdf_close,cdfid
end ; ncout_view

;-------------------------------------------------------------------------------
function create_cdf, path

command = "echo 'netcdf " + path + "{}' | ncgen -o " + path
spawn, command
cdfid = ncdf_open(path, /WRITE)
IF cdfid LT 0 THEN message, 'ncdf_open failed: '+path
NCDF_CONTROL,cdfid,/redef
return, cdfid

end

;-------------------------------------------------------------------------------
PRO nctools
common share1, started
if not keyword_set(started) then begin
  started='yes'
endif
END
