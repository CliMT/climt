; *****************************COPYRIGHT*******************************
; (C) Crown copyright Met Office. All rights reserved.
; For further details please refer to the file COPYRIGHT.txt
; which you should have received as part of this distribution.
; *****************************COPYRIGHT*******************************
PRO ncprofiles, basename

; Creates many of the files needed by the radiation code using
; standard default values.
; Note: Files will be overwritten - comment out the calls to 
; ncout* that are not required.

ext     = '.t'                         ; Template file extension
                                       
stoa    = 1365.0                       ; Solar Irradiance
szen    = 0.0                          ; Solar Zenith Angle
                                       
lwalb   = 0.0                          ; Surface Albedos
swalb   = 0.06                         ; Ocean~0.06, Desert~0.3, New Snow~0.8
basis   = 1                            ;

co2     = 5.3938e-4 ; CO2 (355 ppmv)   ;
ch4     = 9.637e-7  ; CH4 (1.74 ppmv)  ; Main Gases
n2o     = 4.7255e-7 ; N2O (311 ppbv)   ; (Mass mixing ratios)
o2      = 0.23139   ; O2  (Oxygen)     ;

cfc11   = 0.0                          ;
cfc12   = 0.0                          ;
cfc113  = 0.0                          ; CFCs
hcfc22  = 0.0                          ; (Mass mixing ratios)
hfc125  = 0.0                          ;
hfc134a = 0.0                          ;  
                                       
re      = 10.0e-6                      ; Cloud Droplets Effective Radii
ire     = 30.0e-6                      ; Ice Crystal Effective Radii
recv    = 10.0e-6                      ; Convective ql Effective Radii
irecv   = 30.0e-6                      ; Convective qi Effective Radii

;------------------------------------------------------------------------------

nctools
strip_ncdf, basename+ext, tfile
lon=tfile.val.lon
lat=tfile.val.lat
p=tfile.val.plev

; Solar beam:
ncout2d, basename + '.stoa', lon, lat, stoa, longname='Solar Irradiance', units='WM-2'
ncout2d, basename + '.szen', lon, lat, szen, longname='Solar Zenith Angle', units='degrees'

; Surface albedos:
; Note: <basename>.surf should be linked to <basename>.surfsw or
; <basename>.surflw  for the SW and LW calls to the radiation code
; respectively.
ncout_surf, basename + '.surfsw', lon, lat, basis, swalb
ncout_surf, basename + '.surflw', lon, lat, basis, lwalb

; Main Gases:
ncout3d, basename + '.co2', lon, lat, p, co2, longname='Carbon Dioxide MMR'
ncout3d, basename + '.ch4', lon, lat, p, ch4, longname='Methane MMR'
ncout3d, basename + '.n2o', lon, lat, p, n2o, longname='Dinitrogen oxide MMR'
ncout3d, basename + '.o2',  lon, lat, p, o2,  longname='Oxygen MMR'

; CFCs:
ncout3d, basename + '.cfc11',   lon, lat, p, cfc11,   longname='CFC-11 MMR'
ncout3d, basename + '.cfc12',   lon, lat, p, cfc12,   longname='CFC-12 MMR'
ncout3d, basename + '.cfc113',  lon, lat, p, cfc113,  longname='CFC-113 MMR'
ncout3d, basename + '.hcfc22',  lon, lat, p, hcfc22,  longname='HCFC-22 MMR'
ncout3d, basename + '.hfc125',  lon, lat, p, hfc125,  longname='HFC-125 MMR'
ncout3d, basename + '.hfc134a', lon, lat, p, hfc134a, longname='HFC-134A MMR'

; Cloud particle effective radii:
ncout3d, basename + '.re',  lon, lat, p, re, longname='Cloud Droplets Effective Radii', units='M'
ncout3d, basename + '.ire', lon, lat, p, ire, longname='Ice Crystal Effective Radii', units='M'
ncout3d, basename + '.recv',  lon, lat, p, recv, longname='Convective ql Effective Radii', units='M'
ncout3d, basename + '.irecv', lon, lat, p, irecv, longname='Convective qi Effective Radii', units='M'

; Call ncout_tl to create .tl and .tstar files:
ncout_tl, basename

END
