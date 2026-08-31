; *****************************COPYRIGHT*******************************
; (C) Crown copyright Met Office. All rights reserved.
; For further details please refer to the file COPYRIGHT.txt
; which you should have received as part of this distribution.
; *****************************COPYRIGHT*******************************
pro circ2nc, arg1

; Converts the text files from http://circ.gsfc.nasa.gov
; into netCDF input files for the radiation code.
; Should be called with the case number as an argument, eg:
; > circ2nc, 2
;
; Make sure you have sourced set_rad_env first.

nctools
spectral_file=['$RAD_DATA/spectra/old/sp_sw_hadgem1_3r'   $
              ,'$RAD_DATA/spectra/ses_sw_jm1_1'    $
              ,'$RAD_DATA/spectra/sp_sw_220_r'    $
              ]
sp_file_name=['hadgem','ses','220']

nspfs=n_elements(spectral_file)
ncase=strtrim(string(arg1),2)
Tsfc_sza_nlev='Tsfc_sza_nlev_case'+ncase+'.txt'
aerosol_input='aerosol_input_case'+ncase+'.txt'
cloud_input='cloud_input_case'+ncase+'.txt'
layer_input='layer_input_case'+ncase+'.txt'
level_input='level_input_case'+ncase+'.txt'
sfcalbedo_input='sfcalbedo_input_case'+ncase+'.txt'
SW_charts='SW_charts_1cm-1_case'+ncase+'.txt'
LW_lblrtm_bb='LW_lblrtm_bb_case'+ncase+'.txt'

a1=''
a2=strarr(2)
a3=strarr(3)
a4=strarr(4)
a5=strarr(5)
a6=strarr(6)

openr,1,Tsfc_sza_nlev
  readf,1,a5
  readf,1,nlevs
  readf,1,tstar
  readf,1,szen
  readf,1,stoa
close,1

print, 'Number of levels: ', nlevs
print, 'Surface temperature: ',tstar
print, 'Solar zenith angle: ',szen
print, 'Solar irradiance: ',stoa
stoa=stoa/cos(szen*!pi/180.)
print, 'Solar irradiance (corrected): ',stoa

aerosol=fltarr(4,nlevs-1)
openr,1,aerosol_input
  readf,1,a3
  readf,1,alpha
  readf,1,a2
  readf,1,aerosol
close,1

cloud=fltarr(6,nlevs-1)
openr,1,cloud_input
  readf,1,a3
  readf,1,cloud
close,1

layer=fltarr(13,nlevs-1)
legend=''
openr,1,layer_input
  readf,1,a2
  readf,1,legend
  readf,1,layer
close,1

level=fltarr(4,nlevs)
openr,1,level_input
  readf,1,a1
  readf,1,level
close,1

albedo=fltarr(4,49180)
openr,1,sfcalbedo_input
  readf,1,a6
  readf,1,albedo
close,1
wvl=albedo(0,*)

charts=fltarr(4,49180)
openr,1,SW_charts
  readf,1,a5
  readf,1,charts
close,1

lwout=fltarr(5,nlevs)
openr,1,LW_lblrtm_bb
  readf,1,a4
  readf,1,lwout
close,1

;-------------------------------------

basename='case'+ncase
lon=[0.0]
lat=[0.0]
pl=fltarr(nlevs)
tl=fltarr(nlevs)
z=fltarr(nlevs) ; height on levels in m

lwflxup=fltarr(nlevs)
lwflxdn=fltarr(nlevs)

lwalb   = 0.0
basis   = 1

pstar=[level(2,0)]
for i=0,nlevs-1 do begin
  z(i)=level(1,i)*1000.0
  pl(i)=level(2,i)*100.0
  tl(i)=level(3,i)
  lwflxup(i)=lwout(2,nlevs-1-i)
  lwflxdn(i)=lwout(3,nlevs-1-i)
endfor

lwHR=fltarr(nlevs-1)

p=fltarr(nlevs-1)
t=fltarr(nlevs-1)
h2o=fltarr(nlevs-1)
co2=fltarr(nlevs-1)
o3=fltarr(nlevs-1)
n2o=fltarr(nlevs-1)
co=fltarr(nlevs-1)
ch4=fltarr(nlevs-1)
o2=fltarr(nlevs-1)
ccl4=fltarr(nlevs-1)
cfc11=fltarr(nlevs-1)
cfc12=fltarr(nlevs-1)

clfr=fltarr(nlevs-1)
lwc=fltarr(nlevs-1)
iwc=fltarr(nlevs-1)
re=fltarr(nlevs-1)
ire=fltarr(nlevs-1)

for i=0,nlevs-2 do begin
  lwHR(i)=lwout(4,nlevs-1-i)
  p(i)=layer(1,i)*100.0
  t(i)=layer(2,i)
  h2o(i)=layer(3,i)*18.0154/28.97
  co2(i)=layer(4,i)*44.0098/28.97
  o3(i)=layer(5,i)*47.9982/28.97
  n2o(i)=layer(6,i)*44.0128/28.97
  co(i)=layer(7,i)*28.0101/28.97
  ch4(i)=layer(8,i)*16.0428/28.97
  o2(i)=layer(9,i)*31.9988/28.97
  ccl4(i)=layer(10,i)*153.8215/28.97
  cfc11(i)=layer(11,i)*137.370/28.97
  cfc12(i)=layer(12,i)*120.910/28.97

  clfr(i)=cloud(1,i)
  lwc(i)=cloud(2,i)/((z(i+1)-z(i))*1000.0)
  iwc(i)=cloud(3,i)/((z(i+1)-z(i))*1000.0)
  
  re(i)=cloud(4,i)*1.0e-6
  ire(i)=cloud(5,i)*1.0e-6
endfor

re(where(re eq 0.0))=6.0e-6
ire(where(ire eq 0.0))=30.0e-6

; Temperature:
ncout3d, basename + '.tstar', lon, lat, pstar, tstar, longname='Surface temperature', units='K'
ncout3d, basename + '.t', lon, lat, p, t, longname='Temperature', units='K'
ncout3d, basename + '.tl', lon, lat, pl, tl, longname='Temperature on levels', units='K'

; Solar beam:
ncout2d, basename + '.stoa', lon, lat, stoa, longname='Solar Irradiance', units='WM-2'
ncout2d, basename + '.szen', lon, lat, szen, longname='Solar Zenith Angle', units='degrees'

; Surface albedos:
; Note: <basename>.surf should be linked to <basename>.surfsw or
; <basename>.surflw  for the SW and LW calls to the radiation code
; respectively.
ncout_surf, basename + '.surflw', lon, lat, basis, lwalb

; Aerosols:
ammr=fltarr(nlevs-1)
for ii=0,nlevs-2 do begin
  ammr(ii)=aerosol(1,ii)/(z(ii+1)-z(ii))
endfor

sel=where((aerosol(1,*) gt 0.0), count)
if (count gt 0) then begin
  omega=aerosol(2,sel(0))
  asym=aerosol(3,sel(0))
endif else begin
  omega=0.0
  asym=0.0
endelse

for i=0,nspfs-1 do begin
  bands=0
  openr,1,spectral_file(i)
  while not eof(1) do begin
    readf, 1, a1
    if strmid(a1,1,23) eq 'umber of spectral bands' or $
       strmid(a1,1,23) eq 'UMBER OF SPECTRAL BANDS' then begin
      bands=fix(strmid(a1,27,6))
      band_limits=fltarr(3,bands)
    endif
    if strmid(a1,1,34) eq 'pecification of spectral intervals' or $
       strmid(a1,1,34) eq 'PECIFICATION OF SPECTRAL INTERVALS' then begin
      readf, 1, a2
      readf,1,band_limits
    endif
  endwhile
  close,1
  
  kabs=fltarr(bands)
  kscat=fltarr(bands)
  absp=fltarr(nlevs-1,bands)
  scat=fltarr(nlevs-1,bands)

  swalb=fltarr(bands)
  wswalb=fltarr(bands)
  w2swalb=fltarr(bands)
  sp_file=sp_file_name(i)
  for band=0,bands-1 do begin
    ll=1.0/(100.0*band_limits(1,band))
    ul=1.0/(100.0*band_limits(2,band))
    sel=where((wvl ge ul) and (wvl le ll))
    swalb(band)=SUM(albedo(1,sel)*albedo(3,sel))/SUM(albedo(3,sel))
;    wswalb(band)=SUM(albedo(2,sel))/SUM(charts(1,sel))
    w2swalb(band)=SUM(albedo(1,sel)*charts(1,sel))/SUM(charts(1,sel))

    kext=SUM((1.0e4/wvl(sel))^(-alpha))/n_elements(sel)
    if i eq 0 then print,band+1,' Tau: ',kext*aerosol(1,0),kext*aerosol(1,1), $
      kext*aerosol(1,2),kext*aerosol(1,3),kext*aerosol(1,4),kext*aerosol(1,5)
    kscat(band)=omega*kext
    kabs(band)=kext-kscat(band)
    absp(*,band)=ammr*kabs(band)
    scat(*,band)=ammr*kscat(band)
  endfor
  ncout_spectral_surf, basename + '.surfsw_' + sp_file, lon, lat, bands, swalb  
;  ncout_spectral_surf, basename + '.surfwsw_' + sp_file, lon, lat, bands, wswalb  
  ncout_spectral_surf, basename + '.surfw2sw_' + sp_file, lon, lat, bands, w2swalb  
  ncout_opt_prop, basename + '.op_soot_' + sp_file, lon, lat, p, bands, absp, scat, asym
endfor

; Spectrally constant albedos:
CASE ncase OF
'1': ncout_surf, basename + '.surfsw_con', lon, lat, 1, 0.196
'2': ncout_surf, basename + '.surfsw_con', lon, lat, 1, 0.188
'3': ncout_surf, basename + '.surfsw_con', lon, lat, 1, 0.171
'4': ncout_surf, basename + '.surfsw_con', lon, lat, 1, 0.670
'5': ncout_surf, basename + '.surfsw_con', lon, lat, 1, 0.670
'6': ncout_surf, basename + '.surfsw_con', lon, lat, 1, 0.136
'7': ncout_surf, basename + '.surfsw_con', lon, lat, 1, 0.164
ENDCASE

; Main Gases:
ncout3d, basename + '.q',   lon, lat, p, h2o, longname='Specific Humidity'
ncout3d, basename + '.co2', lon, lat, p, co2, longname='Carbon Dioxide MMR'
ncout3d, basename + '.ch4', lon, lat, p, ch4, longname='Methane MMR'
ncout3d, basename + '.n2o', lon, lat, p, n2o, longname='Dinitrogen oxide MMR'
ncout3d, basename + '.o2',  lon, lat, p, o2,  longname='Oxygen MMR'
ncout3d, basename + '.o3',  lon, lat, p, o3,  longname='Ozone MMR'

; CFCs:
ncout3d, basename + '.cfc11',   lon, lat, p, cfc11,   longname='CFC-11 MMR'
ncout3d, basename + '.cfc12',   lon, lat, p, cfc12,   longname='CFC-12 MMR'

; Gases currently unused:
ncout3d, basename + '.co',  lon, lat, p, h2o, longname='Carbon Monoxide MMR'
ncout3d, basename + '.ccl4',lon, lat, p, ccl4,longname='CCL4 MMR'

; Cloud fields:
ncout3d, basename + '.clfr',lon, lat, p, clfr,longname='Cloud fraction'
ncout3d, basename + '.lwc', lon, lat, p, lwc, longname='Liquid water content', units='KGM-3'
ncout3d, basename + '.iwc', lon, lat, p, iwc, longname='Ice water content', units='KGM-3'
ncout3d, basename + '.re',  lon, lat, p, re,  longname='Droplet effective radius', units='M'
ncout3d, basename + '.ire', lon, lat, p, ire, longname='Ice effective radius', units='M'

spawn,'Cdentomix -q '+basename+'.q -t '+basename+'.t -o '+basename+'.lwm '+basename+'.lwc'
spawn,'Cdentomix -q '+basename+'.q -t '+basename+'.t -o '+basename+'.iwm '+basename+'.iwc'

ncout3d, basename + '.ccfr',lon, lat, p, 0.0, longname='Convective Cloud fraction'
ncout3d, basename +'.lwmcv',lon, lat, p, 0.0, longname='Convective liquid water content', units='KGM-3'
ncout3d, basename +'.iwmcv',lon, lat, p, 0.0, longname='Convective ice water content', units='KGM-3'
ncout3d, basename + '.recv',lon, lat, p, re,  longname='Droplet effective radius', units='M'
ncout3d, basename +'.irecv',lon, lat, p, ire, longname='Ice effective radius', units='M'

;LW output files:
ncout3d, basename + '_lwref.dflx',lon,lat,pl,lwflxdn,longname='downward flux',units='WM-2'
ncout3d, basename + '_lwref.uflx',lon,lat,pl,lwflxup,longname='upward flux',units='WM-2'
ncout3d, basename + '_lwref.hrts',lon,lat,p,lwHR,longname='heating rates',units='K.D-1'

end
