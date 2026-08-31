; *****************************COPYRIGHT*******************************
; (C) Crown copyright Met Office. All rights reserved.
; For further details please refer to the file COPYRIGHT.txt
; which you should have received as part of this distribution.
; *****************************COPYRIGHT*******************************
PRO ncout_tl, basename

; Program to create netCDF file of temperature on levels (.tl file).
; This is extrapolated from the temperature in layers (.t file),
; and the surface temperature (.tstar) if available.
; Later modified to work with tidl  

nctools
strip_ncdf, basename+'.t', tfile

lon=tfile.val.lon
lat=tfile.val.lat
p=tfile.val.plev
t=tfile.val.t
tdiml=tfile.diml.t
t=reform(t,tdiml)

n_lon = n_elements(lon)
n_lat = n_elements(lat)
layers= n_elements(p)
pl=fltarr(layers+1)
tlev=fltarr(n_lon,n_lat,layers+1)

tdims=tfile.dims.t
londim=where(tdims eq 'lon')
if londim(0) ne 0 then begin
   d=size(t)
   d=d(0)
   i=indgen(d)
   i(0)=londim(0)
   i(londim(0))=0
   t=transpose(t,[i])
   tdiml(londim(0))=tdiml(0)
   tdiml(0)=n_lon
   t=reform(t,tdiml)
endif
latdim=where(tdims eq 'lat')
if latdim(0) ne 1 then begin
   d=size(t)
   d=d(0)
   i=indgen(d)
   i(1)=londim(0)
   i(londim(0))=1
   t=transpose(t,[i])
   PRINT,t
   tdiml(latdim(0))=tdiml(1)
   tdiml(1)=n_lat
   t=reform(t,tdiml)
endif

order=reverse(sort(p))
p=p(order)
t=t(*,*,order)

x=findfile(basename+'.tstar',count=exists)
if exists then begin
    strip_ncdf, basename+'.tstar', tsfile
    tstar=tsfile.val.tstar
    pstar=tsfile.val.plev
endif else begin
    pstar=p(0)
endelse

pl(0)=pstar
pl(1:layers-1)=0.5*(p(0:layers-2)+p(1:layers-1))
pl(layers)=p(layers-1)*pl(layers-1)/p(layers-2)

xs=alog(p)
ts=alog(pl)

FOR i = 0, n_lon-1 DO BEGIN
    FOR j = 0, n_lat-1 DO BEGIN
        ys=reform( t(i,j,*) )
        tlev(i,j,*)=INTERPOL(ys,xs,ts)
    ENDFOR
ENDFOR

if exists then begin
   tlev(*,*,0)=tstar
endif else begin
   ncout3d, basename + '.tstar', lon, lat, pstar, tlev(*,*,0), $
             longname='Surface Temperature', units='K'
endelse

const=287.*250./(9.80665*1000.0)
plot,t(0,0,*),-alog(p/pl(0))*const,xtitle='Temperature (K)', $
 ytitle='Approx height (km)',title='Interpolation of temperature to levels'
oplot,tlev(0,0,*),-alog(pl/pl(0))*const,psym=1

ncout3d, basename + '.tl', lon, lat, pl, tlev, $
          longname='Temperature on levels', units='K'

;print, 'ncout_tl - file: ',basename+'.tl'

end ; ncout_tl
