; *****************************COPYRIGHT*******************************
; (C) Crown copyright Met Office. All rights reserved.
; For further details please refer to the file COPYRIGHT.txt
; which you should have received as part of this distribution.
; *****************************COPYRIGHT*******************************
pro ncplot, file

; Plots a mean profile of the variable in the supplied file
; against height (calculated from pressure assuming isothermal atmos.)

strip_ncdf, file, dgs

parts  = strsplit(file,'\.',/extract)
name = parts(1)

lon=dgs.val.lon
lat=dgs.val.lat
p=dgs.val.plev

n_lon = n_elements(lon)
n_lat = n_elements(lat)
layers= n_elements(p)

cmd="var=dgs.val."+name
err=execute(cmd)
cmd="vardiml=dgs.diml."+name
err=execute(cmd)
var=reform(var,vardiml)

cmd="dims=dgs.dims."+name
err=execute(cmd)
londim=where(dims eq 'lon')
if londim(0) ne 0 then begin
   d=size(var)
   d=d(0)
   i=indgen(d)
   i(0)=londim(0)
   i(londim(0))=0
   var=transpose(var,[i])
   vardiml(londim(0))=vardiml(0)
   vardiml(0)=n_lon
   var=reform(var,vardiml)
endif
latdim=where(dims eq 'lat')
if latdim(0) ne 1 then begin
   d=size(var)
   d=d(0)
   i=indgen(d)
   i(1)=londim(0)
   i(londim(0))=1
   var=transpose(var,[i])
   vardiml(latdim(0))=vardiml(1)
   vardiml(1)=n_lat
   var=reform(var,vardiml)
endif

vmean=fltarr(layers)
for i=0,layers-1 do begin
   vmean(i)=TOTAL(var(*,*,i))/(n_lon*n_lat)
endfor

const=287.*250./(9.80665*1000.0)
plot,vmean,-alog(p/max(p))*const,xtitle=name,ytitle='Approx height (km)'

;print,'Pressure, '+name
;for i=0,n_elements(p)-1 do begin
;   print,p(i),vmean(i)
;endfor

end
