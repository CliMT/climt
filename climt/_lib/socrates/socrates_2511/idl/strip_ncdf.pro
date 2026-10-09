; *****************************COPYRIGHT*******************************
; (C) Crown copyright Met Office. All rights reserved.
; For further details please refer to the file COPYRIGHT.txt
; which you should have received as part of this distribution.
; *****************************COPYRIGHT*******************************
;
; Procedure to strip the diagnostics from a netCDF file and place
; the data in a structure.
;
; Input:  file - the name of the netCDF file
; Output: data  - the structure containing the stripped data
; 
; Call as
;        IDL> strip_ncdf, '~/path/to/my/file.nc', data
; The data structure contains 3 substructure; val containing the
; values of the data, diml containing the dimensions of the data and
; dims containing the names of the dimensions. 
; Then to use the data:
;        IDL> x=data.val.my_variable
; x now contains the data from "my_variable".
;
;;=================================================================
PRO strip_ncdf, file, data


;---------------------
; Initialization stuff
;---------------------

message, /informational, "processing "+file

;--------------------------
; open the netcdf file
;--------------------------
; There is a problem expanding shell variables ($HOME,~,...)
; when using ncid, so here's a quick work around...
spawn, /sh, 'echo '+file, spawnout
 file=spawnout(0)
ncid=ncdf_open(file,/nowrite)

;--------------------------
; Inquire file
;--------------------------
result=ncdf_inquire(ncid)

;--------------------------
; Place the desired 
; variables in local arrays.
;--------------------------

names=strarr(result.nvars)
for ivar=0, result.Nvars-1 do begin
  vardata = NCDF_VARINQ(ncid, ivar)
  valid_varname = vardata.Name
  varname=valid_varname
  ndims=vardata.ndims>1
  names(ivar)=vardata.Name
;--------------------------
;get the dimensions
;--------------------------
  dimnames=strarr(ndims)
  dimlengths=lonarr(ndims)
  FOR idim=0,ndims-1 DO BEGIN
    ncdf_diminq,ncid,vardata.dim(idim),dimname,length
    IF dimname EQ 'TIME' THEN dimname='TIMES' ; Rather foolishly, the LEM
                                ; outputs the dimension as TIME but
                                ; the variable as TIMES.
    dimnames(idim)=dimname
    dimlengths(idim)=length
  ENDFOR
  NCDF_VARGET, ncid, ivar, varname      
  if ivar eq 0 then begin
    val = create_struct(valid_varname,varname)
    dims = create_struct(valid_varname,dimnames)
    diml = create_struct(valid_varname,dimlengths)
  endif else begin
    val = create_struct(val, valid_varname,varname)
    dims = create_struct(dims,valid_varname,dimnames)
    diml = create_struct(diml,valid_varname,dimlengths)
  endelse
endfor
data = create_struct(['val','diml','dims'],val,diml,dims)

;--------------------------
;Close the netcdf file
;--------------------------

NCDF_CLOSE,ncid

end
