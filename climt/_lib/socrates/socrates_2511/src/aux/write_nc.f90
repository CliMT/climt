! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine to write to a netCDF file.
! Its intended application is to write additional fields 
! from an offline run to a diagnostic file produced in an LFRic run.
!
! Method: In a call to this routine a single field (typically a flux)
!         for a single timestep is added to the netCDF file
!
! ---------------------------------------------------------------------

SUBROUTINE write_nc(filename,var_name,i_time,dataout)

USE realtype_rd, ONLY : RealExt

USE netcdf

IMPLICIT NONE

CHARACTER(LEN=*),INTENT(in) :: filename
CHARACTER(LEN=*),INTENT(in) :: var_name
INTEGER,INTENT(in) :: i_time
REAL(RealExt),ALLOCATABLE,INTENT(in) :: dataout(:)

REAL(RealExt),ALLOCATABLE :: dataint1(:)
REAL(RealExt),ALLOCATABLE :: dataint2(:,:)
REAL(RealExt),ALLOCATABLE :: dataint3(:,:,:)

INTEGER :: ierr
INTEGER :: ierr_var_exist
INTEGER :: ncid
INTEGER :: varid
INTEGER :: varid_out
INTEGER :: ndims
INTEGER,ALLOCATABLE :: dimids(:)
CHARACTER(LEN=200), ALLOCATABLE :: dimension_name(:)
INTEGER,ALLOCATABLE :: dimension_size(:)
INTEGER :: vartype
INTEGER :: total_dims
INTEGER :: total_size
INTEGER :: i
INTEGER :: j
INTEGER :: k

! Open the file for reading and writing
ierr=NF90_OPEN(Trim(filename),NF90_WRITE,ncid)
IF (ierr.NE.0) STOP

ierr_var_exist=NF90_INQ_VARID(ncid, Trim(var_name), varid)

IF (ierr_var_exist.EQ.0) THEN

   ierr=NF90_INQUIRE_VARIABLE(ncid,varid,ndims=ndims,xtype=vartype)
   IF (ierr.NE.0) STOP

   ALLOCATE(dimids(ndims))
   ALLOCATE(dimension_name(ndims))
   ALLOCATE(dimension_size(ndims))
   ierr=NF90_INQUIRE_VARIABLE(ncid,varid,dimids=dimids)
   IF (ierr.NE.0) STOP

   DO i=1,ndims
      ierr=NF90_INQUIRE_DIMENSION(ncid,dimids(i),dimension_name(i),dimension_size(i))
      IF (ierr.NE.0) STOP
   END DO

   total_size=product(dimension_size)
   total_dims=ndims

   DO i=1,ndims
      IF (trim(dimension_name(i))=='time') THEN
         total_size=total_size/dimension_size(i)
         total_dims=total_dims-1
      END IF
   END DO


   ierr=NF90_DEF_VAR(ncid, var_name//'_offline', NF90_FLOAT, dimids(1:ndims), varid_out)
   IF (ierr.NE.0) STOP
   
   ierr=NF90_ENDDEF(ncid)
   IF (ierr.NE.0) STOP
   
   SELECT CASE (total_dims)

      CASE (1)

         ALLOCATE(dataint1(dimension_size(1)))

         DO j=1,dimension_size(1)

            dataint1(j)=dataout(j)

         END DO ! dimension_size(1)

         ierr=NF90_PUT_VAR(ncid,varid_out,dataint1,start=(/1,i_time/))
         IF (ierr.NE.0) STOP

         DEALLOCATE(dataint1)

      CASE (2)

         ALLOCATE(dataint2(dimension_size(1),dimension_size(2)))

         DO i=1,dimension_size(1)
            DO j=1,dimension_size(2)

               dataint2(i,j)=dataout((i-1)*dimension_size(2)+j)

            END DO ! dimension_size(2)
         END DO ! dimension_size(1)

         ierr=NF90_PUT_VAR(ncid,varid_out,dataint2,start=(/1,1,i_time/))
         IF (ierr.NE.0) STOP

         DEALLOCATE(dataint2)

      CASE (3)

         ALLOCATE(dataint3(dimension_size(1),dimension_size(2),dimension_size(3)))

         DO i=1,dimension_size(1)
            DO j=1,dimension_size(2)
               DO k=1,dimension_size(3)

                  dataint3(i,j,k)=dataout(((i-1)*dimension_size(2)+(j-1))*dimension_size(3)+k) 

               END DO ! dimension_size(3)
            END DO ! dimension_size(2)
         END DO ! dimension_size(1)

         ierr=NF90_PUT_VAR(ncid,varid_out,dataint3,start=(/1,1,1,i_time/))
         IF (ierr.NE.0) STOP

         DEALLOCATE(dataint3)

      CASE DEFAULT

      WRITE(*,*) 'treatment of field with dimension ' ,total_dims, '+1 has not been implemented'

   END SELECT

END IF

! Close the file
ierr=NF90_CLOSE(ncid)
IF (ierr.NE.0) STOP

END SUBROUTINE write_nc
