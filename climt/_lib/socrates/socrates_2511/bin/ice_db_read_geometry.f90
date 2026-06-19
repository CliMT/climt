! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read the geometry of ice crystals for a database
!
SUBROUTINE ice_db_read_geometry &
!
(ierr, DBGeom) 
!
! Method:
!   An auxiliary file is read to get the relation between the
!   crystal size and its volume and projected area. This is
!   rather hard-wired at the moment and some generalziation
!   may be required.
!
!
!
! Modules used
  USE realtype_rd
  USE def_std_io_icf
  USE error_pcf
  USE def_db_crystal_geometry
  USE db_type_ucf
  USE spline_fit_mod, ONLY: spline_fit
!
!
  IMPLICIT NONE
!
!
!
! Dummy arguments.
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
  TYPE (STR_db_cryst_geom), Intent(OUT) :: DBGeom
!   Geometry for specifying crystals in a database

!
! Local arguments.
! 
  CHARACTER(LEN=80) :: file_geom_info
!   Filename of geometrical information consistent with the
!   database
  INTEGER :: iu_geom_input
!   Unit number for I/O from the geometrical file
  INTEGER  ::  i
!   Loop variable
  REAL, Allocatable, Dimension(:, :) :: geom_data_temp
!   Temporary array to read in data from 
!   rough aggreagte geometry info file
!   (Note: various files of geometries have been produced
!   with different number of columns and different orders
!   of columns: care is required to determine which to use.
!   This program is written on the assumption that the
!   ordering is:-
!     1. Maximum dimension
!     2. Projected area
!     3. Volume
!     4. Projected area equivalent radius
!     5. Volume equivalent radius)
!
!
!
! Open the rough aggregate geometry file
  CALL get_free_unit(ierr, iu_geom_input)
  IF (ierr /= I_NORMAL) STOP
  WRITE(iu_stdout, '(A)') &
    'Enter the name of the geometrical file.'
  READ(iu_stdin, '(A)') file_geom_info
  OPEN(UNIT=iu_geom_input, FILE=file_geom_info, STATUS='unknown')
!
! Read details from the geometrical file. A hard-wired size is
! used at the moment.
!
  DBGeom%n_geom=24
!
  ALLOCATE(geom_data_temp(24,5))
!
  DO i=1, DBGeom%n_geom
    READ(iu_geom_input,*) geom_data_temp(i,1:5)
  ENDDO
!
! Set up splined fits to the projected area and volume in the
! square and cube of the maximum dimension respectively. For
! irregular particles we do not expect precise linearity in 
! these powers, but they are a reasonable starting point.
  ALLOCATE(DBGeom%dm2(DBGeom%n_geom))
  ALLOCATE(DBGeom%proj_area(DBGeom%n_geom))
  ALLOCATE(DBGeom%d2_proj_area(DBGeom%n_geom))
  ALLOCATE(DBGeom%dm3(DBGeom%n_geom))
  ALLOCATE(DBGeom%volume(DBGeom%n_geom))
  ALLOCATE(DBGeom%d2_volume(DBGeom%n_geom))
!
! Transfer from the temporary array, converting to SI units.
  DO i=1,DBGeom%n_geom
    DBGeom%dm2(i)       = 1.0E-12_RealK * &
      REAL(geom_data_temp(i,1), RealK)**2
    DBGeom%proj_area(i) = 1.0E-12_RealK * &
      REAL(geom_data_temp(i,2), RealK)
    DBGeom%dm3(i)       = 1.0E-18_RealK * &
      REAL(geom_data_temp(i,1), RealK)**3
    DBGeom%volume(i)    = 1.0E-18_RealK * &
      REAL(geom_data_temp(i,3), RealK)
  ENDDO
!
  DEALLOCATE(geom_data_temp)
  CLOSE(iu_geom_input)
!
! Calculate second derivatives to carry out a splined fit.
  ALLOCATE(DBGeom%d2_proj_area(DBGeom%n_geom))
  ALLOCATE(DBGeom%d2_volume(DBGeom%n_geom))
  DBGeom%d2_proj_area = 0.0
  DBGeom%d2_volume    = 0.0
  CALL spline_fit(DBGeom%n_geom, DBGeom%dm2, DBGeom%proj_area, &
    DBGeom%d2_proj_area)
  CALL spline_fit(DBGeom%n_geom, DBGeom%dm3, DBGeom%volume, &
    DBGeom%d2_volume)
!
!
!
   RETURN
!
END SUBROUTINE ice_db_read_geometry
