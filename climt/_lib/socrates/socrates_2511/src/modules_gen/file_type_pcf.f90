! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module defining types of intermediate files.

MODULE file_type_pcf

! Description:
!
! This defines the identification numbers used to define intermediate
! files generated during preprocessing.

  IMPLICIT NONE

  INTEGER, Parameter :: it_file_line_trans_form         =  2
!           Formatted line transmission data
  INTEGER, Parameter :: it_file_cont_trans_form         =  3
!           Formatted continuum transmission
  INTEGER, Parameter :: it_file_line_trans              = 12
!           File of line transmission data
  INTEGER, Parameter :: it_file_cont_trans              = 13
!           File of continuum transmission data
  INTEGER, Parameter :: it_file_line_trans_summary      = 14
!           File of line transmission data
  INTEGER, Parameter :: it_file_cont_trans_summary      = 15
!           File of continuum transmission data
  INTEGER, Parameter :: it_file_line_fit                =  4
!           File of fits to line data
  INTEGER, Parameter :: it_file_cont_fit                =  5
!           File of fits to continuum data
  INTEGER, Parameter :: it_file_mie_dry                 =  6
!           File of dry Mie data
  INTEGER, Parameter :: it_file_mie_humid               =  7
!           File of moist Mie data
  INTEGER, Parameter :: it_file_ave_mie_dry             =  8
!           File of dry averaged Mie data
  INTEGER, Parameter :: it_file_ave_mie_humid           =  9
!           File of moist averaged Mie data
  INTEGER, Parameter :: it_file_cloud_fit               = 10
!           File of fits to cloud data
  INTEGER, Parameter :: it_file_cloud_obs               = 11
!           File of observational cloud data
  INTEGER, Parameter :: it_file_adt_dry                 = 16
!           File of dry scattering properties from ADT
  INTEGER, Parameter :: it_file_adt_humid               = 20
!           File of dry scattering properties from ADT
  INTEGER, Parameter :: it_file_nsice_ave               = 17
!           File of averaged scattering properties for 
!           non-spherical ice
  INTEGER, Parameter :: it_file_phf_mie_dry             = 18
!           File of Mie scattering data with more than one moment
!           of the phase function
  INTEGER, Parameter :: it_file_phf_mie_humid           = 19
!           File of Mie scattering data with more than one moment
!           of the phase function
  INTEGER, Parameter :: it_file_ave_phf_mie_dry         = 21
!           File of averaged dry single scattering data including
!           higher moments of the phase function
  INTEGER, Parameter :: it_file_ave_phf_mie_humid       = 22
!           File of averaged moist single scattering data including
!           higher moments of the phase function
  INTEGER, Parameter :: it_file_scat_database           = 23
!           File of scattering properties from scattering
!           database with more than one moment
!           of the phase function
  INTEGER, Parameter :: it_file_cloud_fit_phf           = 24
!           File of fitted single scattering properties
!           including higher momnets of the phase function
  INTEGER, Parameter :: it_file_cont_gen_fit            = 25
!           File of fits to generalised continuum data
  INTEGER, Parameter :: it_file_line_fit_self           = 26
!           File of fits to line data that includes self-broadening
  INTEGER, Parameter :: it_file_line_fit_id             = 27
!           File of fits to line data using the gas identifier rather than index
  INTEGER, Parameter :: it_file_line_fit_self_id        = 28
!           File of fits to line data including self-broadening using the gas id
  INTEGER, Parameter :: it_file_scat_mass               = 29
!           File of single scattering properties as a function of
!           mass mixing ratio and mean particle mass
  INTEGER, Parameter :: it_file_ave_scat_mass           = 30
!           File of averaged single scattering properties as a
!           function of mass mixing ratio and mean particle mass

END MODULE file_type_pcf
