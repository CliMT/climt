! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to methods of merging vertical profiles.
!
MODULE method_merge_pcf
!
! Description:
!
!   Vertical profiles for the radiation code must be specified over
!   the whole atmospheric column. Input from different sources may
!   be available onluy over restricted ranges and must therefore
!   be merged. This modules defines methods for doing that.
!
!- End of header
!
!
!
  IMPLICIT NONE
!
!
  INTEGER, Parameter :: P_method_merge   = 3
!   Size allocated for methods of merging data
!
  INTEGER, Parameter :: IP_merge_direct  = 1
!   Overlay data directly
  INTEGER, Parameter :: IP_merge_linear  = 2
!   Merge data in linearly
  INTEGER, Parameter :: IP_merge_zero    = 3
!   Overlap and zero background
!
!
  CHARACTER  (LEN=50), Dimension(P_method_merge) :: &
                                     name_method_merge=(/ &
    'Overlay data directly above background.           ', &
    'Merge data linearly over background.              ', &
    'Merge data directly with zero background.         ' /)
!
!
!
END MODULE method_merge_pcf
