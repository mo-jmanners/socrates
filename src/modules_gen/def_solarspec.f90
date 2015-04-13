! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to declare a structure for solar spectral data.

MODULE def_solarspec

  USE realtype_rd

  IMPLICIT NONE


  TYPE StrSolarSpec

    INTEGER :: n_points
!     Number of points in the spectrum
    REAL  (RealK), Pointer :: wavelength(:)
!     Wavelengthe at which the spectral irradiance is specified
    REAL  (RealK), Pointer :: irrad(:)
!     Solar spectral irradiance in units of Wm-2.m-1

  END TYPE StrSolarSpec

END MODULE def_solarspec
