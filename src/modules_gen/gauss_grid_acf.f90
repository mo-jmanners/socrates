! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to set parameters for calculating a Gaussian grid.

MODULE gauss_grid_acf

! Description:
!
! This module defines parameters for calculating a Gaussian grid
! for use in scattering calculations.

  IMPLICIT NONE

  INTEGER, Parameter  ::  NP_Gauss_maxit = 100
!   Maximum number of iterations to find Gaussian weights

END MODULE gauss_grid_acf
