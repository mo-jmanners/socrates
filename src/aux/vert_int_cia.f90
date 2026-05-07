! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Program to integrate a CIA column vertically through the atmosphere.
!
! Method:
!   Two CDL-files are read in. The contents are assumed to be mixing
!   ratios on pressure levels. These pressure levels give the
!   masses in the layers so defined. The mean mixing ratio in the
!   layer is taken as the average of the values at the boundaries.
!   This enables us to calculate the amount of the quantity
!   integrating downwards. The results are written to a CDL-file.
!
!-----------------------------------------------------------------------
PROGRAM vert_int_cia

! Modules to set types of variables:
  USE realtype_rd, ONLY: RealK
  USE dimensions_field_ucf, ONLY: npd_profile, npd_latitude, &
    npd_longitude, npd_layer
  USE dimensions_cdl_ucf, ONLY: npd_cdl_dimen, npd_cdl_dimen_size, &
    npd_cdl_data, npd_cdl_var
  USE def_std_io_icf, ONLY: iu_stdin, iu_stdout
  USE error_pcf, ONLY: i_normal
  USE rad_ccf, ONLY: set_socrates_constants, &
    grav_acc, r_gas_dry

  IMPLICIT NONE

! Declaration of variables.
  CHARACTER (LEN=256) :: file_mix1, file_mix2
!       Name of files with mixing ratios
  CHARACTER (LEN=256) :: file_t
!       Name of file with temperatures
  CHARACTER (LEN=256) :: nml_file = ""
!       Name of namelist file
  CHARACTER (LEN=256) :: file_out
!       Name of file for output
  CHARACTER (LEN=24) :: name_vert_coord
!       Name of vertical coordinate
  LOGICAL :: l_vert_coord
!       Flag asserting that vertical coordinate is set
  INTEGER :: ierr, ios
!       Error flag
  INTEGER :: n_latitude
!       Number of latitudes
  INTEGER :: n_longitude
!       Number of longitudes
  INTEGER :: n_profile
!       Number of profiles
  INTEGER :: n_level
!       Number of levels
  INTEGER :: iu_nml
!       Namelist file unit
  INTEGER :: i, l
!       Loop variables
  REAL (RealK) :: latitude(npd_latitude)
!       Latitude
  REAL (RealK) :: longitude(npd_longitude)
!       Longitude
  REAL (RealK) :: p(npd_profile, npd_layer+1)
!       Pressure levels
  REAL (RealK) :: mix_ratio1(npd_profile, npd_layer+1)
!       Mixing ratio of gas 1
  REAL (RealK) :: mix_ratio2(npd_profile, npd_layer+1)
!       Mixing ratio of gas 2
  REAL (RealK) :: t(npd_profile, npd_layer+1)
!       Temperature
  REAL (RealK) :: amount(npd_profile, npd_layer+1)
!       Amounts of quantity
  REAL (RealK) :: mass, density
!       Working mass and density of layer

! Subroutines called:
  EXTERNAL :: assign_input_vert_cdl, output_vert_cdl


  data l_vert_coord/.false./
  data ierr/i_normal/
  data n_latitude/0/
  data n_longitude/0/

  WRITE(iu_stdout, '(/a)') &
    'Filename of gas 1 mixing ratios:'
  READ(iu_stdin, '(a)') file_mix1
  CALL assign_input_vert_cdl(ierr,                             &
    file_mix1, 'mixing ratios', l_vert_coord, name_vert_coord, &
    .true., n_level, .NOT.l_vert_coord,                        &
    n_latitude, latitude, n_longitude, longitude,              &
    1,                                                         &
    n_profile, n_level,                                        &
    p, mix_ratio1,                                             &
    npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1,  &
    npd_cdl_dimen, npd_cdl_dimen_size,                         &
    npd_cdl_data, npd_cdl_var)
  IF (ierr /= i_normal) STOP

  WRITE(iu_stdout, '(/a)') &
    'Filename of gas 2 mixing ratios:'
  READ(iu_stdin, '(a)') file_mix2
  CALL assign_input_vert_cdl(ierr,                             &
    file_mix2, 'mixing ratios', l_vert_coord, name_vert_coord, &
    .true., n_level, .NOT.l_vert_coord,                        &
    n_latitude, latitude, n_longitude, longitude,              &
    1,                                                         &
    n_profile, n_level,                                        &
    p, mix_ratio2,                                             &
    npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1,  &
    npd_cdl_dimen, npd_cdl_dimen_size,                         &
    npd_cdl_data, npd_cdl_var)
  IF (ierr /= i_normal) STOP

  WRITE(iu_stdout, '(/a)') &
    'Filename of temperatures:'
  READ(iu_stdin, '(a)') file_t
  CALL assign_input_vert_cdl(ierr,                             &
    file_t, 'mixing ratios', l_vert_coord, name_vert_coord,    &
    .true., n_level, .NOT.l_vert_coord,                        &
    n_latitude, latitude, n_longitude, longitude,              &
    1,                                                         &
    n_profile, n_level,                                        &
    p, t,                                                      &
    npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1,  &
    npd_cdl_dimen, npd_cdl_dimen_size,                         &
    npd_cdl_data, npd_cdl_var)
  IF (ierr /= i_normal) STOP

  WRITE(iu_stdout, "(a)") "Enter the name of the namelist file."
  READ(iu_stdin, "(a)") nml_file
  IF (nml_file /= " ") THEN
    OPEN(NEWUNIT=iu_nml, FILE=nml_file, IOSTAT=ios, &
     STATUS='OLD', ACTION='READ')
    IF (ios /= 0) THEN
      WRITE(iu_stdout, "(a)") 'Namelist file could not be opened.'
      STOP
    END IF
    CALL set_socrates_constants(iu_nml)
    CLOSE(iu_nml)
  END IF

  DO l=1, n_profile
    amount(l, 1)=0.0_RealK
  END DO
  DO i=2, n_level
    DO l=1, n_profile
      mass = (p(l, i)-p(l, i-1))/grav_acc
      density = ( p(l, i)/(r_gas_dry*t(l, i)) + &
        p(l, i-1)/(r_gas_dry*t(l, i-1)) ) / 2.0_RealK
      amount(l, i)=amount(l, i-1)                      &
        + (0.5_RealK*(mix_ratio1(l, i-1)+mix_ratio1(l, i))) &
        * (0.5_RealK*(mix_ratio2(l, i-1)+mix_ratio2(l, i))) &
        * mass * density
    END DO
  END DO

  WRITE(iu_stdout, '(/a)') 'enter the name of the output file.'
  READ(iu_stdin, '(a)') file_out
  CALL output_vert_cdl(ierr,                                       &
    file_out,                                                      &
    n_latitude, latitude, n_longitude, longitude, n_profile,       &
    n_level, trim(name_vert_coord), len(trim(name_vert_coord)), p, &
    'vint', 'float', 'kg2 m-5',                                     &
    'cia column', amount,                                          &
    npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1,      &
    npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var)
  IF (ierr /= i_normal) STOP

END PROGRAM vert_int_cia
