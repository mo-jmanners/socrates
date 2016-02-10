! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to make spectral blocks of type 3.
!
! Method:
!	A solar spectrum is read if necessary. The monochromatic
!	Rayleigh scattering coefficients are calculated and 
!	weighted with the solar spectrum (in make_block_3_1).
!
!- ---------------------------------------------------------------------
SUBROUTINE make_block_3(Sp, Sol, ierr)

  USE realtype_rd
  USE def_spectrum
  USE rad_pcf
  USE dimensions_pp_ucf
  USE def_std_io_icf
  USE def_solarspec
  USE def_refract, ONLY: StrRefract, allocate_refract, deallocate_refract

  IMPLICIT NONE


  TYPE (StrSpecData), INTENT(INOUT) :: Sp
!   Spectral file to be assigned
  TYPE (StrSolarSpec), INTENT(INOUT) :: Sol
!   Solar spectrum
  INTEGER, Intent(INOUT) :: ierr
!   Error flag

! Local variables
  INTEGER :: i
!   Loop variable
  TYPE (StrRefract) :: Refract_H2
!   Refractive index of H2
  INTEGER :: iu_n_h2
!   File unit for H2 refractive index data
  INTEGER :: ios
!   I/O error flag
  CHARACTER(LEN=1) :: char_h2he_atm
!   Response to air/h2he question
  LOGICAL :: l_h2he_atm = .FALSE.
!   Calculate rayleigh scattering for H2/He atmosphere


  IF (ALLOCATED(Sp%Rayleigh%rayleigh_coeff)) &
      DEALLOCATE(Sp%Rayleigh%rayleigh_coeff)
  ALLOCATE(Sp%Rayleigh%rayleigh_coeff(Sp%Dim%nd_band))


! If a solar spectrum is already present that is used, otherwise
! one is read in.
  IF (Sol%n_points > 0) THEN
    WRITE(*, '(/a)') &
      'Rayleigh scattering coefficients will be averaged using '
    WRITE(*, '(a/)') &
      'the solar spectrum read in earlier.'
  ELSE
    CALL read_solar_spectrum(Sol, ierr)
    IF (ierr /= i_normal) RETURN
  ENDIF

! Is the atmosphere composed of air or H2-He gas?
  WRITE(*, '(/A)') 'Is the atmosphere composed of air or H2-He gas (A/H)?'
  DO
    READ(*, *, IOSTAT=ios) char_h2he_atm
    IF (ios /= 0) THEN
      WRITE(*, '(A)') '***error: unrecognized response'
      WRITE(*, '(A)') 'Please re-enter.'
    ELSE
      EXIT
    ENDIF
  ENDDO

! If gas is H2/He read in the refractive index for H2.
  IF (char_h2he_atm=='H' .OR. char_h2he_atm=='h') THEN

!   Obtain the file containing the H2 refractive index data.
    CALL get_free_unit(ierr, iu_n_h2)
    CALL open_file_in(ierr, iu_n_h2, & 
      'Enter the name of the file containing the H2 refractive index data.')
    IF (ierr /= i_normal) RETURN

!   Read first to find the number of points in the spectrum.
    Refract_H2%n_points = 0
    ios = 0
    DO WHILE (ios == 0)
      READ(iu_n_h2,'(f22.15)',IOSTAT=ios)
      Refract_H2%n_points = Refract_H2%n_points + 1
    END DO
    Refract_H2%n_points = Refract_H2%n_points - 2

    CALL allocate_refract(Refract_H2)

!   Read in the file.
    REWIND(iu_n_h2)
    READ(iu_n_h2, '(A)', IOSTAT=ios) ! Skip header line
    DO i = 1, Refract_H2%n_points
      READ(iu_n_h2, *, IOSTAT=ios) &
        Refract_H2%wavelength(i), &
        Refract_H2%re_part(i)
      Refract_H2%wavelength(i) = &
        Refract_H2%wavelength(i)*1e-6_RealK ! Convert to metre
      IF (ios /= 0) THEN
        WRITE(iu_err, '(/A)') '*** Error: Corrupt refractive index data.'
        ierr=i_err_fatal
        RETURN
      ENDIF
    ENDDO

!   Check that reading was carried out: failure to read the data
!   will cause n_points to be 0.
    IF (Refract_H2%n_points == 0) THEN
      WRITE(iu_err, '(/a)') &
        '*** Error: No data were read. ' // &
        'Check format of file of refractive index data.'
      ierr=i_err_fatal
      RETURN
    ENDIF

    l_h2he_atm = .TRUE.

  END IF

  CALL make_block_3_1(Sp, Sol, Refract_H2, l_h2he_atm)
  IF (l_h2he_atm) CALL deallocate_refract(Refract_H2)
  Sp%Basic%l_present(3)=.TRUE.

END SUBROUTINE make_block_3
