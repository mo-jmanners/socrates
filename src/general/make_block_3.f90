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
!	Rayeligh sattering coefficients are calculated and 
!	weighted with the solar spectrum.
!
!- ---------------------------------------------------------------------
SUBROUTINE make_block_3(Spectrum, l_solar_spectrum, n_solar_points,  &
  solar_wavelength, solar_irrad, ierr)

  USE realtype_rd
  USE def_spectrum
  USE rad_pcf
  USE dimensions_pp_ucf
  USE def_std_io_icf

  IMPLICIT NONE


  TYPE (StrSpecData), Intent(INOUT), TARGET :: Spectrum
!   Spectral file to be assigned
  LOGICAL, Intent(INOUT) :: l_solar_spectrum
!   Solar spectral flag
  INTEGER, Intent(INOUT) :: n_solar_points
!   Number of points in spectrum
  REAL (RealK), Intent(INOUT) :: solar_wavelength(npd_solar_points)
!   Wavelengths of solar spectrum
  REAL (RealK), Intent(INOUT) :: solar_irrad(npd_solar_points)
!   Solar irradiance at toa
  INTEGER, Intent(INOUT) :: ierr
!   Error flag

! Local variables
  INTEGER :: n_int
!   Number of intervals for integration
  INTEGER :: i_begin
!   Beginning of spectrum in band
  INTEGER :: i_end
!   End of spectrum in band
  INTEGER :: i_dummy
!   Dummy reading variable
  INTEGER :: i, k
!   Loop variables
  REAL (RealK) :: wave_int(npd_solar_points)
!   Wavelngths for integration
  REAL (RealK) :: weight_int(npd_solar_points)
!   Weighting solar irradiances
  REAL (RealK) :: product_int(npd_solar_points)
!   Irradiance times rayleigh coef.
  CHARACTER(LEN=1) :: l_h2he_atm_char
!   Response to air/h2he question
  INTEGER :: iu_n_h2
!   File unit for H2 refractive index data
  INTEGER :: n_points
!   Number of points in refractive index data
  INTEGER :: ios
!   I/O error flag
  REAL (RealK), DIMENSION(:), ALLOCATABLE :: wavelength_refract_index_H2
!   Wavelength at which H2 refractive index data is evaluated
  REAL (RealK), DIMENSION(:), ALLOCATABLE :: refract_index_H2
!   Refractive index of H2

! Functions called:
  REAL (RealK), EXTERNAL :: rayleigh_scatter_air
  REAL (RealK), EXTERNAL :: rayleigh_scatter_h2he
!   Functions for rayleigh scattering
  REAL (RealK), EXTERNAL :: solar_intensity
!   Function for intensity
  REAL (RealK), EXTERNAL :: trapezoid
!   Trapezoidal integration function


  IF (ALLOCATED(Spectrum%Rayleigh%rayleigh_coeff)) &
      DEALLOCATE(Spectrum%Rayleigh%rayleigh_coeff)
  ALLOCATE(Spectrum%Rayleigh%rayleigh_coeff(Spectrum%Dim%nd_band))


! If a solar spectrum is already present that is used, otherwise
! one is read in.
  IF (l_solar_spectrum) THEN
    WRITE(*, '(/a)') &
      'rayleigh scattering coefficients will be averaged using '
    WRITE(*, '(a/)') &
      'the solar spectrum read in earlier.'
  ELSE
    CALL read_solar_spectrum(ierr, l_solar_spectrum, &
      n_solar_points, solar_wavelength, solar_irrad)
    IF (ierr /= i_normal) RETURN
  ENDIF

! Is the atmosphere composed of air or H2-He gas?
  WRITE(*, '(/A)') 'Is the atmosphere composed of air or H2-He gas (A/H)?'
  DO
    READ(*, *, IOSTAT=ios) l_h2he_atm_char
    IF (ios /= 0) THEN
      WRITE(*, '(A)') '***error: unrecognized response'
      WRITE(*, '(A)') 'Please re-enter.'
    ELSE
      EXIT
    ENDIF
  ENDDO

! In each band the scattering coefficient is averaged at the points
! of the solar spectrum and the limits of the band.
  IF (l_h2he_atm_char=='H' .OR. l_h2he_atm_char=='h') THEN

!   Obtain the file containing the H2 refractive index data.
    CALL get_free_unit(ierr, iu_n_h2)
    CALL open_file_in(ierr, iu_n_h2, & 
      'Enter the name of the file containing the H2 refractive index data.')
    IF (ierr /= i_normal) RETURN

!   Read first to find the number of points in the spectrum.
    n_points = 0
    ios = 0
    DO WHILE (ios == 0)
      READ(iu_n_h2,'(f22.15)',IOSTAT=ios)
      n_points = n_points + 1
    END DO
    n_points = n_points - 2

    ALLOCATE(wavelength_refract_index_H2(n_points))
    ALLOCATE(refract_index_H2(n_points))

!   Read in the file.
    REWIND(iu_n_h2)
    READ(iu_n_h2, '(A)', IOSTAT=ios) ! Skip header line
    DO i = 1, n_points
      READ(iu_n_h2, *, IOSTAT=ios) &
        wavelength_refract_index_H2(i), &
        refract_index_H2(i)
      wavelength_refract_index_H2(i) = &
        wavelength_refract_index_H2(i)*1e-6_RealK ! Convert to metre
      IF (ios /= 0) THEN
        WRITE(iu_err, '(/A)') '*** Error: Corrupt refractive index data.'
        ierr=i_err_fatal
        RETURN
      ENDIF
    ENDDO

!   Check that reading was carried out: failure to read the data
!   will cause n_solar_points to be 0.
    IF (n_points == 0) THEN
      WRITE(iu_err, '(/a)') &
        '*** Error: No data were read. ' // &
        'Check format of file of refractive index data.'
      ierr=i_err_fatal
      RETURN
    ENDIF
  
    DO i=1, Spectrum%Basic%n_band
      CALL point_bracket(Spectrum%Basic%wavelength_short(i), &
        n_solar_points, solar_wavelength, i_dummy, i_begin)
      CALL point_bracket(Spectrum%Basic%wavelength_long(i), &
        n_solar_points, solar_wavelength, i_end, i_dummy)
!     I_BEGIN and I_END are the points of the solar spectrum
!     just within the band. form an array of wavelength points for
!     integration
      n_int=1
      wave_int(n_int)=Spectrum%Basic%wavelength_short(i)
      weight_int(n_int)=solar_intensity(wave_int(n_int), &
        n_solar_points, solar_wavelength, solar_irrad)
      product_int(n_int)=weight_int(n_int) &
        *rayleigh_scatter_h2he(wave_int(n_int), &
        wavelength_refract_index_H2, refract_index_H2, n_points)
      DO k=i_begin, i_end
        n_int=n_int+1
        wave_int(n_int)=solar_wavelength(k)
        weight_int(n_int)=solar_intensity(wave_int(n_int), &
          n_solar_points, solar_wavelength, solar_irrad)
        product_int(n_int)=weight_int(n_int) &
          *rayleigh_scatter_h2he(wave_int(n_int), &
          wavelength_refract_index_H2, refract_index_H2, n_points)
      ENDDO
      n_int=n_int+1
      wave_int(n_int)=Spectrum%Basic%wavelength_long(i)
      weight_int(n_int)=solar_intensity(wave_int(n_int), &
        n_solar_points, solar_wavelength, solar_irrad)
      product_int(n_int)=weight_int(n_int) &
        *rayleigh_scatter_h2he(wave_int(n_int), &
        wavelength_refract_index_H2, refract_index_H2, n_points)
      Spectrum%Rayleigh%rayleigh_coeff(i)= &
        trapezoid(n_int, wave_int, product_int)
      Spectrum%Rayleigh%rayleigh_coeff(i)= &
        Spectrum%Rayleigh%rayleigh_coeff(i) &
        /trapezoid(n_int, wave_int, weight_int)
    ENDDO
  ELSE
    DO i=1, Spectrum%Basic%n_band
      CALL point_bracket(Spectrum%Basic%wavelength_short(i), &
        n_solar_points, solar_wavelength, i_dummy, i_begin)
      CALL point_bracket(Spectrum%Basic%wavelength_long(i), &
        n_solar_points, solar_wavelength, i_end, i_dummy)
!     I_BEGIN and I_END are the points of the solar spectrum
!     just within the band. form an array of wavelength points for
!     integration
      n_int=1
      wave_int(n_int)=Spectrum%Basic%wavelength_short(i)
      weight_int(n_int)=solar_intensity(wave_int(n_int), &
        n_solar_points, solar_wavelength, solar_irrad)
      product_int(n_int)=weight_int(n_int) &
        *rayleigh_scatter_air(wave_int(n_int))
      DO k=i_begin, i_end
        n_int=n_int+1
        wave_int(n_int)=solar_wavelength(k)
        weight_int(n_int)=solar_intensity(wave_int(n_int), &
          n_solar_points, solar_wavelength, solar_irrad)
        product_int(n_int)=weight_int(n_int) &
          *rayleigh_scatter_air(wave_int(n_int))
      ENDDO
      n_int=n_int+1
      wave_int(n_int)=Spectrum%Basic%wavelength_long(i)
      weight_int(n_int)=solar_intensity(wave_int(n_int), &
        n_solar_points, solar_wavelength, solar_irrad)
      product_int(n_int)=weight_int(n_int) &
        *rayleigh_scatter_air(wave_int(n_int))
      Spectrum%Rayleigh%rayleigh_coeff(i)= &
        trapezoid(n_int, wave_int, product_int)
      Spectrum%Rayleigh%rayleigh_coeff(i)= &
        Spectrum%Rayleigh%rayleigh_coeff(i) &
        /trapezoid(n_int, wave_int, weight_int)
    ENDDO
  END IF

  Spectrum%Basic%l_present(3)=.TRUE.

END SUBROUTINE make_block_3
