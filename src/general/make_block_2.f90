! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to make spectral blocks of type 2.
!
! Method:
!	A solar spectrum, giving the irradiance against wavelength,
!	is read in. The total flux in the spectrum is determined,
!	using a Rayleigh-Jeans distribution for the far infra-red.
!	The fraction of this flux in each band is calculated. The
!	band will probably not cover the whole spectrum. If 
!	simulating observations this is as required, since the
!	spectral cut-offs of the observing instruments must be 
!	matched. In atmospheric models the whole flux should be
!	calculated and there is therefore an option to enhance
!	the fluxes in the outside bands to include the full flux.
!
!- ---------------------------------------------------------------------
SUBROUTINE make_block_2(Spectrum, l_solar_spectrum, n_solar_points,  &
  solar_wavelength, solar_irrad, ierr)

  USE realtype_rd
  USE def_spectrum
  USE rad_pcf
  USE dimensions_pp_ucf

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

! Local Variables
  INTEGER :: ios
!   IO status
  CHARACTER (LEN=1) :: l_filter
!   Character flag for filter function


  IF (ALLOCATED(Spectrum%Solar%solar_flux_band)) &
      DEALLOCATE(Spectrum%Solar%solar_flux_band)
  ALLOCATE(Spectrum%Solar%solar_flux_band(Spectrum%Dim%nd_band))

  WRITE(*, '(/A)') 'Is a filter function required (Y/N)?'
  DO
    READ(*, *, IOSTAT=ios) l_filter
    IF (ios /= 0) THEN
      WRITE(*, '(A)') '***error: unrecognized response'
      WRITE(*, '(A)') 'Please re-enter.'
    ELSE
      EXIT
    END IF
  END DO

  IF ((l_filter == 'Y').OR.(l_filter == 'y')) THEN
    CALL make_block_2_2(ierr, &
      Spectrum%Basic%n_band, &
      Spectrum%Basic%wavelength_short, &
      Spectrum%Basic%wavelength_long, &
      Spectrum%Basic%l_present(14), &
      Spectrum%Basic%n_band_exclude, &
      Spectrum%Basic%index_exclude, &
      l_solar_spectrum, n_solar_points, &
      solar_wavelength, solar_irrad, &
      Spectrum%Solar%solar_flux_band, &
      Spectrum%Basic%l_present(2) )
  ELSE
    CALL make_block_2_1(ierr, &
      Spectrum%Basic%n_band, &
      Spectrum%Basic%wavelength_short, &
      Spectrum%Basic%wavelength_long, &
      Spectrum%Basic%l_present(14), &
      Spectrum%Basic%n_band_exclude, &
      Spectrum%Basic%index_exclude, &
      l_solar_spectrum, n_solar_points, &
      solar_wavelength, solar_irrad, &
      Spectrum%Solar%solar_flux_band, &
      Spectrum%Basic%l_present(2) )
  END IF

END SUBROUTINE make_block_2
