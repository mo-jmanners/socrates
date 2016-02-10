! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to make spectral blocks of type 17.
!
! Description:
!   This routine creates a look-up table of spectral variability data
!
!------------------------------------------------------------------------------
SUBROUTINE make_block_17(Sp, Sol, ierr)

  USE realtype_rd, ONLY: RealK
  USE def_spectrum, ONLY: StrSpecData, StrSpecVar
  USE def_solarspec, ONLY: StrSolarSpec
  USE def_refract, ONLY: StrRefract

  IMPLICIT none

  TYPE(StrSpecData), Intent(INOUT) :: Sp
!   Spectral file data
  TYPE (StrSolarSpec), Intent(INOUT) :: Sol
!   Mean Solar spectrum
  INTEGER, INTENT(INOUT) :: ierr
!   Error flag

  TYPE(StrSpecData) :: SubSp
!   Temporary spectral file data for sub-bands
  TYPE(StrSpecVar) :: SpVarTmp
!   Temporary spectral variability data
  TYPE (StrSolarSpec) :: VSol
!   Varying Solar spectrum
  TYPE (StrRefract) :: Refract
!   Refractive index of atmosphere

  INTEGER :: ios
!   Reading error flag
  INTEGER :: i
!   Loop variable
  INTEGER :: iu_solar
!   Unit number for the solar spectrum data file
  LOGICAL :: l_count
!   Flag for counting points
  REAL (RealK) :: scale_wv, scale_irr
!   Scaling for wavelength and irradiance to correct units

  CHARACTER  (LEN=1) :: char
  CHARACTER  (LEN=80) :: line
  INTEGER :: band, sub_band, number_term
  INTEGER :: sub_bands(Sp%Dim%nd_band)
  INTEGER :: n_times
  LOGICAL :: l_monthly
  REAL (RealK) :: wavelength(2, Sp%Dim%nd_k_term, Sp%Dim%nd_band)
  REAL (RealK) :: wave_inc


  IF ( .NOT. Sp%Basic%l_present(17) ) THEN
    sub_bands=0
    sub_bands(1:Sp%Basic%n_band)=1
    DO
      WRITE(*, '(/a)') 'Enter band to be sub-divided (0 to finish): '
      READ(*, *, IOSTAT=ios) band
      IF (band == 0) EXIT
      number_term = Sp%Gas%i_band_k(band, Sp%Gas%index_absorb(1, band))
      sub_bands(band)=number_term
      WRITE(*, '(a,i5,a)') &
        'There are ', number_term, ' major gas k-terms in this band.'
      WRITE(*, '(a)') 'Do you want to divide equally in wavelength (E),'
      WRITE(*, '(a)') 'or provide band limits (L) for each k-term?'
      READ(*, '(a)') char
      IF ( (char.eq.'e').OR.(char.eq.'E') ) THEN
        wave_inc = ( Sp%Basic%wavelength_long(band) - &
                     Sp%Basic%wavelength_short(band) ) / number_term
        DO i=1, number_term
          wavelength(1, i, band) = Sp%Basic%wavelength_short(band) + &
                                   wave_inc*(i-1)
          wavelength(2, i, band) = Sp%Basic%wavelength_short(band) + &
                                   wave_inc*(i)
        END DO
      ELSE
        WRITE(*, '(a)') 'Enter band limits (metres): '
        DO i=1, number_term
          READ(*, *, IOSTAT=ios) wavelength(1, i, band), &
                                        wavelength(2, i, band)
        END DO 
      END IF
    END DO
    
    Sp%Var%n_sub_band = SUM(sub_bands)
    Sp%Dim%nd_sub_band = Sp%Var%n_sub_band
    
    IF (ALLOCATED(Sp%Var%index_sub_band)) &
       DEALLOCATE(Sp%Var%index_sub_band)
    IF (ALLOCATED(Sp%Var%wavelength_sub_band)) &
       DEALLOCATE(Sp%Var%wavelength_sub_band)
    ALLOCATE(Sp%Var%index_sub_band( 2, Sp%Dim%nd_sub_band ))
    ALLOCATE(Sp%Var%wavelength_sub_band( 2, Sp%Dim%nd_sub_band ))
    
    sub_band=1
    DO band=1, Sp%Basic%n_band
      IF (sub_bands(band) == 1) THEN
        Sp%Var%index_sub_band(1, sub_band) = band
        Sp%Var%index_sub_band(2, sub_band) = 0
        Sp%Var%wavelength_sub_band(1, sub_band) = &
          Sp%Basic%wavelength_short(band)
        Sp%Var%wavelength_sub_band(2, sub_band) = &
          Sp%Basic%wavelength_long(band)
        sub_band = sub_band + 1
      ELSE
        DO i=1, sub_bands(band)
          Sp%Var%index_sub_band(1, sub_band) = band
          Sp%Var%index_sub_band(2, sub_band) = i
          Sp%Var%wavelength_sub_band(1, sub_band) = wavelength(1, i, band)
          Sp%Var%wavelength_sub_band(2, sub_band) = wavelength(2, i, band)
          sub_band = sub_band + 1
        END DO
      END IF
    END DO

    Sp%Var%n_times  = 0
    Sp%Dim%nd_times = 0
    Sp%Var%n_repeat_times  = 0
    IF (ALLOCATED(Sp%Var%time)) &
       DEALLOCATE(Sp%Var%time)
    IF (ALLOCATED(Sp%Var%total_solar_flux)) &
       DEALLOCATE(Sp%Var%total_solar_flux)
    IF (ALLOCATED(Sp%Var%solar_flux_sub_band)) &
       DEALLOCATE(Sp%Var%solar_flux_sub_band)
    IF (ALLOCATED(Sp%Var%rayleigh_coeff)) &
       DEALLOCATE(Sp%Var%rayleigh_coeff)
    ALLOCATE(Sp%Var%time( 4, Sp%Dim%nd_times ))
    ALLOCATE(Sp%Var%total_solar_flux( Sp%Dim%nd_times ))
    ALLOCATE(Sp%Var%solar_flux_sub_band( Sp%Dim%nd_sub_band, Sp%Dim%nd_times ))
    ALLOCATE(Sp%Var%rayleigh_coeff( Sp%Dim%nd_sub_band, 0:Sp%Dim%nd_times ))
  END IF

! Fill temporary spectral type to hold sub-bands as full bands
  SubSp%Basic%n_band = Sp%Var%n_sub_band
  SubSp%Dim%nd_band  = Sp%Var%n_sub_band
  ALLOCATE(SubSp%Basic%wavelength_short(SubSp%Basic%n_band))
  SubSp%Basic%wavelength_short = Sp%Var%wavelength_sub_band(1, :)
  ALLOCATE(SubSp%Basic%wavelength_long(SubSp%Basic%n_band))
  SubSp%Basic%wavelength_long = Sp%Var%wavelength_sub_band(2, :)
  ALLOCATE(SubSp%Basic%l_present(0:14))
  SubSp%Basic%l_present(14) = .FALSE.
  ALLOCATE(SubSp%Solar%solar_flux_band(SubSp%Basic%n_band))
  ALLOCATE(SubSp%Rayleigh%rayleigh_coeff(SubSp%Basic%n_band))

  IF ( .NOT. Sp%Basic%l_present(17) ) THEN
    ! Fill sub-band Rayleigh coefficients for mean solar spectrum
    WRITE(*, '(/a)') 'A mean solar spectrum is needed for the mean'
    WRITE(*, '(a)')  'Rayleigh coefficients per sub-band:'
    CALL make_block_3(SubSp, Sol, ierr)
    Sp%Var%rayleigh_coeff(:,0) = SubSp%Rayleigh%rayleigh_coeff
  END IF

  DO
    IF (Sp%Var%n_times == 0) THEN
      ! If there are currently no times in the look-up table we can
      ! choose the number of varying Rayleigh coefficients to set 
      WRITE(*, '(a)') &
        'How many sub-bands will require a varying Rayleigh coefficient:'
      READ(*, *, IOSTAT=ios) Sp%Var%n_rayleigh_coeff
    END IF

    WRITE(*, '(a)') 'Number of times / dates to add to spectral data'
    WRITE(*, '(a)') '(0 to finish, -1 to clear current data):'
    READ(*, *, IOSTAT=ios) n_times
    IF (n_times == 0 .OR. ios /= 0) EXIT
    IF (n_times < 0) THEN
      Sp%Var%n_times  = 0
      Sp%Dim%nd_times = 0
      CYCLE
    END IF
    Sp%Dim%nd_times = Sp%Dim%nd_times + n_times

!   Save current variability data
    SpVarTmp = Sp%Var

!   Reallocate arrays to extend times
    DEALLOCATE(Sp%Var%time)
    DEALLOCATE(Sp%Var%total_solar_flux)
    DEALLOCATE(Sp%Var%solar_flux_sub_band)
    DEALLOCATE(Sp%Var%rayleigh_coeff)
    ALLOCATE(Sp%Var%time( 4, Sp%Dim%nd_times ))
    ALLOCATE(Sp%Var%total_solar_flux( Sp%Dim%nd_times ))
    ALLOCATE(Sp%Var%solar_flux_sub_band( Sp%Dim%nd_sub_band, Sp%Dim%nd_times ))
    ALLOCATE(Sp%Var%rayleigh_coeff( Sp%Dim%nd_sub_band, 0:Sp%Dim%nd_times ))
    Sp%Var%time(:, 1:Sp%Var%n_times) = &
      SpVarTmp%time(:, 1:Sp%Var%n_times)
    Sp%Var%total_solar_flux(1:Sp%Var%n_times) = &
      SpVarTmp%total_solar_flux(1:Sp%Var%n_times)
    Sp%Var%solar_flux_sub_band(:, 1:Sp%Var%n_times) = &
      SpVarTmp%solar_flux_sub_band(:, 1:Sp%Var%n_times)
    Sp%Var%rayleigh_coeff(:, 0:Sp%Var%n_times) = &
      SpVarTmp%rayleigh_coeff(:, 0:Sp%Var%n_times)

!   Read data from CMIP5 format spectral variability file
    CALL get_free_unit(ios, iu_solar)
    CALL open_file_in(ios, iu_solar, & 
      'Enter location of data file (see http://solarisheppa.geomar.de/cmip5):')

!   Read first to find the number of points in the spectrum.
    VSol%n_points = 0
    l_count=.FALSE.
    DO
      READ(iu_solar, '(A)', IOSTAT=ios) line
      IF (ios /= 0) THEN
        EXIT
      ELSE IF (INDEX(line, 'wavelength (nm) grid') > 1) THEN
        l_count=.TRUE.
      ELSE IF (INDEX(line, 'wavelength bands') > 1) THEN
        EXIT
      ELSE IF (l_count) THEN
        VSol%n_points=VSol%n_points+1
      ENDIF
    ENDDO
    VSol%n_points = VSol%n_points * 5 ! Five wavelengths per line

    ALLOCATE(VSol%wavelength(VSol%n_points))
    ALLOCATE(VSol%irrad(     VSol%n_points))

    scale_wv=1.0E-09_RealK
    scale_irr=0.9965E+06_RealK ! Wm-3 in TIM scale
    l_monthly = .FALSE.
    REWIND(iu_solar)
    DO
      READ(iu_solar, '(a)', IOSTAT=ios) line
      IF (ios /= 0) THEN
        print*, 'end of file'
        STOP
      END IF
      IF (INDEX(line, 'wavelength (nm) grid') > 1) THEN
        DO i = 1, VSol%n_points/5
          READ(iu_solar, *, IOSTAT=ios) VSol%wavelength(i*5-4:i*5)
        END DO
!       Read down to the start of the irradaince data
        DO
          READ(iu_solar, '(a)', IOSTAT=ios) line
          IF (INDEX(line, 'MONTH') > 1) THEN
            l_monthly = .TRUE.
          END IF
          IF (INDEX(line, 'TSI') > 1) EXIT
          IF (ios /= 0) STOP
        END DO
        EXIT
      ENDIF
    END DO
!   Scale the values to the correct units:
    VSol%wavelength = VSol%wavelength * scale_wv

!   Read in each time / date and calculate normalised spectrum
    DO i = Sp%Var%n_times + 1, Sp%Var%n_times + n_times
      IF (l_monthly) THEN
        READ(iu_solar, '(2i6,f17.6)', IOSTAT=ios) &
          Sp%Var%time(1:2, i), Sp%Var%total_solar_flux(i)
      ELSE
        READ(iu_solar, '(i6,f17.6)', IOSTAT=ios) &
          Sp%Var%time(1, i), Sp%Var%total_solar_flux(i)
        Sp%Var%time(2, i) = 1 ! Hardwire month to 1
      END IF

      ! Multiply by 0.9965 for TIM scale    
      Sp%Var%total_solar_flux(i) = Sp%Var%total_solar_flux(i) * 0.9965

      Sp%Var%time(3, i) = 1 ! Hardwire day of month to 1
      Sp%Var%time(4, i) = 0 ! Hardwire seconds since midnight to 0

      READ(iu_solar, *, IOSTAT=ios) VSol%irrad
      VSol%irrad = VSol%irrad * scale_irr

!     Calculate the normalised solar flux in each sub-band
      CALL make_block_2_1(SubSp, VSol, .TRUE., .FALSE., ierr)
      Sp%Var%solar_flux_sub_band(:,i) = SubSp%Solar%solar_flux_band

!     Calculate Rayleigh scattering coefficients in each sub-band
      CALL make_block_3_1(SubSp, VSol, Refract, .FALSE.)
      Sp%Var%rayleigh_coeff(:,i) = SubSp%Rayleigh%rayleigh_coeff
    END DO
    Sp%Var%n_times = Sp%Var%n_times + n_times

    DEALLOCATE(VSol%irrad)
    DEALLOCATE(VSol%wavelength)
    CLOSE(iu_solar)
  END DO

  WRITE(*, '(a)') 'How many of the final times / dates should be '
  WRITE(*, '(a)') 'periodically repeated into the future:'
  READ(*, *, IOSTAT=ios) Sp%Var%n_repeat_times

  Sp%Basic%l_present(17)=.TRUE.

END SUBROUTINE make_block_17
