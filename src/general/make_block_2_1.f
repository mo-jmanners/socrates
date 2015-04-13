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
      SUBROUTINE make_block_2_1(ierr
     &  , n_band, wave_length_short, wave_length_long
     &  , l_exclude, n_band_exclude, index_exclude
     &  , l_solar_spectrum, n_solar_points
     &  , solar_wavelength, solar_irrad
     &  , solar_flux_band, l_present_2
     &  )

      USE realtype_rd
      USE rad_pcf
      USE dimensions_spec_ucf
      USE dimensions_pp_ucf
      USE def_std_io_icf

      IMPLICIT NONE


!     Dummy arguments
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
      INTEGER, Intent(IN) ::
     &    n_band
!           Number of bands
      INTEGER, Intent(INOUT) ::
     &    n_solar_points
!           Number of points in spectrum
      LOGICAL, Intent(INOUT) ::
     &    l_solar_spectrum
!           Solar spectral flag
      LOGICAL, Intent(IN) ::
     &    l_exclude
!           Flag for exclusion of bands
      INTEGER, Intent(IN) ::
     &    n_band_exclude(n_band)
!           Number of bands excluded from specific bands
     &  , index_exclude(npd_exclude, n_band)
!           Indices of excluded bands
      LOGICAL, Intent(OUT) ::
     &    l_present_2
!           Presence flag
      REAL  (RealK), Intent(OUT) ::
     &    solar_flux_band(n_band)
!           Solar flux in each band
      REAL  (RealK), Intent(IN) ::
     &    wave_length_short(n_band)
!           Bottom of each band
     &  , wave_length_long(n_band)
!           Top of each band
      REAL  (RealK), Intent(INOUT) ::
     &    solar_wavelength(npd_solar_points)
!           Wavelengths of solar spectrum
     &  , solar_irrad(npd_solar_points)
!           Solar irradiance at toa
!
!     Local variables
!
      CHARACTER
     &    char_yn
!           Character response variable
      INTEGER
     &    i_short
!           Beginning of band
     &  , i_long
!           End of band
     &  , n_points_band
!           Number of points within band
     &  , i
!           Loop variable
     &  , j
!           Loop variable
     &  , i_last_band
!           Index of last band
      LOGICAL
     &    l_enhance
!           Enhance outer bands
      REAL  (RealK) ::
     &    total_solar_irradiance
!           Total irradiance
     &  , x(npd_solar_points)
!           Absicissae of wavelength
     &  , y(npd_solar_points)
!           Ordinates of spectrum
     &  , irradiance_tail
!           Irradiance in tail
     &  , wave_length_begin
!           Beginning of range
     &  , wave_length_end
!           End of range
!
!     Subroutines called:
      EXTERNAL
     &    read_solar_spectrum, inner_bracket
!     Functions called:
      REAL  (RealK) ::
     &    solar_intensity
!           Solar intensity at given w.leng.
     &  , trapezoid
!           Integrating function
     &  , rayleigh_jeans_tail
!           Tail of rayleigh jeans function
      EXTERNAL
     &    solar_intensity, trapezoid, rayleigh_jeans_tail
!
!
!
!     Obtain the solar spectrum if data are not already present.
      IF (l_solar_spectrum) THEN
        WRITE(iu_stdout, '(/a/)')
     &    'Previous solar spectrum will be used.'
      ELSE
        CALL read_solar_spectrum(ierr, l_solar_spectrum
     &    , n_solar_points, solar_wavelength, solar_irrad)
        IF (ierr /= i_normal) RETURN
      ENDIF
!
!     The radiance outside the nominal limits of the spectrum can be
!     assigned to the edging bands.
      WRITE(iu_stdout, '(/a)') 'Assign solar flux outside given bands '
     &  //'to outside bands? (y/n)'
1     READ(iu_stdin, '(a)') char_yn
      IF ( (char_yn == 'Y').OR.(char_yn == 'y') ) THEN
        l_enhance=.true.
      ELSE IF ( (char_yn == 'N').OR.(char_yn == 'n') ) THEN
        l_enhance=.false.
      ELSE
        WRITE(iu_err, '(a)') '+++ Unrecognised response: '
        WRITE(iu_stdout, '(a)') 'Please re-type.'
        goto 1
      ENDIF
!
!     Calculate the total integrated irradiance of this spectrum.
      total_solar_irradiance=trapezoid(n_solar_points
     &  , solar_wavelength, solar_irrad)
!     Add on the tail of the distribution using a rayleigh-jeans law.
      irradiance_tail=
     &  rayleigh_jeans_tail(solar_wavelength(n_solar_points))
      total_solar_irradiance=total_solar_irradiance+irradiance_tail
!
!     For each band find the arrays of wavelength and irradiance.
      DO i=1, n_band
!
!       Find the points of the spectrum just within the band.
        IF (((i /= 1).AND.(i /= n_band)).OR.(.NOT.l_enhance)) THEN
          wave_length_begin=wave_length_short(i)
          wave_length_end=wave_length_long(i)
        ELSE IF (i == 1) THEN
          wave_length_begin=min(solar_wavelength(1)
     &      , wave_length_short(1))
          wave_length_end=wave_length_long(1)
        ELSE IF (i == n_band) THEN
          wave_length_begin=wave_length_short(n_band)
          wave_length_end=max(solar_wavelength(n_solar_points)
     &       , wave_length_long(n_band))
        ENDIF
!
        CALL inner_bracket(ierr
     &    , wave_length_begin, wave_length_end
     &    , n_solar_points, solar_wavelength, i_short, i_long)
        IF (ierr /= i_normal) RETURN
        n_points_band=3+i_long-i_short
        x(1)=wave_length_begin
        y(1)=solar_intensity(x(1)
     &    , n_solar_points, solar_wavelength, solar_irrad)
        x(n_points_band)=wave_length_end
        y(n_points_band)=solar_intensity(x(n_points_band)
     &    , n_solar_points, solar_wavelength, solar_irrad)
        DO j=i_short, i_long
          x(j-i_short+2)=solar_wavelength(j)
          y(j-i_short+2)=solar_irrad(j)
        ENDDO
!       Integrate across the band and normalize by 
!       the total irradiance.
        solar_flux_band(i)=trapezoid(n_points_band, x, y)
     &    /total_solar_irradiance
      ENDDO
!
!     Add the tail to the last band if required. the setting of limits
!     above dealt with the contributions within the region of solar
!     data.
      IF (l_enhance) THEN
        i_last_band=1
        DO i=1, n_band
          IF (wave_length_long(i) > wave_length_long(i_last_band))
     &      i_last_band=i
        ENDDO
        solar_flux_band(i_last_band)=solar_flux_band(i_last_band)
     &    +irradiance_tail/total_solar_irradiance
      ENDIF
!
!     Reduce the solar flux in bands which contain exclusions if
!     necessary.
      IF (l_exclude) THEN
        DO i=1, n_band
          DO j=1, n_band_exclude(i)
            solar_flux_band(i)=solar_flux_band(i)
     &        -solar_flux_band(index_exclude(j, i))
          ENDDO
        ENDDO
      ENDIF
!
      l_present_2=.true.
!
!
!
      RETURN
      END
