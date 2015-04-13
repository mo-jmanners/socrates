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
!       This routine does the same thing as make_block_2_1 except
!       that it takes into account an experimental filter function.
!       This is mainly used for calculating radiances.
!
!- --------------------------------------------------------------------
      SUBROUTINE make_block_2_2(ierr
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
      USE def_inst_flt

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
     &    n_band_exclude(npd_band)
!           Number of bands excluded from specific bands
     &  , index_exclude(npd_exclude, npd_band)
!           Indices of excluded bands
      LOGICAL, Intent(OUT) ::
     &    l_present_2
!           Presence flag
      REAL  (RealK), Intent(OUT) ::
     &    solar_flux_band(npd_band)
!           Solar flux in each band
      REAL  (RealK), Intent(IN) ::
     &    wave_length_short(npd_band)
!           Bottom of each band
     &  , wave_length_long(npd_band)
!           Top of each band
      REAL  (RealK), Intent(INOUT) ::
     &    solar_wavelength(npd_solar_points)
!           Wavelengths of solar spectrum
     &  , solar_irrad(npd_solar_points)
!           Solar irradiance at toa
!
!     Local variables
!
      LOGICAL
     &    L_LAMBDA
!           Flag is set true if wavelength limits is filter 
!           and spectral agree
     &  , L_EXIST
!           Check whether file exists
!
      CHARACTER
     &    char_yn
!           Character response variable
      INTEGER
     &    I_SHORT
!           Beginning of band
     &  , I_LONG
!           End of band
     &  , I
!           Loop variable
     &  , J
!           Loop variable
     &  , i_last_band
!           Index of last band
     &  , N_SOLAR_INTERPOLATED
!            Number of points the solar spectrum is interpolated onto.     
      LOGICAL
     &    l_enhance
!           Enhance outer bands
      REAL  (RealK) ::
     &     LAMBDA_MIN
!           Smallest Wavelength
     &  ,  LAMBDA_MAX 
!           Largest Wavelength
     &  ,  D_LAMBDA 
!           Wavelength Increment
     &  ,  EPS1,EPS2
!           Small numbers
     &  ,  SOLAR_LAMBDA_MAX
!            Largest wavelength in solar spectrum     
     &  ,  SOLAR_LAMBDA_MIN
!            Smallest wavelength in solar spectrum     
     &  ,  BAND_IRRADIANCE
!            irradiance in band     
     &  ,  TOTAL_SOLAR_IRRADIANCE
!           Total irradiance
     &  , irradiance_tail
!           Irradiance in tail
     &  , C1
!           Dummy Variable     
      REAL  (RealK) ::
     &    solar_intensity
!           Solar intensity at given w.leng.
     &  , trapezoid
!           Integrating function
     &  , rayleigh_jeans_tail
!           Tail of rayleigh jeans function
     
       TYPE  (StrFiltResp) :: filter
!   Instrumental response function

      REAL (KIND=KIND(1.0D0)), ALLOCATABLE ::
     &    WEIGHTS(:)
!           Weights contained in the filter function
     &  , LAMBDA(:)
!           Wavelengths corresponding to these waves  
     &  , IRRAD_INTERPOL(:)
!            Interpolated solar spectrum  
     &  , LAMBDA_INTERPOL(:)
!            and the corresponding wavelength     

!     Subroutines called:
      EXTERNAL
     &    read_solar_spectrum, inner_bracket
     &  , solar_intensity, trapezoid, rayleigh_jeans_tail
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
2     READ(iu_stdin, '(a)') char_yn
      IF ( (char_yn == 'Y').OR.(char_yn == 'y') ) THEN
        l_enhance=.true.
      ELSE IF ( (char_yn == 'N').OR.(char_yn == 'n') ) THEN
        l_enhance=.false.
      ELSE
        WRITE(iu_err, '(a)') '+++ Unrecognised response: '
        WRITE(iu_stdout, '(a)') 'Please re-type.'
        goto 2
      ENDIF
      
      DO I=1,N_BAND
     
! Read in filter function from file

         CALL READ_INSTRUMENT_RESPONSE_90(FILTER,IERR)

! and transform wavenumber back into wavelength. Here we need to
! be slightly careful since READ_INSTRUMENT_RESPONSE_90 reads in
! the wavelength transforms them into wavenumbers (m-1) and then
! reorders them in order of increasing wavenumber! 

! Since all the files are in wavelength (and I prefer working in
! wavelength) we change back to wavelength.

         ALLOCATE(LAMBDA(FILTER%N_PTS))
         ALLOCATE(WEIGHTS(FILTER%N_PTS))
         
         DO J=1,FILTER%N_PTS
            LAMBDA(FILTER%N_PTS-J+1)=1.0/FILTER%WAVENUMBER(J)
         ENDDO
!
! Compare wavelength speciefied in the spectral file to those specified
! in the filter file.

         LAMBDA_MIN=WAVE_LENGTH_SHORT(I) 
         LAMBDA_MAX=WAVE_LENGTH_LONG(I)
      
         D_LAMBDA=(LAMBDA_MAX-LAMBDA_MIN)/FLOAT(FILTER%N_PTS-1)
!         print*,'D_LAMBDA=',D_LAMBDA    
          
           
         L_LAMBDA=.FALSE.
         
         EPS1=ABS(LAMBDA_MIN-MINVAL(LAMBDA))
         EPS2=ABS(LAMBDA_MAX-MAXVAL(LAMBDA))
         IF ((EPS1.LE.1.0E-08).AND.(EPS2.LE.1.0E-08)) THEN       
             L_LAMBDA=.TRUE.
         ENDIF
         
! If the limits agree then interpolate filter function onto regular 
! grid. This step is required in order to do the integration. If the
! linits do not agree the weights are all set to 1.

! C1 is used as a dummy variable     
      
         IF (L_LAMBDA) THEN

            DO J=1,FILTER%N_PTS
               C1=(LAMBDA_MIN+(J-1)*D_LAMBDA)*1.0E+06
               CALL LINEAR_INTERPOLATION(LAMBDA*1.0E+06,
     &                FILTER%RESPONSE,
     &                FILTER%N_PTS,C1,
     &                WEIGHTS(FILTER%N_PTS-J+1))  
            ENDDO
   
!            DO J=1,FILTER%N_PTS
!               WRITE(20,*) WEIGHTS(J),
!     &                     FILTER%RESPONSE(FILTER%N_PTS-J+1)
!            ENDDO
  
         ELSE
            WRITE(IU_STDOUT, '(A)') ' WAVELENGTH LIMITS IN THE '
     &        // 'SPECTRAL AND IN THE FILTER FILE DO NOT AGREE!'
            WRITE(IU_STDOUT, '(A)') ' ALL WEIGHTS HAVE BEEN SET TO 1'
            WEIGHTS=1.0   
         ENDIF
      
! Interpolate the solar spectrum onto a regular grid with the same d_lambda
! as the phase function

         SOLAR_LAMBDA_MIN=SOLAR_WAVELENGTH(1)
         SOLAR_LAMBDA_MAX=SOLAR_WAVELENGTH(N_SOLAR_POINTS)
      
         N_SOLAR_INTERPOLATED=(SOLAR_LAMBDA_MAX-
     &             SOLAR_LAMBDA_MIN)/D_LAMBDA
     
         ALLOCATE(IRRAD_INTERPOL(0:N_SOLAR_INTERPOLATED))
         ALLOCATE(LAMBDA_INTERPOL(0:N_SOLAR_INTERPOLATED))
      
         DO J=0,N_SOLAR_INTERPOLATED
            LAMBDA_INTERPOL(J)=SOLAR_LAMBDA_MIN+J*D_LAMBDA
            CALL LINEAR_INTERPOLATION(SOLAR_WAVELENGTH
     &          , SOLAR_IRRAD
     &          , N_SOLAR_POINTS
     &          , LAMBDA_INTERPOL(J)
     &          , IRRAD_INTERPOL(J))
!            WRITE(21,*) J, LAMBDA_INTERPOL(J),IRRAD_INTERPOL(J)
         ENDDO
      
! Calculate total energy input by integrating over the whole
! solar sectrum

         TOTAL_SOLAR_IRRADIANCE=0.0
         CALL SIMPSONS_RULE(N_SOLAR_INTERPOLATED,D_LAMBDA
     &    ,IRRAD_INTERPOL, TOTAL_SOLAR_IRRADIANCE)
 
! Add on the tail of the distribution using a rayleigh-jeans law
       
         IRRADIANCE_TAIL=RAYLEIGH_JEANS_TAIL(
     &                    SOLAR_WAVELENGTH(N_SOLAR_POINTS))
     
         TOTAL_SOLAR_IRRADIANCE=TOTAL_SOLAR_IRRADIANCE
     &                         +IRRADIANCE_TAIL

! Find energy in the required band

         CALL INNER_BRACKET(IERR
     &       , LAMBDA_MIN,LAMBDA_MAX
     &       , N_SOLAR_INTERPOLATED
     &       , LAMBDA_INTERPOL
     &       , I_SHORT,I_LONG)
     
! Integrate across the band and normalize by 
! the total irradiance.   

         DO J=1,FILTER%N_PTS
            WEIGHTS(J)=IRRAD_INTERPOL(I_SHORT-1+J)*WEIGHTS(J)
         ENDDO  
          
         BAND_IRRADIANCE=0.0
         CALL SIMPSONS_RULE(FILTER%N_PTS,D_LAMBDA,WEIGHTS,
     &                      BAND_IRRADIANCE)

!        PRINT*,'TOTAL ENERGY=',TOTAL_SOLAR_IRRADIANCE
!        PRINT*,'ENERGY',BAND_IRRADIANCE
!        PRINT*,'RATIO',BAND_IRRADIANCE/TOTAL_SOLAR_IRRADIANCE

         SOLAR_FLUX_BAND(I)=BAND_IRRADIANCE/TOTAL_SOLAR_IRRADIANCE

      ENDDO
      
! Add the tail to the last band if required. the setting of limits
! above dealt with the contributions within the region of solar
! data.

      IF (L_ENHANCE) THEN
        I_LAST_BAND=1
        DO I=1, N_BAND
          IF (WAVE_LENGTH_LONG(I) > WAVE_LENGTH_SHORT(I_LAST_BAND))
     &      I_LAST_BAND=I
        ENDDO
        SOLAR_FLUX_BAND(I_LAST_BAND)=SOLAR_FLUX_BAND(I_LAST_BAND)
     &    +IRRADIANCE_TAIL/TOTAL_SOLAR_IRRADIANCE
      ENDIF
      
! Reduce the solar flux in bands which contain exclusions if
! necessary.

      IF (L_EXCLUDE) THEN
        DO I=1, N_BAND
          DO J=1, N_BAND_EXCLUDE(I)
            SOLAR_FLUX_BAND(I)=SOLAR_FLUX_BAND(I)
     &        -SOLAR_FLUX_BAND(INDEX_EXCLUDE(J, I))
          ENDDO
        ENDDO
      ENDIF      
      
      DEALLOCATE(LAMBDA)  
      DEALLOCATE(WEIGHTS)
      
      L_PRESENT_2=.TRUE.

      RETURN
      END
      
      
