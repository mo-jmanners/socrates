! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read a solar spectrum.
!
! Method:
!	The solar spetrum is read a a file of irradiances against
!	wavelengths.
!
!- ---------------------------------------------------------------------
      SUBROUTINE read_solar_spectrum(ierr, l_solar_spectrum
     &  , n_solar_points, solar_wavelength, solar_irrad)

      USE realtype_rd
      USE rad_pcf
      USE dimensions_pp_ucf
      USE def_std_io_icf

      IMPLICIT NONE


!     Dummy arguments.
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
     &  , n_solar_points
!           Number of solar data points
      LOGICAL, Intent(OUT) ::
     &    l_solar_spectrum
!           Solar spectral flag
      REAL  (RealK), Intent(OUT) ::
     &    solar_wavelength(npd_solar_points)
!           Wavelengths for irradiance
     &  , solar_irrad(npd_solar_points)
!           Solar irradiance at toa
!
!     Local variables.
      CHARACTER
     &    line*80
!           Line of input data
      INTEGER
     &    ios
!           I/O error flag
     &  , i
!           Loop variable
      LOGICAL
     &    l_data_region
!           Flag for region of data within file.
!
      INTEGER :: iu_solar
!           Unit number for the solar spectrum file

      REAL (RealK) :: scale_wv, scale_irr


!     obtain the file containing the solar spectrum.
      CALL get_free_unit(ierr, iu_solar)
      CALL open_file_in(ierr, iu_solar 
     &  , 'enter the name of the file containing '
     &   //'the solar irradiance data.'
     &  )
      IF (ierr /= i_normal) RETURN
!
!     Read in the file.
      l_data_region=.false.
      scale_wv=1.0_RealK
      scale_irr=1.0_RealK
      i=0
3     read(iu_solar, '(a)', end=4) line
        IF (l_data_region) THEN
          IF (line(1:4) /= '*END') THEN
            backspace(iu_solar)
            i=i+1
            IF (i >= npd_solar_points) THEN
              WRITE(iu_err, '(a)')
     &          '*** error: there are too many points in the'
     &          //'solar_spectrum.'
              ierr=i_err_fatal
              RETURN
            ENDIF
            READ(iu_solar, *, iostat=ios)
     &        solar_wavelength(i), solar_irrad(i)
            IF (ios /= 0) THEN
              WRITE(iu_err, '(a)')
     &          '*** error: the solar_spectrum could not be read.'
              ierr=i_err_io
              RETURN
            ENDIF
          ELSE
             l_data_region=.false.
          ENDIF
        ELSE
          IF (line(1:11) == '*BEGIN_DATA') THEN
            l_data_region=.true.
          ENDIF
          IF (line(1:11) == '*SCALE_DATA') THEN
            READ(iu_solar, *, iostat=ios) scale_wv, scale_irr
          ENDIF
        ENDIF
!
      goto 3
!
4     close(iu_solar)
      n_solar_points=i

!     Scale the values to the correct units:
      solar_wavelength=solar_wavelength*scale_wv
      solar_irrad=solar_irrad*scale_irr

!     Check that reading was carried out: failure to read the data
!     will cause n_solar_points to be 0.
      IF (n_solar_points == 0) THEN
        WRITE(iu_err, '(/a, /a)') '*** error: no data were read.'
     &    , 'check format of file of irradiances.'
      ENDIF
!
!     Now the spectrum has been read, set the reading flag.
      l_solar_spectrum=.true.
!
!
!
      RETURN
      END
