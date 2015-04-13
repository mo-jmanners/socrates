! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to calculate the solar irradiance at a given wavelength.
!
! Method:
!	The wavelength is bracketed between two points of the solar
!	spectrum and linear interpolation is used in the interval.
!	At long wavelngths a Rayleigh-Jeans tail is used.
!
!- ---------------------------------------------------------------------
      FUNCTION solar_intensity(lambda
     &  , n_solar_points, solar_wavelength, solar_irrad
     &  )

      USE realtype_rd
      USE dimensions_pp_ucf
      USE rad_ccf, ONLY: t_effective_solar, d_earth_sun, sun_radius

      IMPLICIT NONE


!     Dummy arguments
!
      REAL  (RealK) ::
     &    solar_intensity
!           Returned intensity
!
      INTEGER, Intent(IN) ::
     &    n_solar_points
!           Number of spectral points
      REAL  (RealK), Intent(IN) ::
     &    lambda
!           Wavelength
     &  , solar_wavelength(npd_solar_points)
!           Wavelengths of spectrum
     &  , solar_irrad(npd_solar_points)
!           Solar irradiance at toa
!
!     Local variables.
      INTEGER
     &    i_short
!           Array point just below lambda
     &  , i_long
!           Array point just above lambda
      REAL  (RealK) ::
     &    fraction
!           Fraction of interval covered
!
!     Subroutines called:
      EXTERNAL
     &    point_bracket
!
!     Functions called:
      REAL  (RealK) ::
     &    planck
!           Planckian function
      EXTERNAL
     &    planck
!
!
!
!     Find the points of the wavelength array bracketting lambda.
      CALL point_bracket(lambda, n_solar_points, solar_wavelength
     &  , i_short, i_long)
!
!     If lambda lies within the range covered we interpolate: if not
!     use a simple extrapolation
      IF (i_short == 0) THEN
!       Lambda is at a wavelength shorter than any in the spectrum.
        solar_intensity=0.0_RealK
      ELSE IF (i_long >  n_solar_points) THEN
!       Use a black body fit at the effective temperature.
        solar_intensity=planck(t_effective_solar, lambda)
     &    *(sun_radius/d_earth_sun)**2
      ELSE
!       Within the spectrum use linear interpolation.
        fraction=(lambda-solar_wavelength(i_short))
     &    /(solar_wavelength(i_long)-solar_wavelength(i_short))
        solar_intensity=fraction*solar_irrad(i_long)
     &    +(1.0_RealK-fraction)*solar_irrad(i_short)
      ENDIF
!
!
!
      RETURN
      END
