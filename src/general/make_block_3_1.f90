! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to make spectral blocks of type 3.
!
! Method:
!   The monochromatic Rayleigh scattering coefficients are calculated and 
!   weighted with the solar spectrum.
!
!- ---------------------------------------------------------------------
SUBROUTINE make_block_3_1(Sp, Sol, Refract, l_h2he_atm)

  USE realtype_rd,   ONLY: RealK
  USE def_spectrum,  ONLY: StrSpecData
  USE def_solarspec, ONLY: StrSolarSpec
  USE def_refract,   ONLY: StrRefract

  IMPLICIT NONE


  TYPE (StrSpecData), INTENT(INOUT) :: Sp
!   Spectral file to be assigned
  TYPE (StrSolarSpec), INTENT(IN) :: Sol
!   Solar spectrum
  TYPE (StrRefract), INTENT(IN) :: Refract
!   Refractive index of gas
  LOGICAL :: l_h2he_atm
!   Calculate rayleigh scattering for H2/He atmosphere

! Local variables
  INTEGER :: n_int
!   Number of intervals for integration
  INTEGER :: i_begin
!   Beginning of spectrum in band
  INTEGER :: i_end
!   End of spectrum in band
  INTEGER :: i_dummy
!   Dummy variable
  INTEGER :: i, k
!   Loop variables
  REAL (RealK) :: wave_int(Sol%n_points+2)
!   Wavelngths for integration
  REAL (RealK) :: weight_int(Sol%n_points+2)
!   Weighting solar irradiances
  REAL (RealK) :: product_int(Sol%n_points+2)
!   Irradiance times rayleigh coef.

! Functions called:
  REAL (RealK), EXTERNAL :: rayleigh_scatter_air
  REAL (RealK), EXTERNAL :: rayleigh_scatter_h2he
!   Functions for rayleigh scattering
  REAL (RealK), EXTERNAL :: solar_intensity
!   Solar intensity at a given wavelength
  REAL (RealK), EXTERNAL :: trapezoid
!   Trapezoidal integration function



! In each band the scattering coefficient is averaged at the points
! of the solar spectrum and the limits of the band.
  DO i=1, Sp%Basic%n_band
    CALL point_bracket(Sp%Basic%wavelength_short(i), &
      Sol%n_points, Sol%wavelength, i_dummy, i_begin)
    CALL point_bracket(Sp%Basic%wavelength_long(i), &
      Sol%n_points, Sol%wavelength, i_end, i_dummy)
!   I_BEGIN and I_END are the points of the solar spectrum
!   just within the band. Form an array of wavelength points for
!   integration:
    n_int=1
    wave_int(n_int)=Sp%Basic%wavelength_short(i)
    weight_int(n_int)=solar_intensity(wave_int(n_int), Sol)
    IF (l_h2he_atm) THEN
      product_int(n_int)=weight_int(n_int) &
        *rayleigh_scatter_h2he(wave_int(n_int), &
        Refract%wavelength, Refract%re_part, Refract%n_points)
      DO k=i_begin, i_end
        n_int=n_int+1
        wave_int(n_int)=Sol%wavelength(k)
        weight_int(n_int)=Sol%irrad(k)
        product_int(n_int)=weight_int(n_int) &
          *rayleigh_scatter_h2he(wave_int(n_int), &
          Refract%wavelength, Refract%re_part, Refract%n_points)
      END DO
      n_int=n_int+1
      wave_int(n_int)=Sp%Basic%wavelength_long(i)
      weight_int(n_int)=solar_intensity(wave_int(n_int), Sol)
      product_int(n_int)=weight_int(n_int) &
        *rayleigh_scatter_h2he(wave_int(n_int), &
        Refract%wavelength, Refract%re_part, Refract%n_points)
    ELSE
      product_int(n_int)=weight_int(n_int) &
        *rayleigh_scatter_air(wave_int(n_int))
      DO k=i_begin, i_end
        n_int=n_int+1
        wave_int(n_int)=Sol%wavelength(k)
        weight_int(n_int)=Sol%irrad(k)
        product_int(n_int)=weight_int(n_int) &
          *rayleigh_scatter_air(wave_int(n_int))
      END DO
      n_int=n_int+1
      wave_int(n_int)=Sp%Basic%wavelength_long(i)
      weight_int(n_int)=solar_intensity(wave_int(n_int), Sol)
      product_int(n_int)=weight_int(n_int) &
        *rayleigh_scatter_air(wave_int(n_int))
    END IF
    Sp%Rayleigh%rayleigh_coeff(i)= trapezoid(n_int, wave_int, product_int)
    Sp%Rayleigh%rayleigh_coeff(i)= &
      Sp%Rayleigh%rayleigh_coeff(i)/trapezoid(n_int, wave_int, weight_int)
  END DO

END SUBROUTINE make_block_3_1
