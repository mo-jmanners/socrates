! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to calculate the Rayleigh-Jeans approximation.
!
! Method:
!	The standard Rayleigh-Jeans formula is used.
!
! Note:
!	This routine is used to provide a tail to the Planck function.
!
!- ---------------------------------------------------------------------
      FUNCTION rayleigh_jeans_tail(lambda)
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE rad_ccf

      IMPLICIT NONE


!     Dummy variables.
      REAL  (RealK), Intent(IN) ::
     &    lambda
!           Initial wavelength
      REAL  (RealK) ::
     &    rayleigh_jeans_tail
!           Tail irradiance.
!
!
!
!     Evaluate black body flux.
      rayleigh_jeans_tail=2.0_RealK*c_light*k_boltzmann
     &   *t_effective_solar/(3.0_RealK*lambda**3)
!     Scale to get the flux at the top of the atmosphere.
      rayleigh_jeans_tail=rayleigh_jeans_tail
     &   *(sun_radius/d_earth_sun)**2
!
!
!
      RETURN
      END
