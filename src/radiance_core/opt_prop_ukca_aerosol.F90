! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to account for the optical properties of UKCA aerosols.
!
! Method:
!   UKCA aerosols are only provided in a Henyey-Greenstein
!   framework, where the phase function is only described by
!   its asymmetry parameter.
!   Atmospheric optical properties are updated by accounting
!   for the optical properties of UKCA aerosols which have
!   been computed by the UKCA_RADAER routines.
!   The routine is based on opt_prop_aerosol.F90 where only
!   the Henyey-Greenstein case is retained.
!
!- ---------------------------------------------------------------------
MODULE opt_prop_ukca_aerosol_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'OPT_PROP_UKCA_AEROSOL_MOD'
CONTAINS
SUBROUTINE opt_prop_ukca_aerosol(control, aer, aer_band                 &
    , n_profile, first_layer, last_layer                                &
    , n_order_phase                                                     &
    , k_ext_tot, k_ext_scat, phase_fnc, forward_scatter                 &
    , nd_profile, id_lt, id_lb                                          &
    , nd_max_order                                                      &
    )


  USE realtype_rd, ONLY: RealK
  USE def_control, ONLY: StrCtrl
  USE def_aer,     ONLY: StrAer

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl), INTENT(IN) :: control

! Aerosol properties:
  TYPE(StrAer),  INTENT(IN) :: aer

  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for profiles
    , id_lt                                                             &
!       Topmost declared layer for output optical properties
    , id_lb                                                             &
!       Bottom declared layer for output optical properties
    , nd_max_order
!       Size allocated for orders of sperical harmonics

! Dummy variables.
  INTEGER, INTENT(IN) ::                                                &
      aer_band                                                          &
!       Aerosol band
    , n_profile                                                         &
!       Number of profiles
    , first_layer                                                       &
!       First layer where propeties are required
    , last_layer                                                        &
!       Last layer where propeties are required
    , n_order_phase
!       Number of terms to retain in the phase function

  REAL(RealK), INTENT(INOUT) ::                                         &
      k_ext_scat(nd_profile, id_lt: id_lb)                              &
!       Scattering extinction
    , k_ext_tot(nd_profile, id_lt: id_lb)                               &
!       Total extinction
    , phase_fnc(nd_profile, id_lt: id_lb, nd_max_order)                 &
!       Phase function
    , forward_scatter(nd_profile, id_lt: id_lb)
!       Forward scattering

! Local variables.
  INTEGER :: l, j, i, ls
!       Loop variables

  REAL(RealK) ::                                                        &
      k_scatter(nd_profile)                                             &
!       Scattering of current extinction of the current aerosol
    , asymmetry(nd_profile)                                             &
!       Asymmetry of the current aerosol
    , ks_phf(nd_profile)
!       Scattering coefficient multiplied by a coefficient in the
!       phase function

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='OPT_PROP_UKCA_AEROSOL'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  IF (control%l_aerosol_combined) THEN
    DO j=1, aer%n_mode
      DO i=first_layer, last_layer
        DO l=1, n_profile
          k_ext_tot(l, i) = k_ext_tot(l, i) &
            + aer%mode_mix_ratio(l, i, j) &
            * aer%mode_absorption(l, i, j, aer_band)
          k_ext_scat(l, i) = k_ext_scat(l, i) &
            + aer%mode_mix_ratio(l, i, j) &
            * aer%mode_scattering(l, i, j, aer_band)
          phase_fnc(l, i, 1) = phase_fnc(l, i, 1) &
            + aer%mode_mix_ratio(l, i, j) &
            * aer%mode_asymmetry(l, i, j, aer_band)
        END DO
      END DO
    END DO
    IF (control%l_rescale) THEN
      DO j=1, aer%n_mode
        DO i=first_layer, last_layer
          DO l=1, n_profile
            forward_scatter(l, i) = forward_scatter(l, i) &
              + aer%mode_mix_ratio(l, i, j) &
              * aer%mode_forward_scatter(l, i, j, aer_band)
          END DO
        END DO
      END DO
    END IF
  ELSE
    DO j=1, aer%n_mode
    
      DO i=first_layer, last_layer
    
        DO l=1, n_profile
    
          k_ext_tot(l, i) = k_ext_tot(l, i) &
                            + aer%mode_mix_ratio(l, i, j) * &
                              aer%mode_absorption(l, i, j, aer_band)
          k_scatter(l)    = aer%mode_mix_ratio(l, i, j) * &
                            aer%mode_scattering(l, i, j, aer_band)
          k_ext_scat(l, i) = k_ext_scat(l, i) + k_scatter(l)
    
        END DO ! l
    
    
!       The phase function: Henyey-Greenstein only.
    
    
!       Set the lowest order in the phase function from the asymmetry.
        DO l=1, n_profile
          asymmetry(l) = aer%mode_asymmetry(l, i, j, aer_band)
          phase_fnc(l, i, 1) = phase_fnc(l, i, 1) + &
                               k_scatter(l) * asymmetry(l)
        END DO ! l
    
    
!       Initialise the product of the scattering and the current
!       moment of the phase function.
        IF (control%l_rescale .OR. (n_order_phase >= 2)) THEN
          DO l=1, n_profile
            ks_phf(l) = k_scatter(l) * asymmetry(l)
          END DO ! l
        END IF
    
    
!       Calculate weighted higher moments recursively.
        DO ls=2, n_order_phase
          DO l=1, n_profile
            ks_phf(l) = ks_phf(l) * asymmetry(l)
            phase_fnc(l, i, ls) = phase_fnc(l, i, ls) + ks_phf(l)
          END DO ! l
        END DO ! ls
    
    
!       Usually, we will retain terms as far as the order of
!       truncation, but in the case of two-stream methods the
!       order of truncation will exceed the order of retention
!       by 1.
    
        IF (control%l_rescale) THEN
    
          IF (control%n_order_forward == n_order_phase) THEN
            DO l=1, n_profile
              forward_scatter(l, i) = forward_scatter(l, i) + ks_phf(l)
            END DO ! l
          ELSE IF (control%n_order_forward == n_order_phase + 1) THEN
            DO l=1, n_profile
              forward_scatter(l, i) = forward_scatter(l, i) + &
                                      ks_phf(l) * asymmetry(l)
            END DO ! l
          ELSE
!           This case probably shouldn't arise so we use
!           inefficient explicit exponentiation.
            DO l=1, n_profile
              forward_scatter(l, i) = forward_scatter(l, i) + &
                k_scatter(l) * asymmetry(l)**control%n_order_forward
            END DO ! l
          END IF
    
        END IF
    
      END DO ! i
    
    END DO ! j
  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE opt_prop_ukca_aerosol
END MODULE opt_prop_ukca_aerosol_mod
