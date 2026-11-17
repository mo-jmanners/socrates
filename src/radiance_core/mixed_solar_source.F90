! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to set the solar source terms in a mixed column.
!
! Method:
!   The direct beam is calculated by propagating down through
!   the column. These direct fluxes are used to `define' the
!   source terms in each layer.
!
!- ---------------------------------------------------------------------
MODULE mixed_solar_source_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'MIXED_SOLAR_SOURCE_MOD'
CONTAINS
SUBROUTINE mixed_solar_source(control, bound                            &
    , n_profile, n_layer, n_cloud_top                                   &
    , flux_inc_direct                                                   &
    , l_scale_solar, adjust_solar_ke                                    &
    , trans_0_free_dir, trans_0_free, source_coeff_free                 &
    , g_ff, g_fc, g_cf, g_cc                                            &
    , trans_0_cloud_dir, trans_0_cloud, source_coeff_cloud              &
    , diffuse_albedo, direct_albedo                                     &
    , flux_direct, flux_direct_dir                                      &
    , source_ground_cloud, source_ground_free                           &
    , s_up_free, s_down_free                                            &
    , s_up_cloud, s_down_cloud                                          &
    , nd_profile, nd_layer, id_ct, nd_source_coeff                      &
    )


  USE realtype_rd, ONLY: RealK
  USE def_control, ONLY: StrCtrl
  USE def_bound,   ONLY: StrBound
  USE rad_pcf, ONLY: ip_direct_csr_scaling, ip_direct_noscaling,        &
                     ip_scf_solar_down, ip_scf_solar_up
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl),      INTENT(IN)    :: control

! Boundary conditions:
  TYPE(StrBound),     INTENT(IN)    :: bound

! Sizes of dummy arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , id_ct                                                             &
!       Topmost declared cloudy layer
    , nd_source_coeff
!       Size allocated for coefficients in the source function


! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , n_cloud_top
!       Top cloudy layer

! Special arrays for equivalent extinction:
  LOGICAL, INTENT(IN) ::                                                &
      l_scale_solar
!       Scaling applied to solar flux
  REAL (RealK), INTENT(IN) ::                                           &
       adjust_solar_ke(nd_profile, nd_layer)
!       Adjustment to solar fluxes with equivalent extinction

  REAL (RealK), INTENT(IN) ::                                           &
      flux_inc_direct(nd_profile)
!       Incident direct solar flux

! Clear-sky optical properties:
  REAL (RealK), INTENT(IN) ::                                           &
      trans_0_free(nd_profile, nd_layer)                                &
!       Free direct transmission
    , trans_0_free_dir(nd_profile, nd_layer)                            &
!       Free unscaled direct transmission
    , source_coeff_free(nd_profile, nd_layer, nd_source_coeff)
!       Clear-sky source coefficients

! Cloudy optical properties:
  REAL (RealK), INTENT(IN) ::                                           &
      trans_0_cloud(nd_profile, nd_layer)                               &
!       Cloudy direct transmission
    , trans_0_cloud_dir(nd_profile, nd_layer)                           &
!       Cloudy unscaled direct transmission
    , source_coeff_cloud(nd_profile, nd_layer, nd_source_coeff)
!       Cloudy reflectance

! Energy transfer coefficients:
  REAL (RealK), INTENT(IN) ::                                           &
      g_ff(nd_profile, id_ct-1: nd_layer)                               &
!       Energy transfer coefficient
    , g_fc(nd_profile, id_ct-1: nd_layer)                               &
!       Energy transfer coefficient
    , g_cf(nd_profile, id_ct-1: nd_layer)                               &
!       Energy transfer coefficient
    , g_cc(nd_profile, id_ct-1: nd_layer)
!       Energy transfer coefficient

! Albedos
  REAL (RealK), INTENT(IN) ::                                           &
      diffuse_albedo(nd_profile)                                        &
!       Diffuse albedo
    , direct_albedo(nd_profile)
!       Direct albedo

! Calculated direct flux and source terms:
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct(nd_profile, 0: nd_layer)                              &
!       Direct flux
    , flux_direct_dir(nd_profile, 0: nd_layer)                          &
!       Direct flux using direct tau
    , source_ground_cloud(nd_profile)                                   &
!       Cloudy source flux at ground
    , source_ground_free(nd_profile)                                    &
!       Free source flux at ground
    , s_up_free(nd_profile, nd_layer)                                   &
!       Free upward source function
    , s_down_free(nd_profile, nd_layer)                                 &
!       Free downward source function
    , s_up_cloud(nd_profile, nd_layer)                                  &
!       Cloudy upward source function
    , s_down_cloud(nd_profile, nd_layer)
!       Cloudy downward source function


! Local variables.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l
!       Loop variable

  REAL (RealK) ::                                                       &
      solar_top_free(nd_profile)                                        &
!       Free solar flux at top of layer
    , solar_top_cloud(nd_profile)                                       &
!       Cloudy solar flux at top of layer
    , solar_base_free(nd_profile)                                       &
!       Free solar flux at base of layer
    , solar_base_cloud(nd_profile)
!       Cloudy solar flux at base of layer

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='MIXED_SOLAR_SOURCE'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! The clear and cloudy direct fluxes are calculated separately
! and added together to form the total direct flux.

! Set incident fluxes.
  DO l=1, n_profile
    flux_direct(l, 0)=flux_inc_direct(l)
  END DO
  IF (control%i_direct_tau == ip_direct_noscaling .OR.                  &
      control%i_direct_tau == ip_direct_csr_scaling) THEN
    DO l=1, n_profile
      flux_direct_dir(l, 0)=flux_inc_direct(l)
    END DO
  END IF


! With equivalent extinction the direct solar flux must be corrected.
  IF (l_scale_solar) THEN

    IF (control%l_diffuse_2stream) THEN
      DO i=1, n_cloud_top-1
        DO l=1, n_profile
          flux_direct(l, i)                                             &
            =flux_direct(l, i-1)*trans_0_free(l, i)                     &
            *adjust_solar_ke(l, i)
          s_up_free(l, i)=source_coeff_free(l, i, ip_scf_solar_up)      &
            *flux_direct(l, i-1)
          s_down_free(l, i)=source_coeff_free(l, i, ip_scf_solar_down)  &
            *flux_direct(l, i-1)
        END DO
      END DO
    ELSE
      DO i=1, n_cloud_top-1
        DO l=1, n_profile
          flux_direct(l, i)                                             &
            =flux_direct(l, i-1)*trans_0_free(l, i)                     &
            *adjust_solar_ke(l, i)
          s_up_free(l, i)=source_coeff_free(l, i, ip_scf_solar_up)      &
            *flux_direct(l, i-1)
          s_down_free(l, i)                                             &
            =(source_coeff_free(l, i, ip_scf_solar_down)                &
            -trans_0_free(l, i))*flux_direct(l, i-1)                    &
            +flux_direct(l, i)
        END DO
      END DO
    END IF

    IF (control%i_direct_tau == ip_direct_noscaling .OR.                &
        control%i_direct_tau == ip_direct_csr_scaling) THEN
      DO i=1, n_cloud_top-1
        DO l=1, n_profile
          flux_direct_dir(l, i)                                         &
            =flux_direct_dir(l, i-1)*trans_0_free_dir(l, i)             & 
             *adjust_solar_ke(l, i)
        END DO
      END DO
    END IF

  ELSE

    DO i=1, n_cloud_top-1
      DO l=1, n_profile
        flux_direct(l, i)                                               &
          =flux_direct(l, i-1)*trans_0_free(l, i)
        s_up_free(l, i)=source_coeff_free(l, i, ip_scf_solar_up)        &
          *flux_direct(l, i-1)
        s_down_free(l, i)                                               &
          =source_coeff_free(l, i, ip_scf_solar_down)                   &
          *flux_direct(l, i-1)
      END DO
    END DO
    IF (control%i_direct_tau == ip_direct_noscaling .OR.                &
        control%i_direct_tau == ip_direct_csr_scaling) THEN
      DO i=1, n_cloud_top-1
        DO l=1, n_profile
          flux_direct_dir(l, i)                                         &
          =flux_direct_dir(l, i-1)*trans_0_free_dir(l, i)
        END DO
      END DO
    END IF

  END IF



! Clear and cloudy region.
! Initialize partial fluxes:
  DO l=1, n_profile
    solar_base_free(l)=flux_direct(l, n_cloud_top-1)
    solar_base_cloud(l)=0.0e+00_RealK
  END DO


  DO i=n_cloud_top, n_layer

!   Transfer fluxes across the interface. The use of only one
!   cloudy flux implicitly forces random overlap of different
!   subclouds within the cloudy parts of the layer.

    DO l=1, n_profile
      solar_top_cloud(l)=g_cc(l, i-1)*solar_base_cloud(l)               &
        +g_fc(l, i-1)*solar_base_free(l)
      solar_top_free(l)=g_ff(l, i-1)*solar_base_free(l)                 &
        +g_cf(l, i-1)*solar_base_cloud(l)
    END DO


!   Propagate the clear and cloudy fluxes through the layer:
    IF (l_scale_solar) THEN

      IF (control%l_diffuse_2stream) THEN
        DO l=1, n_profile
          solar_base_free(l)=solar_top_free(l)                          &
            *trans_0_free(l, i)*adjust_solar_ke(l, i)
          solar_base_cloud(l)=solar_top_cloud(l)                        &
            *trans_0_cloud(l, i)*adjust_solar_ke(l, i)
          s_up_free(l, i)=source_coeff_free(l, i, ip_scf_solar_up)      &
            *solar_top_free(l)
          s_down_free(l, i)=source_coeff_free(l, i, ip_scf_solar_down)  &
            *solar_top_free(l)
          s_up_cloud(l, i)=source_coeff_cloud(l, i, ip_scf_solar_up)    &
            *solar_top_cloud(l)
          s_down_cloud(l, i)=source_coeff_cloud(l, i, ip_scf_solar_down)&
            *solar_top_cloud(l)
        END DO
      ELSE
        DO l=1, n_profile
          solar_base_free(l)=solar_top_free(l)                          &
            *trans_0_free(l, i)*adjust_solar_ke(l, i)
          solar_base_cloud(l)=solar_top_cloud(l)                        &
            *trans_0_cloud(l, i)*adjust_solar_ke(l, i)
          s_up_free(l, i)=source_coeff_free(l, i, ip_scf_solar_up)      &
            *solar_top_free(l)
          s_down_free(l, i)                                             &
            =(source_coeff_free(l, i, ip_scf_solar_down)                &
            -trans_0_free(l, i))*solar_top_free(l)                      &
            +solar_base_free(l)
          s_up_cloud(l, i)                                              &
            =source_coeff_cloud(l, i, ip_scf_solar_up)                  &
            *solar_top_cloud(l)
          s_down_cloud(l, i)                                            &
            =(source_coeff_cloud(l, i, ip_scf_solar_down)               &
            -trans_0_cloud(l, i))*solar_top_cloud(l)                    &
            +solar_base_cloud(l)
        END DO
      END IF

    ELSE

      DO l=1, n_profile
        solar_base_free(l)=solar_top_free(l)                            &
          *trans_0_free(l, i)
        solar_base_cloud(l)=solar_top_cloud(l)                          &
          *trans_0_cloud(l, i)
        s_up_free(l, i)=source_coeff_free(l, i, ip_scf_solar_up)        &
          *solar_top_free(l)
        s_down_free(l, i)                                               &
          =source_coeff_free(l, i, ip_scf_solar_down)                   &
          *solar_top_free(l)
        s_up_cloud(l, i)                                                &
          =source_coeff_cloud(l, i, ip_scf_solar_up)                    &
          *solar_top_cloud(l)
        s_down_cloud(l, i)                                              &
          =source_coeff_cloud(l, i, ip_scf_solar_down)                  &
          *solar_top_cloud(l)
      END DO

    END IF


!   Calculate the total direct flux.
    DO l=1, n_profile
      flux_direct(l, i)=solar_base_free(l)+solar_base_cloud(l)
    END DO

  END DO


! Repeat for direct flux using trans without scaling
  IF (control%i_direct_tau == ip_direct_noscaling .OR.                  &
      control%i_direct_tau == ip_direct_csr_scaling) THEN

    IF (control%l_diffuse_2stream) THEN
      ! Save the last value at the base of the cloud and clear-sky.
      DO l=1, n_profile
        source_ground_cloud(l)=solar_base_cloud(l)
        source_ground_free(l)=solar_base_free(l)
      END DO
    END IF

!   Clear and cloudy region.
!   Initialize partial fluxes:
    DO l=1, n_profile
      solar_base_free(l)=flux_direct_dir(l, n_cloud_top-1)
      solar_base_cloud(l)=0.0e+00_RealK
    END DO


    DO i=n_cloud_top, n_layer

!     Transfer fluxes across the interface. The use of only one
!     cloudy flux implicitly forces random overlap of different
!     subclouds within the cloudy parts of the layer.
      DO l=1, n_profile
        solar_top_cloud(l)=g_cc(l, i-1)*solar_base_cloud(l)             &
          +g_fc(l, i-1)*solar_base_free(l)
        solar_top_free(l)=g_ff(l, i-1)*solar_base_free(l)               &
          +g_cf(l, i-1)*solar_base_cloud(l)
      END DO


!     Propagate the clear and cloudy fluxes through the layer:
      IF (l_scale_solar) THEN

        DO l=1, n_profile
          solar_base_free(l)=solar_top_free(l)                          &
            *trans_0_free_dir(l, i)*adjust_solar_ke(l, i)
          solar_base_cloud(l)=solar_top_cloud(l)                        &
            *trans_0_cloud_dir(l, i)*adjust_solar_ke(l, i)
        END DO

      ELSE

        DO l=1, n_profile
          solar_base_free(l)=solar_top_free(l)                          &
            *trans_0_free_dir(l, i)
          solar_base_cloud(l)=solar_top_cloud(l)                        &
            *trans_0_cloud_dir(l, i)
        END DO

      END IF


!     Calculate the total direct flux.
      DO l=1, n_profile
        flux_direct_dir(l, i)=solar_base_free(l)+solar_base_cloud(l)
      END DO

    END DO

    IF (control%l_diffuse_2stream) THEN
      DO l=1, n_profile
        source_ground_cloud(l)=solar_base_cloud(l)*direct_albedo(l)     &
          -(solar_base_cloud(l)-source_ground_cloud(l))                 &
          *diffuse_albedo(l)
        source_ground_free(l)=solar_base_free(l)*direct_albedo(l)       &
          -(solar_base_free(l)-source_ground_free(l))                   &
          *diffuse_albedo(l)
      END DO
    END IF

  ELSE IF (control%l_diffuse_2stream) THEN

    DO l=1, n_profile
      source_ground_cloud(l)=solar_base_cloud(l)*direct_albedo(l)
      source_ground_free(l)=solar_base_free(l)*direct_albedo(l)
    END DO

  END IF

  IF (.NOT.control%l_diffuse_2stream) THEN
    DO l=1, n_profile
      source_ground_cloud(l)=(direct_albedo(l)-diffuse_albedo(l))       &
        *solar_base_cloud(l)
      source_ground_free(l)=(direct_albedo(l)-diffuse_albedo(l))        &
        *solar_base_free(l)
    END DO
  END IF


! Correct the direct flux at the ground for sloping terrain
  IF (control%l_orog) THEN
    IF (control%i_direct_tau == ip_direct_noscaling .OR. &
        control%i_direct_tau == ip_direct_csr_scaling) THEN
      IF (control%l_diffuse_2stream) THEN
        DO l=1, n_profile
          flux_direct(l, n_layer) = flux_direct(l, n_layer) &
            + flux_direct_dir(l, n_layer) &
            * (bound%orog_corr(l) - 1.0_RealK)
        END DO
      END IF
      DO l=1, n_profile
        flux_direct_dir(l, n_layer) = flux_direct_dir(l, n_layer) &
          * bound%orog_corr(l)
      END DO
    ELSE
      DO l=1, n_profile
        flux_direct(l, n_layer) = flux_direct(l, n_layer) &
          * bound%orog_corr(l)
      END DO
    END IF

    IF (control%l_diffuse_2stream) THEN
      DO l=1, n_profile
        source_ground_cloud(l) = source_ground_cloud(l) &
          + solar_base_cloud(l) * direct_albedo(l) &
          * (bound%orog_corr(l) - 1.0_RealK)

        source_ground_free(l) = source_ground_free(l) &
          + solar_base_free(l) * direct_albedo(l) &
          * (bound%orog_corr(l) - 1.0_RealK)
      END DO
    ELSE
      DO l=1, n_profile
        source_ground_cloud(l) = source_ground_cloud(l) &
          + solar_base_cloud(l) * (direct_albedo(l)-diffuse_albedo(l)) &
          * (bound%orog_corr(l) - 1.0_RealK)

        source_ground_free(l) = source_ground_free(l) &
          + solar_base_free(l) * (direct_albedo(l)-diffuse_albedo(l)) &
          * (bound%orog_corr(l) - 1.0_RealK)

        s_down_free(l, n_layer) = s_down_free(l, n_layer) &
          + solar_base_free(l) &
          * (bound%orog_corr(l) - 1.0_RealK)
       
        s_down_cloud(l, n_layer) = s_down_cloud(l, n_layer) &
          + solar_base_cloud(l) &
          * (bound%orog_corr(l) - 1.0_RealK)
      END DO
    END IF
  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE mixed_solar_source
END MODULE mixed_solar_source_mod
