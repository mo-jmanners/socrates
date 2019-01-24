! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the variables in the Socrates boundary conditions type
!
!------------------------------------------------------------------------------
module socrates_set_bound
implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_SET_BOUND'
contains

subroutine set_bound(bound, control, dimen, spectrum, n_profile, &
  t_ground, cos_zenith_angle, solar_irrad, orog_corr, &
  l_grey_albedo, grey_albedo, &
  l_debug, i_profile_debug)

use def_bound,    only: StrBound, allocate_bound
use def_control,  only: StrCtrl
use def_dimen,    only: StrDim
use def_spectrum, only: StrSpecData
use realtype_rd,  only: RealK
use rad_pcf,      only: &
  ip_solar, ip_infra_red, ip_surf_alb_diff, ip_surf_alb_dir

implicit none


! Boundary properties:
type(StrBound),    intent(out) :: bound

! Control options:
type(StrCtrl),     intent(in)  :: control

! Dimensions:
type(StrDim),      intent(in)  :: dimen

! Spectral data:
type(StrSpecData), intent(in)  :: spectrum

integer, intent(in) :: n_profile
!   Number of atmospheric profiles for radiation calculations

real(RealK), intent(in), optional :: t_ground(n_profile)
!   Effective radiative temperature over whole grid-box
real(RealK), intent(in), optional :: cos_zenith_angle(n_profile)
!   Cosine of solar zenith angle
real(RealK), intent(in), optional :: solar_irrad(n_profile)
!   Solar irradiance at top-of-atmosphere (mean over timestep)
real(RealK), intent(in), optional :: orog_corr(n_profile)
!   Orographic correction factor

logical, intent(in), optional :: l_grey_albedo
!   Set a single grey albedo for the surface
real(RealK), intent(in), optional :: grey_albedo
!   Grey surface albedo

logical, intent(in), optional :: l_debug
integer, intent(in), optional :: i_profile_debug
!   Options for outputting debugging information

! Local variables.
integer :: l, i


! Allocate structure for the core radiation code interface
call allocate_bound(bound, dimen, spectrum)

! Surface albedo
if (l_grey_albedo) then
  bound%rho_alb(1:n_profile, ip_surf_alb_diff, 1:spectrum%basic%n_band) &
    = grey_albedo
  bound%rho_alb(1:n_profile, ip_surf_alb_dir,  1:spectrum%basic%n_band) &
    = grey_albedo
else
  bound%rho_alb(1:n_profile, ip_surf_alb_diff, 1:spectrum%basic%n_band) = 0.0
  bound%rho_alb(1:n_profile, ip_surf_alb_dir,  1:spectrum%basic%n_band) = 0.0
end if

! Set the surface basis functions for a Lambertian surface.
bound%n_brdf_basis_fnc=1
! By defining F_{1,0,0,0} to be 4, rho_alb becomes equal to the diffuse albedo.
bound%f_brdf(1, 0, 0, 0)=4.0


! Surface temperature
if (present(t_ground)) then
  bound%t_ground(1:n_profile) = t_ground(1:n_profile)
end if


if (control%l_tile) then
  ! Set up the surface tiling variables (none for the minute).
  bound%n_tile=3
  bound%n_point_tile=0
end if


! Incident solar flux
if (present(cos_zenith_angle) .and. present(solar_irrad)) then
  do l=1, n_profile
    if (cos_zenith_angle(l) > 0.0) then
      bound%solar_irrad(l)=solar_irrad(l)
      bound%zen_0(l)=1.0/cos_zenith_angle(l)
    else
      bound%solar_irrad(l)=0.0
      bound%zen_0(l)=1.0
    end if
  end do
end if

! Orographic correction factor
if (present(orog_corr)) then
  bound%orog_corr(1:n_profile)=orog_corr(1:n_profile)
end if


if (present(l_debug)) then
  if (l_debug) then
    if (present(i_profile_debug)) then
      l = i_profile_debug
    else
      l = 1
    end if
    write(2000+l,'(A)') 'TEMP(K) SOLAR_IRRAD(WM-2) ZEN_0 OROG_CORR'
    write(2000+l,'(4(1pe16.8))') bound%t_ground(l), &
      bound%solar_irrad(l), bound%zen_0(l), bound%orog_corr(l)
    write(2000+l,'(A)') 'BAND DIFFUSE_ALBEDO DIRECT_ALBEDO'
    do i=1, spectrum%basic%n_band
      write(2000+l,'(i0, 2(1pe16.8))') i, &
        bound%rho_alb(l, ip_surf_alb_diff, i), &
        bound%rho_alb(l, ip_surf_alb_dir, i)
    end do     
  end if
end if


end subroutine set_bound
end module socrates_set_bound
