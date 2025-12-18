program calc_n2o_xsc

! JPL 19-5 recommended cross-section for N2O

! 173nm - 240nm, 194K - 302K: Table 4C-3-3 from JPL 19-5
! Selwyn et al. (1977): Nitrous oxide ultraviolet absorption spectrum
! at stratospheric temperatures. Geophys. Res. Lett. 1977, 4, 427-430,
! doi:10.1029/GL004i010p00427

! 160nm - 173nm, 208K & 298K:
! 240nm - 250nm, 298K:
! Hubrich and Stuhl (1980): Ultraviolet absorption of some halogenated
! methanes and ethanes of atmospheric interest. J. Photochem. 12, 93-107
! doi:10.1016/0047-2670(80)85031-3

implicit none

integer, parameter :: RealK=SELECTED_REAL_KIND(15, 307)
integer :: i, j, io, stat
logical :: exists
character(len=512) :: msg

real(RealK), parameter :: A0 = 68.21023_RealK
real(RealK), parameter :: A1 = -4.071805_RealK
real(RealK), parameter :: A2 = 4.301146e-2_RealK
real(RealK), parameter :: A3 = -1.777846e-4_RealK
real(RealK), parameter :: A4 = 2.520672e-7_RealK

real(RealK), parameter :: B0 = 123.4014_RealK
real(RealK), parameter :: B1 = -2.116255_RealK
real(RealK), parameter :: B2 = 1.111572e-2_RealK
real(RealK), parameter :: B3 = -1.881058e-5_RealK

integer, parameter :: n_hs208_wl = 4
integer, parameter :: n_hs298_wl = 7
real(RealK), parameter :: hs_wls(n_hs298_wl) = [ &
  160.0_RealK, 165.0_RealK, 170.0_RealK, 173.0_RealK, &
  240.0_RealK, 245.0_RealK, 250.0_RealK ]
real(RealK), parameter :: hs208_xsc(n_hs208_wl) = [ &
  4.30e-20_RealK, 5.61e-20_RealK, 8.30e-20_RealK, 1.043e-19_RealK ]
real(RealK), parameter :: hs298_xsc(n_hs298_wl) = [ &
  4.40e-20_RealK, 5.77e-20_RealK, 8.71e-20_RealK, 1.10e-19_RealK, &
  9.732e-24_RealK, 4.25e-24_RealK, 1.36e-24_RealK ]

integer, parameter :: n_t = 4
real(RealK), parameter :: temperature(n_t) = [ &
  208.0_RealK, 238.0_RealK, 268.0_RealK, 298.0_RealK ]
character(*), parameter :: filename(n_t) = [ &
  'n2o_208K.dat', 'n2o_238K.dat', 'n2o_268K.dat', 'n2o_298K.dat' ]

integer, parameter :: n_wl = 68
real(RealK) :: wl(n_wl), hs_wl
real(RealK) :: xsc

integer :: lower
real(RealK) :: wt_lower

! Setup wavelengths from 173nm - 240nm at 1nm resolution
wl(1) = 173.0_RealK
do i=2, n_wl
  wl(i) = wl(i-1) + 1.0_RealK
end do

do j=1, n_t
  ! Open a new .dat file for each temperature
  open(newunit=io, file=filename(j), status='new', action='write', &
    iostat=stat, iomsg=msg)
  if (stat /= 0) then
    stop trim(msg)
  end if

  if (filename(j) == 'n2o_208K.dat' .or. &
      filename(j) == 'n2o_298K.dat') then
    ! For the 208K or 298K files include wavelengths from 160nm - 172nm
    ! using an interpolation from Hubrich and Stuhl (1980)
    hs_wl = hs_wls(1)
    do
      if (hs_wl > 172.9_RealK) exit
      lower = minloc(hs_wl - hs_wls, 1, hs_wl - hs_wls >= 0.0_RealK)
      wt_lower = 1.0_RealK - ( hs_wl - hs_wls(lower) ) &
               / ( hs_wls(lower+1) - hs_wls(lower) )
      if (filename(j) == 'n2o_208K.dat') then
        xsc = wt_lower * hs208_xsc(lower) &
            + (1.0_RealK - wt_lower) * hs208_xsc(lower+1)
      else
        xsc = wt_lower * hs298_xsc(lower) &
            + (1.0_RealK - wt_lower) * hs298_xsc(lower+1)
      end if
      write(io, '(0pf7.2, 2x, 1pe10.3)') hs_wl, xsc
      hs_wl = hs_wl + 1.0_RealK
    end do
  end if

  ! For all temperatures include wavelengths from 173nm - 240nm using
  ! Table 4C-3-3 in the JPL 19-5 report: Recommended Expression for the
  ! Absorption Cross Sections of N2O as a Function of Temperature
  do i=1, n_wl
    xsc = exp( &
          ( A0 + A1*wl(i) + A2*wl(i)**2 + A3*wl(i)**3 + A4*wl(i)**4 ) &
        + ( temperature(j) - 300.0_RealK ) * &
          exp( B0 + B1*wl(i) + B2*wl(i)**2 + B3*wl(i)**3 ) &
        )
    write(io, '(0pf7.2, 2x, 1pe10.3)') wl(i), xsc
  end do

  if (filename(j) == 'n2o_298K.dat') then
    ! For the 298K file include wavelengths from 241nm - 250nm
    ! using an interpolation from Hubrich and Stuhl (1980)
    hs_wl = hs_wls(5) + 1.0_RealK
    do
      if (hs_wl > 250.1_RealK) exit
      lower = minloc(hs_wl - hs_wls, 1, hs_wl - hs_wls >= 0.0_RealK)
      wt_lower = 1.0_RealK - ( hs_wl - hs_wls(lower) ) &
               / ( hs_wls(lower+1) - hs_wls(lower) )
      xsc = wt_lower * hs298_xsc(lower) &
          + (1.0_RealK - wt_lower) * hs298_xsc(lower+1)
      write(io, '(0pf7.2, 2x, 1pe10.3)') hs_wl, xsc
      hs_wl = hs_wl + 1.0_RealK
    end do
  end if
  close(io)
end do

end program calc_n2o_xsc
