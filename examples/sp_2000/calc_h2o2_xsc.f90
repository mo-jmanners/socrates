program calc_h2o2_xsc

! JPL 19-5 recommended cross-section for H2O2

! 190nm - 260nm at 298K: Table 4B-3-3 from JPL 19-5

! 260nm - 350nm, 200K - 320K: Table 4B-3-2 from JPL 19-5
! Nicovich & Wine (1988): Temperature-dependent absorption cross sections
! for hydrogen peroxide vapor. J. Geophys. Res. 1988, 93, 2417-2421,
! doi:10.1029/JD093iD03p02417.

implicit none

integer, parameter :: RealK=SELECTED_REAL_KIND(15, 307)
integer :: i, j, io, stat
logical :: exists
character(len=512) :: msg

real(RealK), parameter :: A0 = 6.4761e4_RealK
real(RealK), parameter :: A1 = -9.2170972e2_RealK
real(RealK), parameter :: A2 = 4.535649_RealK
real(RealK), parameter :: A3 = -4.4589016e-3_RealK
real(RealK), parameter :: A4 = -4.035101e-5_RealK
real(RealK), parameter :: A5 = 1.6878206e-7_RealK
real(RealK), parameter :: A6 = -2.652014e-10_RealK
real(RealK), parameter :: A7 = 1.5534675e-13_RealK
real(RealK), parameter :: B0 = 6.8123e3_RealK
real(RealK), parameter :: B1 = -5.1351e1_RealK
real(RealK), parameter :: B2 = 1.1522e-1_RealK
real(RealK), parameter :: B3 = -3.0493e-5_RealK
real(RealK), parameter :: B4 = -1.0924e-7_RealK

integer, parameter :: n_jpl_wl = 15
real(RealK), parameter :: jpl_wl(n_jpl_wl) = [ &
  190.0_RealK, 195.0_RealK, 200.0_RealK, 205.0_RealK, 210.0_RealK, &
  215.0_RealK, 220.0_RealK, 225.0_RealK, 230.0_RealK, 235.0_RealK, &
  240.0_RealK, 245.0_RealK, 250.0_RealK, 255.0_RealK, 260.0_RealK ]
real(RealK), parameter :: jpl_xsc(n_jpl_wl) = 1.0E-20_RealK * [ &
  67.2_RealK, 56.4_RealK, 47.5_RealK, 40.8_RealK, 35.7_RealK, &
  30.7_RealK, 25.8_RealK, 21.7_RealK, 18.2_RealK, 15.0_RealK, &
  12.4_RealK, 10.2_RealK, 8.3_RealK, 6.7_RealK, 5.315_RealK ]

integer, parameter :: n_t = 7
real(RealK), parameter :: temperature(n_t) = [ &
  200.0_RealK, 220.0_RealK, 240.0_RealK, &
  260.0_RealK, 280.0_Realk, 298.0_RealK, 320.0_RealK ]
character(*), parameter :: filename(n_t) = [ &
  'h2o2_200K.dat', 'h2o2_220K.dat', 'h2o2_240K.dat', &
  'h2o2_260K.dat', 'h2o2_280K.dat', 'h2o2_298K.dat', 'h2o2_320K.dat' ]

integer, parameter :: n_wl = 91
real(RealK) :: wl(n_wl), wl_298
real(RealK) :: xsc
real(RealK) :: chi

integer :: lower
real(RealK) :: wt_lower

! Setup wavelengths from 260nm - 350nm at 1nm resolution
wl(1) = 260.0_RealK
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

  if (filename(j) == 'h2o2_298K.dat') then
    ! For the 298K file include wavelengths from 190nm - 260nm using
    ! an interpolation from Table 4B-3-3 in the JPL 19-5 report
    wl_298 = jpl_wl(1)
    do
      if (wl_298 > 259.9_RealK) exit
      lower = minloc(wl_298 - jpl_wl, 1, wl_298 - jpl_wl >= 0.0_RealK)
      wt_lower = 1.0_RealK - ( wl_298 - jpl_wl(lower) ) &
               / ( jpl_wl(lower+1) - jpl_wl(lower) )
      xsc = wt_lower * jpl_xsc(lower) &
          + (1.0_RealK - wt_lower) * jpl_xsc(lower+1)
      write(io, '(0pf7.2, 2x, 1pe10.3)') wl_298, xsc
      wl_298 = wl_298 + 1.0_RealK
    end do
  end if

  ! For all other temperatures include wavelengths from 260nm - 350nm using
  ! Table 4B-3-2 in the JPL 19-5 report: Mathematical Expression for
  ! Absorption Cross Sections of H2O2 as a Function of Temperature
  do i=1, n_wl
    chi = 1.0_RealK / (1.0_RealK + exp( -1265.0_RealK / temperature(j) ))
    xsc = 1.0e-21_RealK * ( &
          chi * &
          ( A0 + A1*wl(i) + A2*wl(i)**2 + A3*wl(i)**3 + A4*wl(i)**4 &
                          + A5*wl(i)**5 + A6*wl(i)**6 + A7*wl(i)**7 ) &
        + ( 1.0_RealK - chi ) * &
          ( B0 + B1*wl(i) + B2*wl(i)**2 + B3*wl(i)**3 + B4*wl(i)**4 ) &
        )
    write(io, '(0pf7.2, 2x, 1pe10.3)') wl(i), xsc
  end do
  close(io)
end do

end program calc_h2o2_xsc
