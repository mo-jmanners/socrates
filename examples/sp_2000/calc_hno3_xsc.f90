program calc_hno3_xsc

! JPL 19-5 recommended cross-section for HNO3

! 196nm - 350nm, 200K - 360K: Table 4C-8 from JPL 19-5
! Burkholder et al (1993): Temperature dependence of the HNO3 UV absorption
! cross sections, J. Geophys. Res. 98(D12), 22937-22948
! doi: 10.1029/93JD02178

implicit none

integer, parameter :: RealK=SELECTED_REAL_KIND(15, 307)
integer :: i, j, io, stat
logical :: exists
character(len=512) :: msg

integer, parameter :: n_jpl_wl = 78
real(RealK), parameter :: jpl_b(n_jpl_wl) = [ &
 1.70e-3_RealK, &
 1.65e-3_RealK, &
 1.66e-3_RealK, &
 1.69e-3_RealK, &
 1.74e-3_RealK, &
 1.77e-3_RealK, &
 1.85e-3_RealK, &
 1.97e-3_RealK, &
 2.08e-3_RealK, &
 2.17e-3_RealK, &
 2.17e-3_RealK, &
 2.21e-3_RealK, &
 2.15e-3_RealK, &
 2.06e-3_RealK, &
 1.96e-3_RealK, &
 1.84e-3_RealK, &
 1.78e-3_RealK, &
 1.80e-3_RealK, &
 1.86e-3_RealK, &
 1.90e-3_RealK, &
 1.97e-3_RealK, &
 1.97e-3_RealK, &
 1.97e-3_RealK, &
 1.88e-3_RealK, &
 1.75e-3_RealK, &
 1.61e-3_RealK, &
 1.44e-3_RealK, &
 1.34e-3_RealK, &
 1.23e-3_RealK, &
 1.18e-3_RealK, &
 1.14e-3_RealK, &
 1.12e-3_RealK, &
 1.14e-3_RealK, &
 1.14e-3_RealK, &
 1.18e-3_RealK, &
 1.22e-3_RealK, &
 1.25e-3_RealK, &
 1.45e-3_RealK, &
 1.49e-3_RealK, &
 1.56e-3_RealK, &
 1.64e-3_RealK, &
 1.69e-3_RealK, &
 1.78e-3_RealK, &
 1.87e-3_RealK, &
 1.94e-3_RealK, &
 2.04e-3_RealK, &
 2.15e-3_RealK, &
 2.27e-3_RealK, &
 2.38e-3_RealK, &
 2.52e-3_RealK, &
 2.70e-3_RealK, &
 2.92e-3_RealK, &
 3.10e-3_RealK, &
 3.24e-3_RealK, &
 3.52e-3_RealK, &
 3.77e-3_RealK, &
 3.91e-3_RealK, &
 4.23e-3_RealK, &
 4.70e-3_RealK, &
 5.15e-3_RealK, &
 5.25e-3_RealK, &
 5.74e-3_RealK, &
 6.45e-3_RealK, &
 6.70e-3_RealK, &
 7.16e-3_RealK, &
 7.55e-3_RealK, &
 8.16e-3_RealK, &
 9.75e-3_RealK, &
 9.93e-3_RealK, &
 9.60e-3_RealK, &
10.50e-3_RealK, &
10.80e-3_RealK, &
11.80e-3_RealK, &
11.80e-3_RealK, &
 9.30e-3_RealK, &
12.10e-3_RealK, &
11.90e-3_RealK, &
 9.30e-3_RealK ]

real(RealK), parameter :: xsc298(n_jpl_wl) = [ &
9.40e-18_RealK, &
7.70e-18_RealK, &
5.88e-18_RealK, &
4.47e-18_RealK, &
3.28e-18_RealK, &
2.31e-18_RealK, &
1.56e-18_RealK, &
1.04e-18_RealK, &
6.75e-19_RealK, &
4.39e-19_RealK, &
2.92e-19_RealK, &
2.00e-19_RealK, &
1.49e-19_RealK, &
1.18e-19_RealK, &
9.61e-20_RealK, &
8.02e-20_RealK, &
6.82e-20_RealK, &
5.75e-20_RealK, &
4.87e-20_RealK, &
4.14e-20_RealK, &
3.36e-20_RealK, &
2.93e-20_RealK, &
2.58e-20_RealK, &
2.34e-20_RealK, &
2.16e-20_RealK, &
2.06e-20_RealK, &
2.00e-20_RealK, &
1.97e-20_RealK, &
1.96e-20_RealK, &
1.95e-20_RealK, &
1.95e-20_RealK, &
1.93e-20_RealK, &
1.91e-20_RealK, &
1.87e-20_RealK, &
1.83e-20_RealK, &
1.77e-20_RealK, &
1.70e-20_RealK, &
1.62e-20_RealK, &
1.53e-20_RealK, &
1.44e-20_RealK, &
1.33e-20_RealK, &
1.23e-20_RealK, &
1.12e-20_RealK, &
1.01e-20_RealK, &
9.09e-21_RealK, &
8.07e-21_RealK, &
7.09e-21_RealK, &
6.15e-21_RealK, &
5.32e-21_RealK, &
4.53e-21_RealK, &
3.81e-21_RealK, &
3.16e-21_RealK, &
2.63e-21_RealK, &
2.08e-21_RealK, &
1.67e-21_RealK, &
1.33e-21_RealK, &
1.05e-21_RealK, &
8.14e-22_RealK, &
6.28e-22_RealK, &
4.68e-22_RealK, &
3.62e-22_RealK, &
2.71e-22_RealK, &
1.97e-22_RealK, &
1.54e-22_RealK, &
1.08e-22_RealK, &
8.20e-23_RealK, &
6.13e-23_RealK, &
4.31e-23_RealK, &
3.19e-23_RealK, &
2.43e-23_RealK, &
1.96e-23_RealK, &
1.42e-23_RealK, &
1.03e-23_RealK, &
8.61e-24_RealK, &
6.94e-24_RealK, &
5.01e-24_RealK, &
4.15e-24_RealK, &
4.17e-24_RealK ]

integer, parameter :: n_t = 5
real(RealK), parameter :: temperature(n_t) = [ &
  200.0_RealK, 220.0_RealK, 240.0_RealK, &
  260.0_RealK, 280.0_RealK ]
character(*), parameter :: filename(n_t) = [ &
  'hno3_200K.dat', 'hno3_220K.dat', 'hno3_240K.dat', &
  'hno3_260K.dat', 'hno3_280K.dat' ]

real(RealK) :: wl(n_jpl_wl)
real(RealK) :: xsc


! Setup wavelengths from 196nm - 350nm at 2nm resolution
wl(1) = 196.0_RealK
do i=2, n_jpl_wl
  wl(i) = wl(i-1) + 2.0_RealK
end do

do j=1, n_t
  ! Open a new .dat file for each temperature
  open(newunit=io, file=filename(j), status='new', action='write', &
    iostat=stat, iomsg=msg)
  if (stat /= 0) then
    stop trim(msg)
  end if

  do i=1, n_jpl_wl
    xsc = xsc298(i) * exp( jpl_b(i) * (temperature(j)-298.0_RealK) )
    write(io, '(0pf7.2, 2x, 1pe10.3)') wl(i), xsc
  end do
  close(io)
end do

end program calc_hno3_xsc
