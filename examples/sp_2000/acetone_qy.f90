! PT dependent QY for Acetone (CH3COCH3) from Blitz et al 2004,
! DOI: 10.1029/2003GL018793, relaxing to IUPAC values at 266nm.
program acetone_qy

implicit none

integer, parameter :: n_wl=63
integer, parameter :: n_t=4

real, parameter :: r_gas_dry = 287.026
real, parameter :: n_avogadro = 6.022045e+23
real, parameter :: mol_weight_air  = 28.966e-03
real, parameter :: per_m3_to_cm3 = 1.0e-6

integer :: i, j

real :: t(n_t)=[218.0, 248.0, 273.0, 295.0] ! Temperature (K)
real :: p(n_t)=[15400.0, 27380.0, 48700.0, 86600.0] ! Pressure (Pa)
real :: wl(n_wl), wl_param ! Wavelength (nm)
real :: M(n_t) ! Density (molecule cm-3)

real :: qy_co(n_t, n_wl)
real :: qy_ch3co(n_t, n_wl)
real :: qy_total(n_t, n_wl)
real :: a0, b0, A_0
real :: a1, b1, A_1
real :: a2, b2, A_2
real :: a3, b3, c3, A_3
real :: a4, b4, A_4

! Set wavelengths
do i=1, n_wl-1
  wl(i) = 265.0+real(i)
end do
wl(n_wl) = 327.5

! Calculate densities (molecule cm-3) assuming tropospheric P-T values
do j=1, n_t
  M(j) = ( p(j)/(r_gas_dry*t(j)) )*n_avogadro*per_m3_to_cm3/mol_weight_air
end do

do i=1, n_wl
  do j=1, n_t
    if (wl(i) < 279.0) then
      wl_param = 279.0
    else
      wl_param = wl(i)
    end if

    a0 = 0.35*(t(j)/295.0)**(-1.28)
    b0 = 0.068*(t(j)/295.0)**(-2.65)
    A_0 = (a0/(1.0-a0))*exp(b0*(wl(i)-248.0))
    qy_co(j, i) = 1.0/(1.0+A_0)

    if (wl(i) < 302.0) then
      a1 = 1.6e-19*(t(j)/295.0)**(-2.38)
      b1 = 0.55e-3*(t(j)/295.0)**(-3.19)
      A_1 = a1*exp(-b1*((1.0e7/wl(i))-33113.0))
      qy_ch3co(j, i) = (1.0-qy_co(j, i))/(1.0+A_1*M(j))
    else
      a2 = 1.62e-17*(t(j)/295.0)**(-10.03)
      b2 = 1.79e-3*(t(j)/295.0)**(-1.364)
      A_2 = a2*exp(-b2*((1.0e7/wl(i))-30488.0))
    
      a3 = 26.29*(t(j)/295.0)**(-6.59)
      b3 = 5.72e-7*(t(j)/295.0)**(-2.93)
      c3 = 30006.0*(t(j)/295.0)**(-0.064)
      A_3 = a3*exp(-b3*((1.0e7/wl(i))-c3)**2)
    
      a4 = 1.67e-15*(t(j)/295.0)**(-7.25)
      b4 = 2.08e-3*(t(j)/295.0)**(-1.16)
      A_4 = a4*exp(-b4*((1.0e7/wl(i))-30488.0))
    
      qy_ch3co(j, i) = (1.0-qy_co(j, i))*(1.0+A_4*M(j)+A_3) &
                     / ((1.0+A_2*M(j)+A_3)*(1.0+A_4*M(j)))
    end if

    if (wl(i) < 266.0) then
      qy_ch3co(j, i) = 0.55
      qy_co(j, i) = 0.45
    else if (wl(i) < 279.0) then
      qy_ch3co(j, i) = (0.55*(279.0-wl(i)) + qy_ch3co(j, i)*(wl(i)-266.0)) &
                     / (279.0 - 266.0)
      qy_co(j, i) = (0.45*(279.0-wl(i)) + qy_co(j, i)*(wl(i)-266.0)) &
                  / (279.0 - 266.0)
    end if
  end do
end do

write(*,'(a)') 'Densities (molecule cm-3):'
write(*, *) M
write(*, *)
write(*, *) '---- File: acetone_1.qy ----'
write(*,'(a)') '*TEMPERATURES 4'
write(*, '(3x,4f12.1)') t
do i=1, n_wl
  write(*, '(f5.1, 4f12.7)') wl(i), qy_ch3co(:, i)
end do
write(*, *)
write(*, *) '---- File: acetone_2.qy ----'
write(*,'(a)') '*TEMPERATURES 4'
write(*, '(3x,4f12.1)') t
do i=1, n_wl
  write(*, '(f5.1, 4f12.7)') wl(i), qy_co(:, i)
end do
write(*, *)
write(*, *) '---- Total QY ----'
write(*,'(a)') '*TEMPERATURES 4'
write(*, '(3x,4f12.1)') t
do i=1, n_wl
  write(*, '(f5.1, 4f12.7)') wl(i), qy_co(:, i)+qy_ch3co(:, i)
end do

end program acetone_qy
