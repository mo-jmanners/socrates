program calc_n2o5_xsc

! JPL 19-5 recommended cross-section for N2O5: 200-420nm

! 152-198nm: Osborne et al (2000), "Vacuum ultraviolet spectrum of
! dinitrogen pentoxide", J. Quant. Spectrosc. Radiat. Transfer 64, 67-74
! DOI: 10.1016/S0022-4073(99)00104-1

implicit none

integer, parameter :: RealK=SELECTED_REAL_KIND(15, 307)
integer :: i, j, io, stat
logical :: exists
character(len=512) :: msg

integer, parameter :: n_jpl_wl = 127
real(RealK), parameter :: jpl_wls(n_jpl_wl) = [ &
152.0_RealK, &
154.0_RealK, &
156.0_RealK, &
158.0_RealK, &
160.0_RealK, &
162.0_RealK, &
164.0_RealK, &
166.0_RealK, &
168.0_RealK, &
170.0_RealK, &
172.0_RealK, &
174.0_RealK, &
176.0_RealK, &
178.0_RealK, &
180.0_RealK, &
182.0_RealK, &
184.0_RealK, &
186.0_RealK, &
188.0_RealK, &
190.0_RealK, &
192.0_RealK, &
194.0_RealK, &
196.0_RealK, &
198.0_RealK, &
200.0_RealK, &
202.0_RealK, &
204.0_RealK, &
206.0_RealK, &
208.0_RealK, &
210.0_RealK, &
212.0_RealK, &
214.0_RealK, &
216.0_RealK, &
218.0_RealK, &
220.0_RealK, &
222.0_RealK, &
224.0_RealK, &
226.0_RealK, &
228.0_RealK, &
230.0_RealK, &
232.0_RealK, &
234.0_RealK, &
236.0_RealK, &
238.0_RealK, &
240.0_RealK, &
242.0_RealK, &
244.0_RealK, &
246.0_RealK, &
248.0_RealK, &
250.0_RealK, &
252.0_RealK, &
254.0_RealK, &
256.0_RealK, &
258.0_RealK, &
260.0_RealK, &
262.0_RealK, &
264.0_RealK, &
266.0_RealK, &
268.0_RealK, &
270.0_RealK, &
272.0_RealK, &
274.0_RealK, &
276.0_RealK, &
278.0_RealK, &
280.0_RealK, &
282.0_RealK, &
284.0_RealK, &
286.0_RealK, &
288.0_RealK, &
290.0_RealK, &
292.0_RealK, &
294.0_RealK, &
296.0_RealK, &
298.0_RealK, &
300.0_RealK, &
302.0_RealK, &
304.0_RealK, &
306.0_RealK, &
308.0_RealK, &
310.0_RealK, &
312.0_RealK, &
314.0_RealK, &
316.0_RealK, &
318.0_RealK, &
320.0_RealK, &
322.0_RealK, &
324.0_RealK, &
326.0_RealK, &
328.0_RealK, &
330.0_RealK, &
332.0_RealK, &
334.0_RealK, &
336.0_RealK, &
338.0_RealK, &
340.0_RealK, &
342.0_RealK, &
344.0_RealK, &
346.0_RealK, &
348.0_RealK, &
350.0_RealK, &
352.0_RealK, &
354.0_RealK, &
356.0_RealK, &
358.0_RealK, &
360.0_RealK, &
362.0_RealK, &
364.0_RealK, &
366.0_RealK, &
368.0_RealK, &
370.0_RealK, &
372.0_RealK, &
374.0_RealK, &
376.0_RealK, &
378.0_RealK, &
380.0_RealK, &
382.0_RealK, &
384.0_RealK, &
386.0_RealK, &
388.0_RealK, &
390.0_RealK, &
392.0_RealK, &
394.0_RealK, &
396.0_RealK, &
398.0_RealK, &
400.0_RealK, &
410.0_RealK, &
420.0_RealK  ]
real(RealK), parameter :: jpl_xsc(n_jpl_wl) = [ &
3.38e-17_RealK, &
3.65e-17_RealK, &
3.82e-17_RealK, &
3.98e-17_RealK, &
4.03e-17_RealK, &
4.00e-17_RealK, &
3.91e-17_RealK, &
3.72e-17_RealK, &
3.53e-17_RealK, &
3.32e-17_RealK, &
3.10e-17_RealK, &
2.88e-17_RealK, &
2.65e-17_RealK, &
2.44e-17_RealK, &
2.25e-17_RealK, &
2.07e-17_RealK, &
1.92e-17_RealK, &
1.77e-17_RealK, &
1.64e-17_RealK, &
1.49e-17_RealK, &
1.38e-17_RealK, &
1.25e-17_RealK, &
1.13e-17_RealK, &
1.02e-17_RealK, &
9.10e-18_RealK, &
8.42e-18_RealK, &
7.71e-18_RealK, &
6.82e-18_RealK, &
5.85e-18_RealK, &
4.45e-18_RealK, &
3.81e-18_RealK, &
3.22e-18_RealK, &
2.67e-18_RealK, &
2.20e-18_RealK, &
1.81e-18_RealK, &
1.51e-18_RealK, &
1.29e-18_RealK, &
1.13e-18_RealK, &
9.84e-19_RealK, &
8.82e-19_RealK, &
8.05e-19_RealK, &
7.40e-19_RealK, &
6.92e-19_RealK, &
6.46e-19_RealK, &
5.98e-19_RealK, &
5.31e-19_RealK, &
4.93e-19_RealK, &
4.56e-19_RealK, &
4.19e-19_RealK, &
3.86e-19_RealK, &
3.55e-19_RealK, &
3.26e-19_RealK, &
2.99e-19_RealK, &
2.75e-19_RealK, &
2.52e-19_RealK, &
2.31e-19_RealK, &
2.11e-19_RealK, &
1.94e-19_RealK, &
1.78e-19_RealK, &
1.62e-19_RealK, &
1.49e-19_RealK, &
1.37e-19_RealK, &
1.24e-19_RealK, &
1.14e-19_RealK, &
1.05e-19_RealK, &
9.59e-20_RealK, &
8.74e-20_RealK, &
7.94e-20_RealK, &
7.20e-20_RealK, &
6.52e-20_RealK, &
5.88e-20_RealK, &
5.29e-20_RealK, &
4.75e-20_RealK, &
4.26e-20_RealK, &
3.81e-20_RealK, &
3.40e-20_RealK, &
3.03e-20_RealK, &
2.70e-20_RealK, &
2.40e-20_RealK, &
2.13e-20_RealK, &
1.90e-20_RealK, &
1.68e-20_RealK, &
1.49e-20_RealK, &
1.33e-20_RealK, &
1.18e-20_RealK, &
1.05e-20_RealK, &
9.30e-21_RealK, &
8.26e-21_RealK, &
7.35e-21_RealK, &
6.54e-21_RealK, &
5.82e-21_RealK, &
5.18e-21_RealK, &
4.62e-21_RealK, &
4.12e-21_RealK, &
3.68e-21_RealK, &
3.28e-21_RealK, &
2.93e-21_RealK, &
2.62e-21_RealK, &
2.34e-21_RealK, &
2.10e-21_RealK, &
1.88e-21_RealK, &
1.67e-21_RealK, &
1.49e-21_RealK, &
1.33e-21_RealK, &
1.20e-21_RealK, &
1.07e-21_RealK, &
9.58e-22_RealK, &
8.52e-22_RealK, &
7.63e-22_RealK, &
6.85e-22_RealK, &
6.13e-22_RealK, &
5.45e-22_RealK, &
4.84e-22_RealK, &
4.31e-22_RealK, &
3.83e-22_RealK, &
3.41e-22_RealK, &
3.05e-22_RealK, &
2.73e-22_RealK, &
2.42e-22_RealK, &
2.15e-22_RealK, &
1.93e-22_RealK, &
1.72e-22_RealK, &
1.50e-22_RealK, &
1.34e-22_RealK, &
1.4e-22_RealK, &
9.0e-23_RealK, &
5.0e-23_RealK  ]

integer, parameter :: n_jpl_t_wl = 17
real(RealK), parameter :: jpl_t_wls(n_jpl_t_wl) = [ &
260.0_RealK, &
270.0_RealK, &
280.0_RealK, &
290.0_RealK, &
300.0_RealK, &
310.0_RealK, &
320.0_RealK, &
330.0_RealK, &
340.0_RealK, &
350.0_RealK, &
360.0_RealK, &
370.0_RealK, &
380.0_RealK, &
390.0_RealK, &
400.0_RealK, &
410.0_RealK, &
420.0_RealK  ]
real(RealK), parameter :: jpl_a(n_jpl_t_wl) = [ &
-18.27_RealK, &
-18.42_RealK, &
-18.59_RealK, &
-18.72_RealK, &
-18.84_RealK, &
-18.90_RealK, &
-18.93_RealK, &
-18.87_RealK, &
-18.77_RealK, &
-18.71_RealK, &
-18.31_RealK, &
-18.14_RealK, &
-18.01_RealK, &
-18.42_RealK, &
-18.59_RealK, &
-18.13_RealK, &
-18.37_RealK  ]
real(RealK), parameter :: jpl_b(n_jpl_t_wl) = [ &
-0.091_RealK, &
-0.104_RealK, &
-0.112_RealK, &
-0.135_RealK, &
-0.170_RealK, &
-0.226_RealK, &
-0.294_RealK, &
-0.388_RealK, &
-0.492_RealK, &
-0.583_RealK, &
-0.770_RealK, &
-0.885_RealK, &
-0.992_RealK, &
-0.949_RealK, &
-0.966_RealK, &
-1.160_RealK, &
-1.160_RealK  ]

integer, parameter :: n_t = 5
real(RealK), parameter :: temperature(n_t) = [ &
  233.0_RealK, 235.0_RealK, 255.0_RealK, 275.0_RealK, 295.0_RealK ]
character(*), parameter :: filename(n_t) = [ &
  'n2o5_215K.dat', 'n2o5_235K.dat', 'n2o5_255K.dat', &
  'n2o5_275K.dat', 'n2o5_295K.dat' ]

real(RealK) :: wl, xsc, a, b

integer :: lower
real(RealK) :: wt_lower


do j=1, n_t
  ! Open a new .dat file for each temperature
  open(newunit=io, file=filename(j), status='new', action='write', &
    iostat=stat, iomsg=msg)
  if (stat /= 0) then
    stop trim(msg)
  end if

  if (filename(j) == 'n2o5_295K.dat') then
    ! For the 295K file include wavelengths from 152nm - 259nm
    wl = jpl_wls(1)
    do
      if (wl > 259.9_RealK) exit
      lower = minloc(wl - jpl_wls, 1, wl - jpl_wls >= 0.0_RealK)
      wt_lower = 1.0_RealK - ( wl - jpl_wls(lower) ) &
               / ( jpl_wls(lower+1) - jpl_wls(lower) )
      xsc = wt_lower * jpl_xsc(lower) &
          + (1.0_RealK - wt_lower) * jpl_xsc(lower+1)
      write(io, '(0pf7.2, 2x, 1pe10.3)') wl, xsc
      wl = wl + 1.0_RealK
    end do
  end if

  ! For all temperatures include wavelengths from 260nm - 420nm using
  ! coefficients from Harwood et al. (1993), "Temperature-dependent absorption
  ! cross sections of N2O5", J. Photochem. Photobiol. A 73, 167-175
  ! DOI: 10.1016/1010-6030(93)90001-2,
  wl = jpl_t_wls(1)
  do
    if (wl > 420.1_RealK) exit
    lower = minloc(wl - jpl_t_wls, 1, wl - jpl_t_wls >= 0.0_RealK)
    wt_lower = 1.0_RealK - ( wl - jpl_t_wls(lower) ) &
             / ( jpl_t_wls(lower+1) - jpl_t_wls(lower) )
    b = wt_lower * jpl_b(lower) + (1.0_RealK - wt_lower) * jpl_b(lower+1)
    a = wt_lower * jpl_a(lower) + (1.0_RealK - wt_lower) * jpl_a(lower+1)
    xsc = 10.0_RealK**(a + 1000.0_RealK*b/temperature(j))
    write(io, '(0pf7.2, 2x, 1pe10.3)') wl, xsc
    wl = wl + 1.0_RealK
  end do
  close(io)
end do

end program calc_n2o5_xsc
