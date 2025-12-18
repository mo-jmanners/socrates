program create_pt

implicit none

integer, parameter :: RealK=SELECTED_REAL_KIND(15, 307)
integer :: i, j
real(RealK) :: lp, linc(68), tinc
real(RealK) :: pressure(68), temperature(7, 68)
real(RealK), parameter :: uss_t(68) = [ &
394.322_RealK, &
318.684_RealK, &
254.477_RealK, &
220.871_RealK, &
203.618_RealK, &
194.454_RealK, &
189.995_RealK, &
187.773_RealK, &
187.122_RealK, &
187.472_RealK, &
188.339_RealK, &
190.432_RealK, &
194.779_RealK, &
199.151_RealK, &
202.563_RealK, &
205.976_RealK, &
209.601_RealK, &
213.750_RealK, &
217.898_RealK, &
222.759_RealK, &
228.114_RealK, &
233.479_RealK, &
238.001_RealK, &
242.522_RealK, &
247.045_RealK, &
251.838_RealK, &
256.631_RealK, &
261.271_RealK, &
264.890_RealK, &
268.510_RealK, &
270.671_RealK, &
270.413_RealK, &
266.493_RealK, &
262.500_RealK, &
258.405_RealK, &
254.387_RealK, &
250.399_RealK, &
246.493_RealK, &
242.598_RealK, &
238.824_RealK, &
235.057_RealK, &
231.302_RealK, &
228.908_RealK, &
227.236_RealK, &
225.786_RealK, &
224.509_RealK, &
223.272_RealK, &
222.060_RealK, &
220.825_RealK, &
219.579_RealK, &
218.340_RealK, &
217.159_RealK, &
216.700_RealK, &
216.700_RealK, &
216.700_RealK, &
216.700_RealK, &
216.700_RealK, &
216.700_RealK, &
216.705_RealK, &
218.608_RealK, &
226.712_RealK, &
235.108_RealK, &
243.864_RealK, &
252.959_RealK, &
262.362_RealK, &
272.143_RealK, &
282.250_RealK, &
288.200_RealK  ]

linc(1:8) = 1.0_RealK/4.0_RealK
linc(9:14) = 1.0_RealK/6.0_RealK
linc(15:22) = 1.0_RealK/8.0_RealK
linc(23:32) = 1.0_RealK/10.0_RealK
linc(33:68) = 1.0_RealK/12.0_RealK

lp=log10(1.1e-3_RealK)
do i=1, 68
  lp=lp+linc(i)
  pressure(i) = 10.0_RealK**(lp)
  if (i == 1) tinc = 50.0_RealK
  if (i == 2) tinc = 40.0_RealK
  if (i == 3) tinc = 30.0_RealK
  if (i > 3) tinc = 20.0_RealK
  do j=1, 7
    temperature(j, i) = uss_t(i) + (real(j, RealK)*tinc) - (tinc*4.0_RealK)
  end do
  write(*,'(1pe12.6, 7(0pf9.3))') pressure(i), temperature(:, i)
end do

end program create_pt
