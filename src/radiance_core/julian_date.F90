! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! @brief Set the Julian date for a Gregorian calendar date and time

module julian_date_mod
implicit none
contains

function julian_date(year, month, day, second)

use realtype_rd, only: RealK
implicit none

real(RealK) :: julian_date
integer, intent(in) :: year
integer, intent(in), optional :: month, day, second

real(Realk), parameter :: julian_day_one = 1721425.5_RealK
! Julian day number for the beginning of the first day (year 1) of
! the Gregorian calendar.

integer, parameter :: days_to_month(12) = &
   [0, 31, 59, 90,120,151,181,212,243,273,304,334]


julian_date = julian_day_one + real( 365*(year-1) &
  + (year-1)/4 - (year-1)/100 + (year-1)/400, RealK )

if (present(month)) then
  julian_date = julian_date + real( days_to_month(month), RealK )
  if (month > 2) then
    if (mod(year, 400) == 0 .or. &
       (mod(year, 4) == 0 .and. mod(year, 100) /= 0)) then
      ! Leap year
      julian_date = julian_date + 1.0_RealK
    end if
  end if
end if

if (present(day)) then
  julian_date = julian_date + real( day, RealK ) - 1.0_RealK
end if

if (present(second)) then
  julian_date = julian_date + real( second, RealK ) / 86400.0_RealK
end if

end function julian_date
end module julian_date_mod
