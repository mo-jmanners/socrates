! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the weighting function for correlated-k.
!
SUBROUTINE select_weight_ck_90 &
!
(i_weight, SolarSpec, abs_path, l_interactive, ierr)
!
! Description:
!   A list of possible weighting functions is displayed and
!   the user selects one. If solar weighting is used the
!   solar spectrum is read in.
!
!
!
! Modules used:
  USE realtype_rd
  USE def_solarspec
  USE weighting_pcf
  USE error_pcf
  USE def_std_io_icf
!
!
  IMPLICIT NONE
!
!
!
! Dummy variables.
  LOGICAL, Intent(IN) :: l_interactive
!   Flag for interactive operation
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
  INTEGER, Intent(OUT) :: i_weight
!   Method of weighting
  TYPE (StrSolarSpec), Intent(OUT) :: SolarSpec
!   Solar spectral irradiance data
  REAL (RealK), ALLOCATABLE, Intent(INOUT) :: abs_path(:)
!   Cumulative absorber pathlength to each pressure in lookup table
!
! Local variables
!
!
  INTEGER :: i
  INTEGER :: ios
!   I/O error flag
  INTEGER :: iu_path
!   Unit number for cumulative path file
  INTEGER :: n_p_path
!   Number of pressure levels read for cumulative path
  CHARACTER(LEN=256) :: file_in
!   Input file for cumulative absorber path
  CHARACTER (LEN=256) :: line
!   Line of input data
  LOGICAL :: l_count
!   Flag for counting points
  REAL (RealK) :: pressure
!   Dummy pressure variable
!
! Display menu of weightings.
  WRITE(iu_stdout, '(/a)') &
    'Select the method of weighting the transmittances.'
  WRITE(iu_stdout, '(3x, i2, a)' ) ip_weight_planck, &
    '. Planckian weighting at transmission temperature.'
  WRITE(iu_stdout, '(3x, i2, a)' ) ip_weight_d_planck, &
    '. Differential planckian weighting at transmission temperature.'
  WRITE(iu_stdout, '(3x, i2, a)' ) ip_weight_solar, &
    '. TOA solar spectral weighting.'
  WRITE(iu_stdout, '(3x, i2, a)' ) ip_weight_uniform, &
    '. Uniform weighting.'
  WRITE(iu_stdout, '(3x, i2, a)' ) ip_weight_solar_path, &
    '. Solar spectral weighting with absorption along path.'
!
  WRITE(iu_stdout, '(a/)') 'Enter required number.'
  DO
    READ(iu_stdin, *, iostat=ios) i_weight
    IF ( (ios /= 0)                         .OR.   &
         ( (i_weight /= ip_weight_planck)   .AND.  &
           (i_weight /= ip_weight_d_planck) .AND.  &
           (i_weight /= ip_weight_solar)    .AND.  &
           (i_weight /= ip_weight_solar_path).AND. &
           (i_weight /= ip_weight_uniform) ) ) then
      WRITE(iu_err, '(a)') '+++ Erroneous response:'
      IF (l_interactive) then
        WRITE(iu_stdout, '(a)') 'please re-enter.'
      ELSE
        ierr=i_err_fatal
        RETURN
      ENDIF
    ELSE
      EXIT
    ENDIF
  ENDDO

! If solar weighting is used read the data in.
  if (i_weight == ip_weight_solar .OR. i_weight == ip_weight_solar_path) then
    CALL read_solar_spectrum(SolarSpec, ierr)
    if (ierr /= i_normal) RETURN
  end if

  ! If solar path weighting is used read in the cumulative absorber path
  if (i_weight == ip_weight_solar_path) then
    write(*, '(/a)') &
      'Filename for cumulative absorber path for each P in the lookup table:'
    file_in=''
    read(iu_stdin, '(a)', iostat=ios) file_in
    open(newunit=iu_path, file=trim(file_in), &
      status='old', action="read", iostat=ios)
    if (ios == 0) then
      n_p_path = 0
      l_count=.false.
      do
        read(iu_path, '(a)', iostat=ios) line
        if (ios /= 0) then
          exit
        else if (line(1:11) == '*BEGIN_DATA') then
          l_count=.true.
        else if (line(1:16) == '*BEGIN_2COL_DATA') then
          l_count=.true.
        else if (line(1:4) == '*END') then
          l_count=.false.
        else if (l_count) then
          n_p_path = n_p_path + 1
        end if
      end do
      if (size(abs_path) /= n_p_path) then
        write(iu_err, '(a)') &
          'Error: wrong number of pressures in absorber path file'
        ierr=i_err_fatal
        return
      end if
      rewind(iu_path)
      do
        read(iu_path, '(a)', iostat=ios) line
        if (line(1:11) == '*BEGIN_DATA') then
          do i = 1, n_p_path
            read(iu_path, *, iostat=ios) abs_path(i)
          end do
          exit
        else if (line(1:16) == '*BEGIN_2COL_DATA') then
          do i = 1, n_p_path
            read(iu_path, *, iostat=ios) pressure, abs_path(i)
          end do
          exit
        end if
      end do
    else
      write(iu_err, '(a)') 'Error: problem reading absorber path file'
      ierr=i_err_fatal
      return
    end if
    close(iu_path)
  end if
  
END SUBROUTINE select_weight_ck_90
