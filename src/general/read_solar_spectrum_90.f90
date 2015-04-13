! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read a solar spectrum.
!
SUBROUTINE read_solar_spectrum_90 &
!
(SolarSpec, ierr) 
!
! Description:
!   The solar spetrum is read a a file of irradiances against
!   wavelengths.
!
!
! Modules to set types of variables:
  USE realtype_rd
  USE def_std_io_icf
  USE def_solarspec
  USE error_pcf
!
!
  IMPLICIT NONE
!
!
! Include header files.
!
! Dummy arguments.
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
  TYPE (StrSolarSpec), Intent(OUT) :: SolarSpec
!   Solar spectrum
!
! Local variables.
  CHARACTER (LEN=80) :: line
!   Line of input data
  INTEGER :: iu_solar
!   Unit number for the solar spectrum
  INTEGER :: ios
!   I/O error flag
  INTEGER :: i
!   Loop variable
  LOGICAL :: l_count
!   Flag for counting points
!
! Subroutines called:
  EXTERNAL &
      open_file_in
!
!
!
! Obtain the file containing the solar spectrum.
  CALL get_free_unit(ierr, iu_solar)
  CALL open_file_in(ierr, iu_solar, & 
    'Enter the name of the file containing the solar irradiance data.')
  IF (ierr /= i_normal) RETURN
!
! Read first to find the number of points in the spectrum.
  SolarSpec%n_points = 0
  l_count=.FALSE.
  DO
    READ(iu_solar, '(A)', IOSTAT=ios) line
    IF (ios /= 0) THEN
      EXIT
    ELSE IF (line(1:11) == '*BEGIN_DATA') THEN
      l_count=.TRUE.
    ELSE IF (line(1:4) == '*END') THEN
      l_count=.FALSE.
    ELSE IF (l_count) THEN
      SolarSpec%n_points=SolarSpec%n_points+1
    ENDIF
  ENDDO
!
  ALLOCATE(SolarSpec%wavelength(SolarSpec%n_points))
  ALLOCATE(SolarSpec%irrad(SolarSpec%n_points))
!
! Read in the file.
  REWIND(iu_solar)
  DO
    READ(iu_solar, '(A)', IOSTAT=ios) line
    IF (line(1:11) == '*BEGIN_DATA') THEN
      DO i = 1, SolarSpec%n_points
        READ(iu_solar, *, IOSTAT=ios) &
          SolarSpec%wavelength(i), &
          SolarSpec%irrad(i)
        IF (ios /= 0) THEN
          WRITE(iu_err, '(/A)') '*** Error: Corrupt solar spectrum.'
          ierr=i_err_fatal
          RETURN
        ENDIF
      ENDDO
      EXIT
    ENDIF
  ENDDO
!
! Check that reading was carried out: failure to read the data
! will cause n_solar_points to be 0.
  IF (SolarSpec%n_points == 0) THEN
    WRITE(iu_err, '(/a)') &
      '*** Error: No data were read. Check format of file of irradiances.'
    ierr=i_err_fatal
    RETURN
  ENDIF
!
!
  RETURN
END
