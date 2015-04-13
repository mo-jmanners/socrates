! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to make spectral blocks of type 1.
!
SUBROUTINE make_block_1 &
!
(Spectrum, l_interactive, ierr)
!
! Description:
!   This routine defines the wavelengths of the spectral bands
!
! Method:
!   Straightforward.
!
!- End of header
!
!
!
! Modules to set types of variables:
  USE realtype_rd
  USE def_spectrum
  USE def_std_io_icf
  USE error_pcf
!
!
  IMPLICIT none
!
!
!
! Dummy arguments.
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
  LOGICAL, Intent(IN) :: l_interactive
!   Flag for interactive operation
  TYPE(StrSpecData), Intent(INOUT), Target :: Spectrum
!
! Local variables.
!
!
  CHARACTER  (LEN=1) :: char_unit*1
!   Unit type
  INTEGER :: ios
!   Reading error flag
  INTEGER :: i
!   Loop variable
  LOGICAL :: l_metre = .FALSE.
!   Logical for units of metres
  LOGICAL :: l_inverse_cm = .FALSE.
!   Logical for units of inverse cm
  LOGICAL :: l_micron = .FALSE.
!   Logical for units of microns
  REAL  (RealK) :: wavelength_temp
!   Temporary for interchange
  TYPE (StrSpecBasic), Pointer :: SpBasic
!   Pointer to basic components of the spectrum
  TYPE (StrSpecDim), Pointer :: SpDim
!   Pointer to dimensions within the spectrum
!
!
!
  SpDim   => Spectrum%Dim
  ALLOCATE(Spectrum%Basic%wavelength_long(SpDim%nd_band))
  ALLOCATE(Spectrum%Basic%wavelength_short(SpDim%nd_band))
  SpBasic => Spectrum%Basic
!
! Obtain the limits for each band.
  WRITE(iu_stdout, '(/A, /A/)') &
    'For each band in turn specify the limits of the bands in ', &
    'metres, inverse centimetres or microns.'
  WRITE(iu_stdout, '(A)') 'type "m" for metres, "c" for inverse ' &
    //'centimetres, or "u" for microns.'
  DO
    READ(iu_stdin, '(A)') char_unit
    IF ( (char_unit.eq.'m').OR.(char_unit.eq.'M') ) THEN
      l_metre=.TRUE.
      EXIT
    ELSE IF ( (char_unit.eq.'c').OR.(char_unit.eq.'C') ) THEN
      l_inverse_cm=.TRUE.
      EXIT
    ELSE IF ( (char_unit.eq.'u').OR.(char_unit.eq.'U') ) THEN
      l_micron=.TRUE.
      EXIT
    ELSE
      WRITE(iu_err, '(a)') 'Unknown reponse:'
      IF (l_interactive) THEN
        WRITE(iu_stdout, '(a)') 'Please re-enter.'
      ELSE
        ierr=i_err_fatal
        RETURN
      ENDIF
    ENDIF
  ENDDO
!
  DO i=1, SpBasic%n_band
    WRITE(iu_stdout, '(A21, 1X, I5)') 'Enter limits for band ', i
    DO
      READ(iu_stdin, *, IOSTAT=ios) &
        SpBasic%wavelength_short(i), SpBasic%wavelength_long(i)
      IF (ios == 0) THEN
        EXIT
      ELSE IF (l_interactive) THEN
        WRITE(iu_stdout, '(A)') 'Please re-enter.'
      ELSE
        ierr=i_err_fatal
        RETURN
      ENDIF
    ENDDO
!
!   Convert to metres if in inverse cm or microns.
    IF (l_inverse_cm) THEN
      SpBasic%wavelength_short(i) = &
        1.0e-02_RealK / SpBasic%wavelength_short(i)
      SpBasic%wavelength_long(i) = &
        1.0e-02_RealK / SpBasic%wavelength_long(i)
    ELSE IF (l_micron) THEN
      SpBasic%wavelength_short(i) = &
        1.0e-06_RealK * SpBasic%wavelength_short(i)
      SpBasic%wavelength_long(i) = &
        1.0e-06_RealK * SpBasic%wavelength_long(i)
    ENDIF
!   Interchange the limits if they are in the wrong order.
    IF (SpBasic%wavelength_short(i) > SpBasic%wavelength_long(i)) THEN
      wavelength_temp = SpBasic%wavelength_short(i)
      SpBasic%wavelength_short(i) = SpBasic%wavelength_long(i)
      SpBasic%wavelength_long(i) = wavelength_temp
    ENDIF
  ENDDO
!
  WRITE(iu_stdout, '(A//)') 'All bands specified.'
!
!
!
  RETURN
END SUBROUTINE make_block_1
