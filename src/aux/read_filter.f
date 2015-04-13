! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read a filter function.
!
! Method:
!	Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE read_filter(ierr
     &  , n_filter, wavenumber_filter, trans_filter)
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
      USE error_pcf
      USE dimensions_obs_ucf
      USE def_data_in_icf
!
!
      IMPLICIT NONE
!
!
!
!     Dummy arguments.

      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
     &  , n_filter
!           Number of filter points
      REAL  (RealK) ::
     &    wavenumber_filter(npd_filter)
!           Wavenumbers of filter points
     &  , trans_filter(npd_filter)
!           Transmission at filter points
!
!     Local variables
      CHARACTER
     &    line*80
!           Input string
     &  , file_filter*80
!           Name of filter file.
      INTEGER
     &    ios
!           I/O error flag
     &  , i
!           Loop variable
     &  , j
!           Loop variable
      LOGICAL
     &    l_micron
!           True if units are microns
     &  , l_metre
!           True if units are metres
     &  , l_nanometre
!           True if units are nanometres
     &  , l_inverse_cm
!           True if units are inverse cms.
!
!
      WRITE(iu_stdout, '(/a)')
     &  'Enter the name of the filter file.'
      READ(iu_stdin, '(a)') file_filter
      OPEN(unit=iu_filter, file=file_filter, iostat=ios, status='old')
      IF (ios /= 0) THEN
        WRITE(iu_err, '(/a)') '*** Error opening filter file.'
        ierr=i_err_fatal
        RETURN
      ENDIF
!
!     Read line until filter directive is found.
1     read(iu_filter, '(a)') line
      IF (line(1:7) /= '*FILTER') goto 1
!
!     Read header to find unit.
      READ(iu_filter, '(a)') line
      j=0
2       j=j+1
        IF (j <= 78) THEN
          IF (line(j:j) /= '(') goto 2
        ELSE
          WRITE(iu_err, '(/a)') 
     &      '*** Error: Line following *filter directive '
     &      //'is incorrect.'
        ENDIF
!
!     Initialize the unit flags to .FALSE.
      l_metre=.false.
      l_micron=.false.
      l_nanometre=.false.
      l_inverse_cm=.false.
!
      IF (line(j+1:j+1) == 'M') THEN
        l_metre=.true.
      ELSE IF (line(j+1:j+2) == 'UM') THEN
        l_micron=.true.
      ELSE IF (line(j+1:j+2) == 'NM') THEN
        l_nanometre=.true.
      ELSE IF (line(j+1:j+4) == 'CM-1') THEN
        l_inverse_cm=.true.
      ELSE
        WRITE(iu_err, '(/a)') '*** Error: Illegal unit.'
        ierr=i_err_fatal
        RETURN
      ENDIF
!
      n_filter=0
3     read(iu_filter, '(a)') line
      IF (line(1:4) /= '*END') THEN
        n_filter=n_filter+1
        backspace(iu_filter)
        READ(iu_filter, *, iostat=ios)
     &    wavenumber_filter(n_filter), trans_filter(n_filter)
        IF (ios /= 0) THEN
          WRITE(iu_err, '(/a)') '** Error: Unreadable filter given.'
          ierr=i_err_fatal
          RETURN
        ENDIF
        goto 3
      ENDIF
!
!     Convert the wavenumbers to inverse metres.
      IF (l_metre) THEN
        DO i=1, n_filter
          wavenumber_filter(i)=1.0_RealK/wavenumber_filter(i)
        ENDDO
      ELSE IF (l_micron) THEN
        DO i=1, n_filter
          wavenumber_filter(i)=1.0e+06_RealK/wavenumber_filter(i)
        ENDDO
      ELSE IF (l_nanometre) THEN
        DO i=1, n_filter
          wavenumber_filter(i)=1.0e+09_RealK/wavenumber_filter(i)
        ENDDO
      ELSE IF (l_inverse_cm) THEN
        DO i=1, n_filter
          wavenumber_filter(i)=1.0e+02_RealK*wavenumber_filter(i)
        ENDDO
      ENDIF
!
!
!
      RETURN
      END
