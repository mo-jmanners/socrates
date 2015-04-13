! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to tidy a spectral file.
!
PROGRAM tidy_90
!
! Description:
!   This program tidies a spectral file, removing undesirable features
!   and unnecessary contributions.
!
! Method:
!      A menu of operations on a spectral file is displayed. The
!      user selects those which are to be carried out. The
!      modified spectral file is written out.
!
! Modules used
  USE realtype_rd
  USE def_spectrum
  USE def_std_io_icf
  USE error_pcf
!
!
  IMPLICIT NONE
!
!
!
  CHARACTER (LEN=132) :: file_spectral
!   Name of the spectral file
  CHARACTER (LEN=1)   :: char_qna
!   Character response variable
  CHARACTER (LEN=1)   :: char_yn
!   Character response variable
  INTEGER :: ierr = 0
!   Error flag
  LOGICAL :: l_interactive
!   Flag for interactive operation
!
  INTEGER :: i_process
!   Process to be carried out
  INTEGER :: ios
!   Error flag
  INTEGER :: i
!   Loop variable
  INTEGER :: j
!   Loop variable
  INTEGER :: i_gas
!   Species of gas
  LOGICAL  :: l_exist
!   Existence flag
!
  TYPE(StrSpecData) :: Spectrum
!   The spectral configuration to be defined
!
!
!
! Read in the spectral file.
  WRITE(*, "(a)") "Enter the name of the spectral file."
  DO
    READ(*, "(a)") file_spectral
    CALL read_spectrum(file_spectral, Spectrum, ierr)
    IF (ierr == i_normal) THEN
      EXIT
    ELSE IF (l_interactive) THEN
      WRITE(*, "(a)") "Please re-specify"
      ierr=i_normal
    ELSE
      STOP
    ENDIF
  ENDDO
!
!
! Determine whether data are to be appended to the old file,
! or whether a new file is to be written.
  WRITE(iu_stdout, '(/a)') 'Type'
  WRITE(iu_stdout, '(6x, a)') &
    '"o" to overwrite the existing file;'
  WRITE(iu_stdout, '(6x, a)') &
    '"n" to create a new file;'
  WRITE(iu_stdout, '(6x, a)') &
    '"q" to quit.'
  WRITE(iu_stdout, '(/)')
!
  DO
    READ(iu_stdin, '(a)') char_qna
    SELECT CASE(char_qna)
      CASE('o', 'O')
        EXIT
      CASE('n', 'N')
        WRITE(iu_stdout, '(a)') 'Enter the name of the new file'
        DO
          READ(iu_stdin, '(a)') file_spectral
!         Check for existence of file.
          INQUIRE(file=file_spectral, exist=l_exist)
          IF (l_exist) THEN
            IF (l_interactive) THEN
              WRITE(iu_stdout, '(/a)') &
                'This file already exists: do you wish to overwrite? (y/n)'
              DO
                READ(iu_stdin, '(a)') char_yn
                SELECT CASE(char_yn)
                  CASE('n', 'N')
                    STOP
                  CASE('y', 'Y')
                    EXIT
                END SELECT
              ENDDO
            ELSE
              STOP
            ENDIF
          ELSE
            EXIT
          ENDIF
        ENDDO
      CASE('q', 'Q')
        STOP
      CASE DEFAULT
        WRITE(iu_err, '(a)') '+++ erroneous response:'
        IF (.NOT.l_interactive) STOP
    END SELECT
    EXIT
  ENDDO
!
! Now decide which processes are needed.
  DO
!
    WRITE(iu_stdout, '(/a)') &
      'Select from the following list of operations.'
!
    WRITE(iu_stdout, '(/6x, a)') &
      '1.   Remove gases from bands where they are weak.'
    WRITE(iu_stdout, '(6x, a)') &
      '2.   Force the sum of esft weights to be 1.'
    WRITE(iu_stdout, '(6x, a)') &
      '3.   Remove negative pressure scalings from esft data.'
!    WRITE(iu_stdout, '(6x, a)') &
!      '4.   Remove continua from bands where they are weak.'
    WRITE(iu_stdout, '(6x, a)') &
      '5.   Remove negative pressure scalings continuum data.'
    WRITE(iu_stdout, '(6x, a)') &
      '6.   Set the major gas in each band.'
    WRITE(iu_stdout, '(6x, a)') &
      '7.   Reorder the aerosols.'
    WRITE(iu_stdout, '(/5x, a//)') &
      '-1.  to finish.'
!
    DO
      READ(iu_stdin, *, iostat=ios) i_process
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ Erroneous input:'
        IF (.NOT.l_interactive) THEN
          STOP
        ELSE
          WRITE(iu_stdout, '(a)') 'Please re-type.'
        ENDIF
      ELSE
        EXIT
      ENDIF
    ENDDO
!
!
!   For each valid process call the appropriate routine.
    SELECT CASE(i_process)
      CASE(-1)
!       Write out the spectral file.
        CALL out_spectrum(file_spectral, Spectrum, ierr)
        STOP
      CASE(1)
        IF (Spectrum%Basic%l_present(5)) THEN
          CALL remove_weak
          IF (ierr /= i_normal) STOP
        ELSE
          WRITE(iu_err,'(/a)') &
            'The file contains no k-distribution data: the process ' &
            //'cannot be carried out.'
        ENDIF
      CASE(2)
        IF (Spectrum%Basic%l_present(5)) THEN
          DO i=1, Spectrum%Basic%n_band
            DO j=1, Spectrum%Gas%n_band_absorb(i)
              i_gas=Spectrum%Gas%index_absorb(j, i)
              CALL sum_unity(Spectrum%Gas%i_band_k(i, j), &
                Spectrum%Gas%w(1, i, j), Spectrum%Gas%k(1, i, j))
            ENDDO
          ENDDO
        ELSE
          WRITE(iu_err,'(/a)') &
            'The file contains no k-distribution data: the process ' &
            //'cannot be carried out.'
        ENDIF
      CASE(3)
        IF (Spectrum%Basic%l_present(5)) THEN
          CALL remove_negative_gas_90( &
            Spectrum%Dim%nd_band, &
            Spectrum%Dim%nd_species, &
            Spectrum%Dim%nd_k_term, &
            Spectrum%Dim%nd_scale_variable, &
            Spectrum%Basic%n_band, &
            Spectrum%Gas%n_band_absorb, Spectrum%Gas%index_absorb, &
            Spectrum%Gas%type_absorb, Spectrum%Gas%i_scale_fnc, &
            Spectrum%Gas%i_band_k, Spectrum%Gas%scale)
        ELSE
          WRITE(iu_err,'(/a)') &
            'The file contains no esft data: the process ' &
            //'cannot be carried out.'
        ENDIF
      CASE(4)
        IF (Spectrum%Basic%l_present(9)) THEN
!          CALL remove_weak_cont
          IF (ierr /= i_normal) STOP
      ELSE
        WRITE(iu_err,'(/a)') &
          'The file contains no continuum data: the process ' &
          //'cannot be carried out.'
      ENDIF
      CASE(5)
        IF (Spectrum%Basic%l_present(9)) THEN
          CALL remove_negative_cont_90( &
            Spectrum%Dim%nd_band, &
            Spectrum%Dim%nd_continuum, &
            Spectrum%Dim%nd_scale_variable, &
            Spectrum%Basic%n_band, &
            Spectrum%Cont%n_band_continuum, &
            Spectrum%Cont%index_continuum, &
            Spectrum%Cont%i_scale_fnc_cont, &
            Spectrum%Cont%scale_cont)
        ELSE
          WRITE(iu_err,'(/a)') &
            'The file contains no continuum data: the process ' &
            //'cannot be carried out.'
        ENDIF
      CASE(6)
        IF (Spectrum%Basic%l_present(5)) THEN
          CALL set_major_gas
          IF (ierr /= i_normal) STOP
        ELSE
          WRITE(iu_err,'(/a)') &
            'The file contains no k-distribution data: the process ' &
            //'cannot be carried out.'
        ENDIF
      CASE(7)
        IF (Spectrum%Basic%l_present(11)) THEN
          CALL reorder_aerosols
          IF (ierr /= i_normal) STOP
        ELSE
          WRITE(iu_err,'(/a)') &
            'The file contains no aerosol data: the process ' &
            //'cannot be carried out.'
        ENDIF        
      CASE DEFAULT
        WRITE(iu_err, '(a)') '+++ Invalid type of process:'
        IF (l_interactive) THEN
          WRITE(iu_stdout, '(a)') 'Please re-type.'
        ELSE
          STOP
        ENDIF
    END SELECT
    IF (ierr /= i_normal) stop
!
    WRITE(iu_stdout, '(/a/)') 'specify next process.'
!
  ENDDO

CONTAINS

!+ ---------------------------------------------------------------------
! Subroutine to remove weakly absorbing gases from bands.
!
! Method:
!	A tolerance for negelecting gaseous absorption is obtained.
!	A typical amount of each gas in a column is given. For each
!	band a test is made to see whether this amount of absorber
!	gives a transmission close to 1 as defined by the tolerance.
!	If so, the gas is removed from this spectral band.
!
!- ---------------------------------------------------------------------
  SUBROUTINE remove_weak

    USE gas_list_pcf

    LOGICAL, EXTERNAL :: lock_code
!     Logical to forbid interactive looping

    INTEGER :: K
!     Loop variable
    INTEGER :: N_BAND_ABSORB_TEMP
!     Temporary number of absorbers
    REAL (RealK) :: COLUMN(Spectrum%Dim%nd_species)
!     Column amounts for gases
    REAL (RealK) :: TRANS_NEGLECT
!     Transmission for neglecting
    REAL (RealK) :: TRANS_COLUMN
!     Transmission of column


!   Obtain column amounts of absorber.
    WRITE(IU_STDOUT, '(/A, /A)') &
      'For each absorber enter the amounts of the gas to test for', &
      'neglecting the gas in a band.'
    DO i=1, Spectrum%Gas%n_absorb
      WRITE(IU_STDOUT, '(4X, 2A)') 'Column amount for ', &
        name_absorb(Spectrum%Gas%type_absorb(i))
1     READ(IU_STDIN, *, IOSTAT=IOS) column(i)
      IF (IOS.NE.0) THEN
        WRITE(IU_ERR, '(A)') '+++ ERRONEOUS RESPONSE:'
        IF (LOCK_CODE(.TRUE.)) THEN
          STOP
        ELSE
          WRITE(IU_STDOUT, '(A)') 'PLEASE RE-TYPE.'
          GOTO 1
        ENDIF
      ENDIF
    ENDDO
!   
    WRITE(IU_STDOUT, '(/A)') &
      'Enter the transmission for neglecting gases.'
2   READ(IU_STDIN, *, IOSTAT=IOS) trans_neglect
    IF (IOS.NE.0) THEN
      WRITE(IU_ERR, '(A)') '+++ ERRONEOUS RESPONSE:'
      IF (LOCK_CODE(.TRUE.)) THEN
        STOP
      ELSE
        WRITE(IU_STDOUT, '(A)') 'PLEASE RE-TYPE.'
        GOTO 2
      ENDIF
    ENDIF
!   
!   Go through the bands removing gases which are too weak.
    DO i=1, Spectrum%Basic%n_band
      n_band_absorb_temp=Spectrum%Gas%n_band_absorb(i)
      Spectrum%Gas%n_band_absorb(i)=0
      DO j=1, n_band_absorb_temp
        i_gas=Spectrum%Gas%index_absorb(j, i)
        trans_column=0.0E+00_RealK
        DO k=1, Spectrum%Gas%i_band_k(i, i_gas)
          trans_column=trans_column + Spectrum%Gas%w(k, i, i_gas) &
            *EXP(-Spectrum%Gas%k(k, i, i_gas)*column(i_gas))
        ENDDO
        IF (trans_column < trans_neglect) THEN
!         The gas is included.
          Spectrum%Gas%n_band_absorb(i) = &
            Spectrum%Gas%n_band_absorb(i)+1
          Spectrum%Gas%index_absorb( &
            Spectrum%Gas%n_band_absorb(i), i) = i_gas
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE remove_weak


!+ ---------------------------------------------------------------------
! Subroutine to set the major gas in each band.
!
! Method:
!	A typical amount of each gas in a column is given. For each
!	band a test is made to see which gas is the most absorbing.
!	The index numbers are reordered to make this the major gas.
!
!- ---------------------------------------------------------------------
  SUBROUTINE set_major_gas

    USE gas_list_pcf

    LOGICAL, EXTERNAL :: lock_code
!     Logical to forbid interactive looping

    INTEGER :: K
!     Loop variable
    INTEGER :: MAJOR_GAS
!     Index number of major gas
    INTEGER :: N_MAJOR_GAS
!     Order number of major gas
    REAL (RealK) :: COLUMN(Spectrum%Dim%nd_species)
!     Column amounts for gases
    REAL (RealK) :: TRANS_MAJOR
!     Transmission for major gas
    REAL (RealK) :: TRANS_COLUMN
!     Transmission of column


!   Obtain column amounts of absorber.
    WRITE(IU_STDOUT, '(/A, /A)') &
      'For each absorber enter the amounts of the gas to find', &
      'the major gas in each band.'
    DO i=1, Spectrum%Gas%n_absorb
      WRITE(IU_STDOUT, '(4X, 2A)') 'Column amount for ', &
        name_absorb(Spectrum%Gas%type_absorb(i))
1     READ(IU_STDIN, *, IOSTAT=IOS) column(i)
      IF (IOS.NE.0) THEN
        WRITE(IU_ERR, '(A)') '+++ ERRONEOUS RESPONSE:'
        IF (LOCK_CODE(.TRUE.)) THEN
          STOP
        ELSE
          WRITE(IU_STDOUT, '(A)') 'PLEASE RE-TYPE.'
          GOTO 1
        ENDIF
      ENDIF
    ENDDO

!   Go through the bands to find the transmission for each gas
    DO i=1, Spectrum%Basic%n_band
      trans_major=1.0E+00_RealK
      DO J=1, Spectrum%Gas%n_band_absorb(i)
        I_GAS=Spectrum%Gas%index_absorb(J, I)
        trans_column=0.0E+00_RealK
        DO k=1, Spectrum%Gas%i_band_k(i, i_gas)
          trans_column=trans_column + Spectrum%Gas%w(k, i, i_gas) &
            *EXP(-Spectrum%Gas%k(k, i, i_gas)*column(i_gas))
        ENDDO
        IF (trans_column < trans_major) THEN
!         Major gas (so far).
          major_gas=i_gas
          n_major_gas=j
          trans_major=trans_column
        ENDIF
      ENDDO
!     Swap index of major gas with the first gas listed
      Spectrum%Gas%index_absorb(n_major_gas, i) = &
        Spectrum%Gas%index_absorb(1, i)
      Spectrum%Gas%index_absorb(1, i)=major_gas
    ENDDO

  END SUBROUTINE set_major_gas


!+ ---------------------------------------------------------------------
! Subroutine to reorder the aerosol indices
!- ---------------------------------------------------------------------
  SUBROUTINE reorder_aerosols

    INTEGER :: type_aerosol(Spectrum%Aerosol%n_aerosol)
    INTEGER :: map(Spectrum%Aerosol%n_aerosol), i, n

    WRITE(*, '(a)') 'Enter aerosol type numbers in the new order:'
    READ(*, *, IOSTAT=IOS) type_aerosol

    n=Spectrum%Aerosol%n_aerosol

!   Map aerosol indices onto new order
    DO i = 1, n
      map(i)=MINLOC(ABS(Spectrum%Aerosol%type_aerosol-type_aerosol(i)),1)
    END DO
    Spectrum%Aerosol%l_aero_spec(1:n)=Spectrum%Aerosol%l_aero_spec(map)
    Spectrum%Aerosol%type_aerosol(1:n)=Spectrum%Aerosol%type_aerosol(map)
    Spectrum%Aerosol%i_aerosol_parm(1:n)=Spectrum%Aerosol%i_aerosol_parm(map)
    Spectrum%Aerosol%n_aerosol_phf_term(1:n) = &
      Spectrum%Aerosol%n_aerosol_phf_term(map)
    Spectrum%Aerosol%nhumidity(1:n)=Spectrum%Aerosol%nhumidity(map)
    Spectrum%Aerosol%abs(:,1:n,:)=Spectrum%Aerosol%abs(:,map,:)
    Spectrum%Aerosol%scat(:,1:n,:)=Spectrum%Aerosol%scat(:,map,:)
    Spectrum%Aerosol%phf_fnc(:,:,1:n,:)=Spectrum%Aerosol%phf_fnc(:,:,map,:)
    Spectrum%Aerosol%humidities(:,1:n)=Spectrum%Aerosol%humidities(:,map)

    IF (Spectrum%Basic%l_present(15)) THEN
      Spectrum%Aerosol%i_aod_type(1:n)=Spectrum%Aerosol%i_aod_type(map)
      Spectrum%Aerosol%aod_abs(:,1:n,:)=Spectrum%Aerosol%aod_abs(:,map,:)
      Spectrum%Aerosol%aod_scat(:,1:n,:)=Spectrum%Aerosol%aod_scat(:,map,:)
    END IF

  END SUBROUTINE reorder_aerosols

END PROGRAM tidy_90
