! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to caclulate weights for a filter function.
!
! Method:
!	Appropriate weights for each band in the spectral file
!	are calculted from a filter function which is read in.
!	splines are used to represent the function, so many
!	points should be used for functions which are not smooth,
!
!- ---------------------------------------------------------------------
      SUBROUTINE filter_function(ierr, n_band
     &  , wave_length_short, wave_length_long
     &  , weight_band
     &  , nd_band, nd_filter
     &  )
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE error_pcf
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
!     Sizes of dummy arrays:
      INTEGER, Intent(IN) ::
     &    nd_band
!           Size allocated for spectral bands
     &  , nd_filter
!           Size allocated for points in a filter function
!
      INTEGER, Intent(IN) ::
     &    n_band
!           Number of bands
      REAL  (RealK), Intent(IN) ::
     &    wave_length_short(nd_band)
!           Shortwave limits of each band
     &  , wave_length_long(nd_band)
!           Longwave limits of each band
      REAL  (RealK), Intent(OUT) ::
     &    weight_band(nd_band)
!           Weights for each band
!
!     Local variables
      INTEGER
     &    n_filter
!           Number of filter points
     &  , i
!           Loop variable
     &  , pointer(nd_filter)
!           Pointer to order filter
      REAL  (RealK) ::
     &    wavenumber_filter(nd_filter)
!           Wavenumbers of filter points
     &  , wavenumber_filter_ord(nd_filter)
!           Ordered wavenumbers of filter points
     &  , trans_filter(nd_filter)
!           Transmittance of filter
     &  , trans_filter_ord(nd_filter)
!           Ordered transmittance of filter
     &  , second_derivative(nd_filter)
!           Second derivative for spline
     &  , wavenumber_low
!           Low wavenumber of band
     &  , wavenumber_high
!           High wavenumber of band
!
!     Subroutines called:
      EXTERNAL
     &     read_filter, shell_sort, spline_fit, integrate_spline
!
!
!
!     Read the filter function
      CALL read_filter(ierr
     &  , n_filter, wavenumber_filter, trans_filter)
!
!     Order the points of the filter.
!     Initialize the pointer to the present order.
      DO i=1, n_filter
        pointer(i)=i
      ENDDO
      CALL shell_sort(n_filter, pointer, wavenumber_filter)
      DO i=1, n_filter
        wavenumber_filter_ord(i)=wavenumber_filter(pointer(i))
        trans_filter_ord(i)=trans_filter(pointer(i))
      ENDDO
!
!     Establish a spline fit to the filter function.
      CALL spline_fit(n_filter, wavenumber_filter_ord, trans_filter_ord
     &  , second_derivative)
!
!     Calculate the weighting function for each band by integrating
!     the spline fit
      DO i=1, n_band
        wavenumber_low=1.0_RealK/wave_length_long(i)
        wavenumber_high=1.0_RealK/wave_length_short(i)
        CALL integrate_spline(ierr
     &    , wavenumber_low, wavenumber_high
     &    , n_filter, wavenumber_filter_ord, trans_filter_ord
     &    , second_derivative, weight_band(i))
        IF (ierr /= i_normal) THEN
          IF (ierr == i_err_range) THEN
!           If the band lies outside the spline range the weighting 
!           is 0.We recover from this error.
            weight_band(i)=0.0_RealK
            ierr=i_normal
          ELSE
            ierr=i_err_fatal
          ENDIF
        ENDIF
        weight_band(i)=weight_band(i)/(wavenumber_high-wavenumber_low)
      ENDDO
!
!
!
      RETURN
      END
