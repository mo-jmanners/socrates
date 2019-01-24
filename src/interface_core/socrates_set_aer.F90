! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the variables in the Socrates aerosol type
!
!------------------------------------------------------------------------------
module socrates_set_aer
implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_SET_AER'
contains

subroutine set_aer(aer, control, dimen, spectrum, n_profile, n_layer)

use def_aer,      only: StrAer, allocate_aer, allocate_aer_prsc
use def_control,  only: StrCtrl
use def_dimen,    only: StrDim
use def_spectrum, only: StrSpecData

implicit none


! Aerosol properties:
type(StrAer),      intent(out) :: aer

! Control options:
type(StrCtrl),     intent(in)  :: control

! Dimensions:
type(StrDim),      intent(in)  :: dimen

! Spectral data:
type(StrSpecData), intent(in)  :: spectrum

integer, intent(in) :: n_profile
integer, intent(in) :: n_layer


! Allocate structure for the core radiation code interface
call allocate_aer(aer, dimen, spectrum)
call allocate_aer_prsc(aer, dimen, spectrum)

end subroutine set_aer
end module socrates_set_aer
