! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! @brief Define the Socrates input type

module socrates_def_input

use realtype_rd, only: RealExt
use gas_list_pcf, only: npd_gases

implicit none

type :: StrInputGas

real(RealExt) :: scalar_mass_mix_ratio = 0.0_RealExt
! Mass mixing ratio if well mixed, kg/kg

real(RealExt), pointer :: mass_mix_ratio(:,:) => null()
! Mass mixing ratio, kg/kg (n_profile, n_layer)

end type StrInputGas


type :: StrInput
  type (StrInputGas) :: gas(npd_gases)
end type StrInput

end module socrates_def_input
