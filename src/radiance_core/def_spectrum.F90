! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to declare a structure of spectral data.
!
! Description:
!   This module contains the heirarchical declaration of structures
!   of spectral data.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiance Core
!
!------------------------------------------------------------------------------
! CAUTION - Any changes made to this routine need to be mirrored in the
!           setup_spectra_mod module in the UM.
!------------------------------------------------------------------------------
MODULE def_spectrum

USE realtype_rd

IMPLICIT NONE


TYPE StrSpecDim
  INTEGER :: nd_type
!   Size allocated for spectral blocks
  INTEGER :: nd_band
!   Size allocated for spectral bands
  INTEGER :: nd_exclude
!   Size allocated for excluded bands
  INTEGER :: nd_k_term
!   Size allocated for k-terms
  INTEGER :: nd_species
!   Size allocated for gaseous species
  INTEGER :: nd_scale_variable
!   Size allocated for scaling variables
  INTEGER :: nd_continuum
!   Size allocated for continua
  INTEGER :: nd_drop_type
!   Size allocated for drop types
  INTEGER :: nd_ice_type
!   Size allocated for ice crystal types
  INTEGER :: nd_aerosol_species
!   Size allocated for aerosol species
  INTEGER :: nd_aerosol_mr
!   Size allocated for aerosol mixing ratios
  INTEGER :: nd_thermal_coeff
!   Size allocated for thermal coefficients
  INTEGER :: nd_cloud_parameter
!   Size allocated for cloud parameters
  INTEGER :: nd_humidity
!   Size allocated for humidities
  INTEGER :: nd_aod_wavel
!   Number of wavelengths for aerosol optical depths
  INTEGER :: nd_phase_term
!   Size allocated for terms in the phase function
  INTEGER :: nd_tmp
!   Number of reference temperature for k-terms
  INTEGER :: nd_pre
!   Number of reference pressures for k-terms
  INTEGER :: nd_mix
!   Number of eta for mixture absorbing species
  INTEGER :: nd_band_mix_gas
!   Number of bands where mixed species exist
END TYPE StrSPecDim


TYPE StrSpecBasic
  LOGICAL, ALLOCATABLE      :: l_present(:)
!   Blocks of spectral data in the file
  INTEGER                   :: n_band
!   Number of Spectral Band used
  REAL (RealK), ALLOCATABLE :: wavelength_long(:)
!   Lower wavelength limits for the band
  REAL (RealK), ALLOCATABLE :: wavelength_short(:)
!   Higher wavelengths limits for the band
  INTEGER, ALLOCATABLE      :: n_band_exclude(:)
!   Number of exclusions from each band
  INTEGER, ALLOCATABLE      :: index_exclude(:, :)
!   List of excluded bands within each region
END TYPE StrSpecBasic


TYPE StrSpecSolar
  REAL (RealK), ALLOCATABLE :: solar_flux_band(:)
!   Fraction of the solar spectrum in each band
  REAL (RealK), ALLOCATABLE :: solar_flux_band_ses(:, :)
!   Fraction of the solar spectrum for each k-term
  REAL (RealK), ALLOCATABLE :: weight_blue(:)
!   Fraction of the surface flux designated as "blue" in each band
END TYPE StrSpecSolar


TYPE StrSpecRayleigh
  REAL (RealK), ALLOCATABLE :: rayleigh_coeff(:)
!   Rayleigh scattering coefficients in each band
END TYPE StrSpecRayleigh


TYPE StrSpecGas
  INTEGER  :: n_absorb
!   Total number of gaseous absorbers
  INTEGER, ALLOCATABLE      :: n_band_absorb(:)
!   Number of gaseous absorbers in each band
  INTEGER, ALLOCATABLE      :: index_absorb(:, :)
!   Number of gaseous absorbers
  INTEGER, ALLOCATABLE      :: type_absorb(:)
!   Actual types of each gas in the spectral file
  INTEGER, ALLOCATABLE      :: n_mix_gas(:)
!   Number of mixed gases in a band
  INTEGER, ALLOCATABLE      :: index_mix_gas(:, :)
!   Index of mixed absorbers in each band
  INTEGER, ALLOCATABLE      :: num_mix(:)
!   Number of binary parameter for interpolation of absorption
!   coefficient for mixture of two species
  INTEGER, ALLOCATABLE      :: mix_gas_band(:)
!   Sequence band number (not real band number) of mixed species
  INTEGER, ALLOCATABLE      :: num_ref_p(:, :)
!   Number of reference pressures
  INTEGER, ALLOCATABLE      :: num_ref_t(:, :)
!   Number of reference temperatures
  INTEGER, ALLOCATABLE      :: i_band_k(:, :)
!   Number of k-terms in each band for each gas
  INTEGER, ALLOCATABLE      :: i_band_k_ses(:)
!   Number of k-terms in band for each gas
  INTEGER, ALLOCATABLE      :: i_scale_k(:, :)
!   Type of scaling applied to each k-term
  INTEGER, ALLOCATABLE      :: i_scale_fnc(:, :)
!   Type of scaling function
  INTEGER, ALLOCATABLE      :: i_scat(:, :, :)
!   Method of scattering treatment for each k-term

  REAL (RealK), ALLOCATABLE :: k(:, :, :)
!   Absorption coefficients of k-terms
  REAL (RealK), ALLOCATABLE :: w(:, :, :)
!   Weights for k-terms
  REAL (RealK), ALLOCATABLE :: scale(:, :, :, :)
!   Scaling parameters for each absorber and term
  REAL (RealK), ALLOCATABLE :: p_ref(:, :)
!   Reference pressures for scaling functions
  REAL (RealK), ALLOCATABLE :: t_ref(:, :)
!   Reference temperatures for scaling functions

  REAL (RealK), ALLOCATABLE :: p_lookup(:)
  REAL (RealK), ALLOCATABLE :: t_lookup(:, :)
  REAL (RealK), ALLOCATABLE :: k_lookup(:, :, :, :, :)
  REAL (RealK), ALLOCATABLE :: w_ses(:, :)
  REAL (RealK), ALLOCATABLE :: k_mix_gas(:, :, :, :, :)
!   Absorption coefficients for mixture species
  REAL (RealK), ALLOCATABLE :: f_mix(:)
!   Mixing ratio of mixed absorber amount

  LOGICAL, ALLOCATABLE      :: l_doppler(:)
!   Flag for Doppler broadening for each species
  REAL (RealK), ALLOCATABLE :: doppler_cor(:)
!   Doppler correction terms
END TYPE StrSpecGas


TYPE StrSpecPlanck
  INTEGER                   :: n_deg_fit
!   Degree of the fit to the Planckian function
  REAL (RealK), ALLOCATABLE :: thermal_coeff(:, :)
!   Coefficients in polynomial fit to source function
  REAL (RealK), ALLOCATABLE :: theta_planck_tbl(:)
!   Temperatures at which the band-integrated Planck function
!   has been evaluated.
  REAL (RealK)              :: t_ref_planck
!   Reference temperature for the Plackian function
  LOGICAL                   :: l_planck_tbl
!   Flag for using a look-up table instead of a polynomial
END TYPE StrSpecPlanck


TYPE StrSpecCont
  INTEGER, ALLOCATABLE      :: n_band_continuum(:)
!   Number of continua in each band
  INTEGER, ALLOCATABLE      :: index_continuum(:, :)
!   List of continua in each band
  INTEGER                   :: index_water
!   Index of water vapour of continua in each band
  INTEGER, ALLOCATABLE      :: i_scale_fnc_cont(:, :)
!   Types of scaling functions for continua

  REAL (RealK), ALLOCATABLE :: k_cont(:, :)
!   ABsorption coefficients for continuum absorption
  REAL (RealK), ALLOCATABLE :: scale_cont(:, :, :)
!   Reference temperature for the Plackian function
  REAL (RealK), ALLOCATABLE :: p_ref_cont(:, :)
!   Reference pressures for continuum scaling functions
  REAL (RealK), ALLOCATABLE :: t_ref_cont(:, :)
!   Reference temperatures for continuum scaling functions
  REAL (RealK), ALLOCATABLE :: k_cont_ses(:, :, :, :)
  REAL (RealK), ALLOCATABLE :: k_h2oc(:, :, :, :)
!   Absorption coefficient for water vapour continuum
END TYPE StrSpecCont


TYPE StrSpecDrop
  LOGICAL, ALLOCATABLE      :: l_drop_type(:)
!   Flags for types of droplets present
  INTEGER, ALLOCATABLE      :: i_drop_parm(:)
!   Form of parametrization for each type of droplet
  INTEGER, ALLOCATABLE      :: n_phf(:)
!   Number of moments of the phase fuction fitted (N. B. This
!   array is not set for parametrizations which are implicitly
!   restricted to the asymmetry.)

  REAL (RealK), ALLOCATABLE :: parm_list(:, :, :)
!   Parameters used to fit the optical properties of droplets
  REAL (RealK), ALLOCATABLE :: parm_min_dim(:)
!   Minimum dimension permissible in the parametrization
  REAL (RealK), ALLOCATABLE :: parm_max_dim(:)
!   Maximum dimension permissible in the parametrization
END TYPE StrSpecDrop


TYPE StrSpecAerosol
  LOGICAL, ALLOCATABLE      :: l_aero_spec(:)
!   Flags for species of aerosol present

  INTEGER                   :: n_aerosol
!   Number of aerosol species present in spectral file
  INTEGER                   :: n_aerosol_mr
!   Number of aerosol species present in mixing ratio array
  INTEGER, ALLOCATABLE      :: type_aerosol(:)
!   Actual types of aerosols in the spectral file
  INTEGER, ALLOCATABLE      :: i_aerosol_parm(:)
!   Parametrization scheme used for each aerosol
  INTEGER, ALLOCATABLE      :: n_aerosol_phf_term(:)
!   Number of terms in the phase function
  INTEGER, ALLOCATABLE      :: nhumidity(:)
!   Number of values of humidity

  REAL (RealK), ALLOCATABLE :: abs(:, :, :)
!   Absortption by aerosols
  REAL (RealK), ALLOCATABLE :: scat(:, :, :)
!   Scattering by aerosols
  REAL (RealK), ALLOCATABLE :: phf_fnc(:, :, :, :)
!   Phase functions of aerosols
  REAL (RealK), ALLOCATABLE :: humidities(:, :)
!   Humdities of each component

! Fields for aerosol optical depth:
  INTEGER                   :: n_aod_wavel
!   Number of wavelengths
  INTEGER, ALLOCATABLE      :: i_aod_type(:)
!   Relationship between aerosol component and type
  REAL (RealK), ALLOCATABLE :: aod_wavel(:)
!   Wavelengths for the aod
  REAL (RealK), ALLOCATABLE :: aod_abs(:, :, :)
!   Monochromatic specific absorption coefficient
  REAL (RealK), ALLOCATABLE :: aod_scat(:, :, :)
!   Monochromatic specific scattering coefficient
END TYPE StrSpecAerosol


TYPE StrSpecIce
  LOGICAL, ALLOCATABLE      :: l_ice_type(:)
!   Flags for types of ice crystals present
  INTEGER, ALLOCATABLE      :: i_ice_parm(:)
!   Form of parametrization for each type of ice crystal
  INTEGER, ALLOCATABLE      :: n_phf(:)
!   Number of moments of the phase fuction fitted

  REAL (RealK), ALLOCATABLE :: parm_list(:, :, :)
!   Parameters used to fit the optical properties of ice crystals
  REAL (RealK), ALLOCATABLE :: parm_min_dim(:)
!   Minimum dimension permissible in the parametrization
  REAL (RealK), ALLOCATABLE :: parm_max_dim(:)
!   Maximum dimension permissible in the parametrization
END TYPE StrSpecIce


TYPE StrSpecData
  TYPE (StrSpecDim)               :: Dim
  TYPE (StrSpecBasic)             :: Basic
  TYPE (StrSpecSolar)             :: Solar
  TYPE (StrSpecRayleigh)          :: Rayleigh
  TYPE (StrSpecGas)               :: Gas
  TYPE (StrSpecPlanck)            :: Planck
  TYPE (StrSpecCont)              :: Cont
  TYPE (StrSpecDrop)              :: Drop
  TYPE (StrSpecAerosol)           :: Aerosol
  TYPE (StrSpecIce)               :: Ice
END TYPE StrSpecData

END MODULE def_spectrum
