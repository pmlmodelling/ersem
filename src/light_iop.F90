#include "fabm_driver.h"

module ersem_light_iop

   use fabm_types
   use ersem_shared

   implicit none

   private

   type,extends(type_base_model),public :: type_ersem_light_iop
      ! Identifiers for diagnostic variables
      type (type_diagnostic_variable_id)   :: id_EIR, id_parEIR, id_xEPS, id_secchi, id_iopABS, id_iopBBS
      type (type_dependency_id)            :: id_dz, id_abESS ,id_iopABSp, id_iopBBSp
      type (type_horizontal_dependency_id) :: id_I_0, id_zenithA

      ! Parameters
      real(rk) :: abESSX,a0w,b0w,pEIR_eowX
   contains
!     Model procedures
      procedure :: initialize
      procedure :: get_light
   end type type_ersem_light_iop

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_ersem_light_iop),intent(inout),target :: self
      integer,                     intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%a0w,    'a0w',   '1/m',   'absorption coefficient of clear water', default=.036_rk)
      call self%get_parameter(self%b0w,    'b0w',   '1/m',   'backscatter coefficient of clear water', default=.0016_rk)
      call self%get_parameter(self%pEIR_eowX,'pEIR_eow','-', 'photosynthetically active fraction of shortwave radiation', default=.5_rk)

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_EIR,'EIR','W/m^2','shortwave radiation', &
              standard_variable=standard_variables%downwelling_shortwave_flux,source=source_do_column)
      call self%register_diagnostic_variable(self%id_parEIR,'parEIR','W/m^2','photosynthetically active radiation', &
              standard_variable=standard_variables%downwelling_photosynthetic_radiative_flux,source=source_do_column)
      call self%register_diagnostic_variable(self%id_xEPS,'xEPS','1/m','attenuation coefficient of shortwave flux', &
              source=source_do_column)
      call self%register_diagnostic_variable(self%id_secchi,'secchi','m','Secchi depth (1.7/Kd)', &
              standard_variable=secchi_depth, source=source_do_column)
      call self%register_diagnostic_variable(self%id_iopABS,'iopABS','1/m','absorption coefficient of shortwave flux', &
              source=source_do_column)
      call self%register_diagnostic_variable(self%id_iopBBS,'iopBBS','1/m','backscatter coefficient of shortwave flux', &
              source=source_do_column)

      ! Register environmental dependencies (temperature, shortwave radiation)
      call self%register_dependency(self%id_I_0,standard_variables%surface_downwelling_shortwave_flux)
      call self%register_dependency(self%id_dz, standard_variables%cell_thickness)
      call self%register_dependency(self%id_iopABSp,particulate_organic_absorption_coefficient)
      call self%register_dependency(self%id_iopBBSp,particulate_organic_backscatter_coefficient)
      call self%register_dependency(self%id_abESS, type_bulk_standard_variable(name='absorption_of_silt'))
      call self%register_horizontal_dependency(self%id_zenithA, type_horizontal_standard_variable(name='zenith_angle'))
   end subroutine

   subroutine get_light(self,_ARGUMENTS_VERTICAL_)
      class (type_ersem_light_iop),intent(in) :: self
      _DECLARE_ARGUMENTS_VERTICAL_

      real(rk) :: buffer,dz,xEPS,iopABS,iopBBS,xtnc,EIR,abESS,zenithA
      real(rk),parameter :: bpk=.00022_rk

      _GET_HORIZONTAL_(self%id_I_0,buffer)
      _GET_HORIZONTAL_(self%id_zenithA,zenithA)   ! Zenith angle

      buffer = max(buffer, 0.0_rk)

      _VERTICAL_LOOP_BEGIN_
         _GET_(self%id_dz,dz)          ! Layer height (m)
         _GET_(self%id_iopABSp,iopABS) ! Absorption coefficient of shortwave radiation, due to particulate organic material (m-1)
         _GET_(self%id_iopBBSp,iopBBS) ! Backscatter coefficient of shortwave radiation, due to particulate organic material (m-1)
         _GET_(self%id_abESS,abESS)    ! Suspended silt absorption
         iopABS = iopABS+abESS+self%a0w
         iopBBS = iopBBS+bpk+self%b0w
         xEPS = (1._rk+.005_rk*zenithA)*iopABS+4.18_rk*(1._rk-.52_rk*exp(-10.8_rk*iopABS))*iopBBS
         xtnc = xEPS*dz
         EIR = buffer/xtnc*(1.0_rk-exp(-xtnc))  ! Note: this computes the vertical average, not the value at the layer centre.
         buffer = buffer*exp(-xtnc)
         _SET_DIAGNOSTIC_(self%id_EIR,EIR)                     ! Local shortwave radiation
         _SET_DIAGNOSTIC_(self%id_parEIR,EIR*self%pEIR_eowX)   ! Local photosynthetically active radiation
         _SET_DIAGNOSTIC_(self%id_xEPS,xEPS)                   ! Vertical attenuation of shortwave radiation
         _SET_DIAGNOSTIC_(self%id_secchi,1.7_rk/xEPS)          ! Secchi depth estimate as in Poole and Atkins (1929)
         _SET_DIAGNOSTIC_(self%id_iopABS,iopABS)               ! Vertical absorption of shortwave radiation
         _SET_DIAGNOSTIC_(self%id_iopBBS,iopBBS)               ! Vertical backscattering of shortwave radiation
      _VERTICAL_LOOP_END_

   end subroutine get_light

end module
