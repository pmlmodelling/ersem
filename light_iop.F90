#include "fabm_driver.h"

module ersem_light_iop

   use fabm_types
   use ersem_shared

   implicit none

   private

   type,extends(type_base_model),public :: type_ersem_light_iop
      ! Identifiers for diagnostic variables
      type (type_diagnostic_variable_id)   :: id_EIR, id_parEIR, id_xEPS, id_iopADS, id_iopBBS
      type (type_dependency_id)            :: id_dz, id_adESS ,id_iopADSp, id_iopBBSp
      type (type_horizontal_dependency_id) :: id_I_0, id_zenithA

      ! Parameters
      real(rk) :: adESSX,a0w,b0w,pEIR_eowX
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
      call self%get_parameter(self%a0w,    'a0w',   '1/m',   'adsorption coefficient of clear water') 
      call self%get_parameter(self%b0w,    'b0w',   '1/m',   'backscatter coefficient of clear water') 
      call self%get_parameter(self%pEIR_eowX,'pEIR_eow','-',     'photosynthetically active fraction of shortwave radiation') 

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_EIR,'EIR','W/m^2','shortwave radiation', &
              standard_variable=standard_variables%downwelling_shortwave_flux,source=source_do_column)
      call self%register_diagnostic_variable(self%id_parEIR,'parEIR','W/m^2','photosynthetically active radiation', &
              standard_variable=standard_variables%downwelling_photosynthetic_radiative_flux,source=source_do_column)
      call self%register_diagnostic_variable(self%id_xEPS,'xEPS','1/m','attenuation coefficient of shortwave flux', &
              source=source_do_column)
      call self%register_diagnostic_variable(self%id_iopADS,'iopADS','1/m','adsorption coefficient of shortwave flux', &
              source=source_do_column)
      call self%register_diagnostic_variable(self%id_iopBBS,'iopBBS','1/m','backscatter coefficient of shortwave flux', &
              source=source_do_column)

      ! Register environmental dependencies (temperature, shortwave radiation)
      call self%register_dependency(self%id_I_0,standard_variables%surface_downwelling_shortwave_flux)
      call self%register_dependency(self%id_dz, standard_variables%cell_thickness)
      call self%register_dependency(self%id_iopADSp,particulate_organic_adsportion_coefficient)
      call self%register_dependency(self%id_iopBBSp,particulate_organic_backscatter_coefficient)
      call self%register_dependency(self%id_adESS, type_bulk_standard_variable(name='adsorption_of_silt'))
      call self%register_horizontal_dependency(self%id_zenithA, type_horizontal_standard_variable(name='zenith_angle')) 
   end subroutine
   
   subroutine get_light(self,_ARGUMENTS_VERT_)
      class (type_ersem_light_iop),intent(in) :: self
      _DECLARE_ARGUMENTS_VERT_

      real(rk) :: buffer,dz,xEPS,iopADS,iopBBS,xtnc,EIR,adESS,zenithA
      real(rk),parameter :: bpk=.00022_rk

      _GET_HORIZONTAL_(self%id_I_0,buffer)
      _VERTICAL_LOOP_BEGIN_
         _GET_(self%id_dz,dz)     ! Layer height (m)
         _GET_(self%id_iopADSp,iopADS) ! Adsorption coefficient of shortwave radiation, due to particulate organic material (m-1)
         _GET_(self%id_iopBBSp,iopBBS) ! Backscatter coefficient of shortwave radiation, due to particulate organic material (m-1)
         _GET_(self%id_adESS,adESS)   ! Suspended silt adsorption
         _GET_HORIZONTAL_(self%id_zenithA,zenithA)   ! Zenith angle
         iopADS = iopADS+adESS+self%a0w
         iopBBS = iopBBS+bpk+self%b0w
         xEPS = (1._rk+.005_rk*zenithA)*iopADS+4.18_rk*(1._rk-.52_rk*exp(-10.8_rk*iopADS))*iopBBS
         xtnc = xEPS*dz
         EIR = buffer/xtnc*(1.0_rk-exp(-xtnc))  ! Note: this computes the vertical average, not the value at the layer centre.
         buffer = buffer*exp(-xtnc)
         _SET_DIAGNOSTIC_(self%id_EIR,EIR)                     ! Local shortwave radiation
         _SET_DIAGNOSTIC_(self%id_parEIR,EIR*self%pEIR_eowX)   ! Local photosynthetically active radiation
         _SET_DIAGNOSTIC_(self%id_xEPS,xEPS)                   ! Vertical attenuation of shortwave radiation
         _SET_DIAGNOSTIC_(self%id_iopADS,iopADS)                   ! Vertical attenuation of shortwave radiation
         _SET_DIAGNOSTIC_(self%id_iopBBS,iopBBS)                   ! Vertical attenuation of shortwave radiation
      _VERTICAL_LOOP_END_

   end subroutine get_light

end module
