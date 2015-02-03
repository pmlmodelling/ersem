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
      type (type_horizontal_dependency_id) :: id_I_0

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
      call self%get_parameter(self%a0w,    'a0w',   '1/m',   'adsorption coefficient of clear water [1/m]') 
      call self%get_parameter(self%b0w,    'b0w',   '1/m',   'backscatter coefficient of clear water [1/m]') 
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
      !call self%register_dependency(self%zenithA), ??? ) 
   end subroutine
   
   subroutine get_light(self,_ARGUMENTS_VERT_)
      class (type_ersem_light_iop),intent(in) :: self
      _DECLARE_ARGUMENTS_VERT_

      real(rk) :: buffer,dz,xEPS,iopADS,iopBBS,xtnc,EIR,adESS,zenithA

      _GET_HORIZONTAL_(self%id_I_0,buffer)
      _VERTICAL_LOOP_BEGIN_
         _GET_(self%id_dz,dz)     ! Layer height (m)
         _GET_(self%id_iopADSp,iopADS) ! Adsorption coefficient of shortwave radiation, due to particulate organic material (m-1)
         _GET_(self%id_iopBBSp,iopBBS) ! Backscatter coefficient of shortwave radiation, due to particulate organic material (m-1)
         _GET_(self%id_adESS,adESS)   ! Suspended silt adsorption
         !_GET_(self%id_zenithA,zenithA)   ! Zenith angle
         xEPS = (1.+.005*zenithA)*iopADS+4.18*(1.-.52*exp(-10.8*iopADS))*iopBBS
         xtnc = xEPS*dz
         EIR = buffer/xtnc*(1.0_rk-exp(-xtnc))  ! Note: this computes the vertical average, not the value at the layer centre.
         buffer = buffer*exp(-xtnc)
         _SET_DIAGNOSTIC_(self%id_EIR,EIR)                     ! Local shortwave radiation
         _SET_DIAGNOSTIC_(self%id_parEIR,EIR*self%pEIR_eowX)   ! Local photosynthetically active radiation
         _SET_DIAGNOSTIC_(self%id_xEPS,xEPS)                   ! Vertical attenuation of shortwave radiation
         _SET_DIAGNOSTIC_(self%id_iopADS,iopADS)                   ! Adsorption of shortwave radiation
         _SET_DIAGNOSTIC_(self%id_iopBBS,iopBBS)                   ! Backscatter of shortwave radiation
      _VERTICAL_LOOP_END_

   end subroutine get_light

end module
