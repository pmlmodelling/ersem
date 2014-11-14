#include "fabm_driver.h"

module ersem_benthic_column

   use fabm_types

   use ersem_shared

   implicit none

   private

   ! Model that specifies the structure of the three-layer sediment column
   type,extends(type_base_model),public :: type_ersem_benthic_column
      type (type_bottom_state_variable_id)          :: id_D1m,id_D2m
      type (type_horizontal_diagnostic_variable_id) :: id_poro,id_Dtot,id_diff(3),id_EDZ_mixX,id_layer2_thickness

      real(rk) :: EDZ_1X,EDZ_2X,EDZ_3X,d_totX
      real(rk) :: qPWX,EDZ_mixX
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self,configunit)
      class (type_ersem_benthic_column),intent(inout),target :: self
      integer,                                     intent(in)           :: configunit
      
      ! Set time unit to d-1. This implies that all rates (sink/source terms) are given in d-1.
      self%dt = 86400._rk

      call self%get_parameter(self%EDZ_1X,'EDZ_1X','m^2/d','diffusivity in 1st (oxygenated) layer')
      call self%get_parameter(self%EDZ_2X,'EDZ_2X','m^2/d','diffusivity in 2nd (oxidized) layer')
      call self%get_parameter(self%EDZ_3X,'EDZ_3X','m^2/d','diffusivity in 3rd (anoxic) layer')
      call self%get_parameter(self%qPWX,'qPWX','-','fraction of pore water in the sediment')
      call self%get_parameter(self%EDZ_mixX,'EDZ_mixX','d/m','equilibrium diffusive speed between sediment surface water')
      call self%get_parameter(self%d_totX,'d_totX','m','depth of sediment column')

      call self%register_state_variable(self%id_D1m,'D1m','m','depth of bottom interface of 1st layer',standard_variable=depth_of_bottom_interface_of_layer_1)
      call self%register_state_variable(self%id_D2m,'D2m','m','depth of bottom interface of 2nd layer',standard_variable=depth_of_bottom_interface_of_layer_2)
      call self%register_diagnostic_variable(self%id_poro,'poro','-','porosity',standard_variable=sediment_porosity,missing_value=self%qPWX)
      call self%register_diagnostic_variable(self%id_Dtot,'Dtot','m','depth of sediment column',missing_value=self%d_totX,standard_variable=depth_of_sediment_column)
      call self%register_diagnostic_variable(self%id_diff(1),'diff1','m^2/s','diffusivity in layer 1',standard_variable=diffusivity_in_sediment_layer_1,missing_value=self%EDZ_1X)
      call self%register_diagnostic_variable(self%id_diff(2),'diff2','m^2/s','diffusivity in layer 2',standard_variable=diffusivity_in_sediment_layer_2,missing_value=self%EDZ_2X)
      call self%register_diagnostic_variable(self%id_diff(3),'diff3','m^2/s','diffusivity in layer 3',standard_variable=diffusivity_in_sediment_layer_3,missing_value=self%EDZ_3X)
      call self%register_diagnostic_variable(self%id_EDZ_mixX,'cmix','s/m','equilibrium diffusive speed between sediment surface water',standard_variable=pelagic_benthic_transfer_constant,missing_value=self%EDZ_mixX)
      call self%register_diagnostic_variable(self%id_layer2_thickness,'layer2_thickness','m','thickness of second layer',output=output_none)
   end subroutine

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_benthic_column),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: D1m,D2m   

      _HORIZONTAL_LOOP_BEGIN_
          
      _GET_HORIZONTAL_(self%id_D1m,D1m)
      _GET_HORIZONTAL_(self%id_D2m,D2m)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layer2_thickness,D2m-D1m)

      _HORIZONTAL_LOOP_END_
   end subroutine do_bottom

end module
