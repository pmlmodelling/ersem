#include "fabm_driver.h"

! type_ersem_benthic_column
!
! This model specifies the structure of the three-layer sediment column.
!
! This model also computes the diffusivity of solutes in the different layers,
! and a "particulate diffusivity" that represents bioturbation. These variables
! account for variable bioturbation and bioirrigation activity, respectively.

module ersem_benthic_column

   use fabm_types

   use ersem_shared

   implicit none

   private

   type,extends(type_base_model),public :: type_ersem_benthic_column
      type (type_bottom_state_variable_id)          :: id_D1m,id_D2m
      type (type_horizontal_diagnostic_variable_id) :: id_poro,id_Dtot,id_diff(3),id_EDZ_mixX,id_diff_pom,id_layer2_thickness
      type (type_horizontal_dependency_id)          :: id_biotur_tot, id_bioirr_tot

      real(rk) :: d_totX
      real(rk) :: qPWX,EDZ_mixX
      real(rk) :: mturX, hturX, EturX
      real(rk) :: mirrX, hirrX, irr_minX, EDZ_1X, EDZ_2X, EDZ_3X
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self,configunit)
      class (type_ersem_benthic_column),intent(inout),target :: self
      integer,                          intent(in)           :: configunit
      
      ! Set time unit to d-1. This implies that all rates (sink/source terms) are given in d-1.
      self%dt = 86400._rk

      call self%get_parameter(self%qPWX,'qPWX','-','fraction of pore water in the sediment')
      call self%get_parameter(self%EDZ_mixX,'EDZ_mixX','d/m','equilibrium diffusive speed between sediment surface water')
      call self%get_parameter(self%d_totX,'d_totX','m','depth of sediment column')

      ! Bioturbation
      call self%get_parameter(self%mturX,'mturX','-','maximum relative turbation enhancement')
      call self%get_parameter(self%hturX,'hturX','mg C/m^2/d','Michaelis-Menten constant for bioturbation')
      call self%get_parameter(self%EturX,'EturX','m^2/d','basal bioturbation rate')

      ! Bioirrigation
      call self%get_parameter(self%mirrX,   'mirrX','-','maximum relative diffusion enhancement due to bioirrigation')
      call self%get_parameter(self%hirrX,   'hirrX','mg C/m^2/d','Michaelis-Menten constant for bioirrigation')
      call self%get_parameter(self%irr_minX,'irr_minX','-','minimum diffusion enhancement through bioirrigation')
      call self%get_parameter(self%EDZ_1X,  'EDZ_1X','m^2/d','diffusivity in 1st (oxygenated) layer')
      call self%get_parameter(self%EDZ_2X,  'EDZ_2X','m^2/d','diffusivity in 2nd (oxidized) layer')
      call self%get_parameter(self%EDZ_3X,  'EDZ_3X','m^2/d','diffusivity in 3rd (anoxic) layer')

      call self%register_state_variable(self%id_D1m,'D1m','m','depth of bottom interface of 1st layer',standard_variable=depth_of_bottom_interface_of_layer_1)
      call self%register_state_variable(self%id_D2m,'D2m','m','depth of bottom interface of 2nd layer',standard_variable=depth_of_bottom_interface_of_layer_2)
      call self%register_diagnostic_variable(self%id_poro,'poro','-','porosity',standard_variable=sediment_porosity,missing_value=self%qPWX)
      call self%register_diagnostic_variable(self%id_Dtot,'Dtot','m','depth of sediment column',missing_value=self%d_totX,standard_variable=depth_of_sediment_column)
      call self%register_diagnostic_variable(self%id_diff(1),'diff1','m^2/d','diffusivity in layer 1',standard_variable=diffusivity_in_sediment_layer_1)
      call self%register_diagnostic_variable(self%id_diff(2),'diff2','m^2/d','diffusivity in layer 2',standard_variable=diffusivity_in_sediment_layer_2)
      call self%register_diagnostic_variable(self%id_diff(3),'diff3','m^2/d','diffusivity in layer 3',standard_variable=diffusivity_in_sediment_layer_3)
      call self%register_diagnostic_variable(self%id_diff_pom,'diff_pom','m^2/d','particulate diffusivity representing bioturbation',standard_variable=particulate_diffusivity_representing_bioturbation)
      call self%register_diagnostic_variable(self%id_EDZ_mixX,'cmix','s/m','equilibrium diffusive speed between sediment surface water',standard_variable=pelagic_benthic_transfer_constant,missing_value=self%EDZ_mixX)
      call self%register_diagnostic_variable(self%id_layer2_thickness,'layer2_thickness','m','thickness of second layer',output=output_none)

      ! Link to cumulative bioturbation and bioirrigation values, which account for activity of all benthic fauna.
      !call self%register_dependency(self%id_biotur_tot,total_bioturbation_activity)
      !call self%register_dependency(self%id_bioirr_tot,total_bioirrigation_activity)
   end subroutine

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_benthic_column),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: D1m,D2m
      real(rk) :: Ytur, Yirr, Irr_enh, Tur_enh

      _HORIZONTAL_LOOP_BEGIN_

         ! Compute "diffusivity of particulates", which represents bioturbation.
         _GET_HORIZONTAL_(self%id_biotur_tot,Ytur)
         Ytur = 0.0_rk   ! Temporary: link to standard variable not working yet
         Tur_enh = 1.0_rk + self%mturX * Ytur/(Ytur+self%hturX)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_diff_pom,Tur_enh*self%EturX)

         ! Compute diffusivity of solutes that includes bioirrigation enhancement.
         _GET_HORIZONTAL_(self%id_bioirr_tot,Yirr)
         Yirr = 0.0_rk   ! Temporary: link to standard variable not working yet
         Irr_enh = self%irr_minX + self%mirrX * Yirr/(Yirr+self%hirrX)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_diff(1),Irr_enh*self%EDZ_1X)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_diff(2),Irr_enh*self%EDZ_2X)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_diff(3),Irr_enh*self%EDZ_3X)

         ! Compute depth of second layer from the depths of the bottom interfaces of layers 1 and 2.
         _GET_HORIZONTAL_(self%id_D1m,D1m)
         _GET_HORIZONTAL_(self%id_D2m,D2m)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layer2_thickness,D2m-D1m)

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module
