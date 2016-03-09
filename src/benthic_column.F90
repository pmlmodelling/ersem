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
      type (type_horizontal_diagnostic_variable_id) :: id_poro,id_Dtot,id_EDZ_mix,id_layer2_thickness

      real(rk) :: d_tot
      real(rk) :: qPW,EDZ_mix
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

   type,extends(type_base_model),public :: type_ersem_bioturbation
      type (type_horizontal_diagnostic_variable_id) :: id_diff(3),id_diff_pom,id_Dtur
      type (type_horizontal_dependency_id)          :: id_biotur_tot, id_bioirr_tot

      real(rk) :: mtur, htur, Etur, dtur
      real(rk) :: mirr, hirr, irr_min, EDZ_1, EDZ_2, EDZ_3
   contains
      procedure :: initialize => bioturbation_initialize
      procedure :: do_bottom  => bioturbation_do_bottom
   end type

contains

   subroutine initialize(self,configunit)
      class (type_ersem_benthic_column),intent(inout),target :: self
      integer,                          intent(in)           :: configunit

      class (type_ersem_bioturbation), pointer :: bioturbation

      ! Set time unit to d-1. This implies that all rates (sink/source terms) are given in d-1.
      self%dt = 86400._rk

      call self%get_parameter(self%qPW,'qPW','-','sediment porosity')
      call self%get_parameter(self%EDZ_mix,'EDZ_mix','d/m','equilibrium diffusive speed between sediment surface water')
      call self%get_parameter(self%d_tot,'d_tot','m','depth of sediment column')

      call self%register_state_variable(self%id_D1m,'D1m','m','depth of bottom interface of oxygenated layer',standard_variable=depth_of_bottom_interface_of_layer_1)
      call self%register_state_variable(self%id_D2m,'D2m','m','depth of bottom interface of oxidized layer',standard_variable=depth_of_bottom_interface_of_layer_2)
      call self%register_diagnostic_variable(self%id_poro,'poro','-','porosity',standard_variable=sediment_porosity,missing_value=self%qPW,output=output_none,source=source_none)
      call self%register_diagnostic_variable(self%id_Dtot,'Dtot','m','depth of sediment column',missing_value=self%d_tot,standard_variable=depth_of_sediment_column,output=output_none,source=source_none)
      call self%register_diagnostic_variable(self%id_EDZ_mix,'cmix','d/m','equilibrium diffusive speed between sediment surface water',standard_variable=pelagic_benthic_transfer_constant,missing_value=self%EDZ_mix,output=output_none,source=source_none)
      call self%register_diagnostic_variable(self%id_layer2_thickness,'layer2_thickness','m','thickness of second layer',output=output_none,source=source_do_bottom)

      ! Create bioturbation submodel and provide it with parameters
      ! Currently the bioturbation logic must be separate from type_ersem_benthic_column to avoid circular dependencies.
      ! This is because type_ersem_benthic_column provides the max column depth, which is used to compute food for
      ! benthic fauna, which in turn results in the aggrege biturbation/bioirrigation activity.
      allocate(bioturbation)

      ! Bioturbation
      call self%get_parameter(bioturbation%Etur,'Etur','m^2/d','basal bioturbation rate')
      call self%get_parameter(bioturbation%mtur,'mtur','-','maximum relative turbation enhancement')
      call self%get_parameter(bioturbation%htur,'htur','mg C/m^2/d','Michaelis-Menten constant for bioturbation')
      call self%get_parameter(bioturbation%dtur,'dtur','m','bioturbation depth')

      ! Bioirrigation
      call self%get_parameter(bioturbation%EDZ_1,  'EDZ_1','m^2/d','diffusivity in oxygenated layer')
      call self%get_parameter(bioturbation%EDZ_2,  'EDZ_2','m^2/d','diffusivity in oxidized layer')
      call self%get_parameter(bioturbation%EDZ_3,  'EDZ_3','m^2/d','diffusivity in anoxic layer')
      call self%get_parameter(bioturbation%irr_min,'irr_min','-','minimum diffusion enhancement through bioirrigation')
      call self%get_parameter(bioturbation%mirr,   'mirr','-','maximum relative diffusion enhancement due to bioirrigation')
      call self%get_parameter(bioturbation%hirr,   'hirr','mg C/m^2/d','Michaelis-Menten constant for bioirrigation')

      ! Allow bioturbation module to initialize - must be after retrieving its parameters, as they are used when registering diagnostics.
      call self%add_child(bioturbation,'bioturbation',configunit=-1)

   end subroutine

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_benthic_column),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: D1m,D2m

      _HORIZONTAL_LOOP_BEGIN_

         ! Compute depth of second layer from the depths of the bottom interfaces of layers 1 and 2.
         _GET_HORIZONTAL_(self%id_D1m,D1m)
         _GET_HORIZONTAL_(self%id_D2m,D2m)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layer2_thickness,D2m-D1m)

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

   subroutine bioturbation_initialize(self,configunit)
      class (type_ersem_bioturbation),intent(inout),target :: self
      integer,                        intent(in)           :: configunit

      call self%register_diagnostic_variable(self%id_diff(1),'diff1','m^2/d','diffusivity in oxygenated layer',standard_variable=diffusivity_in_sediment_layer_1,source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_diff(2),'diff2','m^2/d','diffusivity in oxidized layer',standard_variable=diffusivity_in_sediment_layer_2,source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_diff(3),'diff3','m^2/d','diffusivity in anoxic layer',standard_variable=diffusivity_in_sediment_layer_3,source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_diff_pom,'diff_pom','m^2/d','particulate diffusivity representing bioturbation',standard_variable=particulate_diffusivity_due_to_bioturbation,source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_Dtur,'Dtur','m','bioturbation depth',standard_variable=bioturbation_depth,missing_value=self%dtur,output=output_none,source=source_none)

      ! Link to cumulative bioturbation and bioirrigation values, which account for activity of all benthic fauna.
      call self%register_dependency(self%id_biotur_tot,total_bioturbation_activity,domain=domain_bottom)
      call self%register_dependency(self%id_bioirr_tot,total_bioirrigation_activity,domain=domain_bottom)
   end subroutine bioturbation_initialize

   subroutine bioturbation_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_bioturbation),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: Ytur, Yirr, Irr_enh, Tur_enh

      _HORIZONTAL_LOOP_BEGIN_

         ! Compute "diffusivity of particulates", which represents bioturbation.
         _GET_HORIZONTAL_(self%id_biotur_tot,Ytur)
         Tur_enh = 1.0_rk + self%mtur * Ytur/(Ytur+self%htur)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_diff_pom,Tur_enh*self%Etur)

         ! Compute diffusivity of solutes that includes bioirrigation enhancement.
         _GET_HORIZONTAL_(self%id_bioirr_tot,Yirr)
         Irr_enh = self%irr_min + self%mirr * Yirr/(Yirr+self%hirr)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_diff(1),Irr_enh*self%EDZ_1)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_diff(2),Irr_enh*self%EDZ_2)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_diff(3),Irr_enh*self%EDZ_3)

      _HORIZONTAL_LOOP_END_

   end subroutine bioturbation_do_bottom

end module
