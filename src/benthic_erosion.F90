#include "fabm_driver.h"

! type_ersem_benthic_column
!
! This model specifies the structure of the three-layer sediment column.
!
! This model also computes the diffusivity of solutes in the different layers,
! and a "particulate diffusivity" that represents bioturbation. These variables
! account for variable bioturbation and bioirrigation activity, respectively.

module ersem_benthic_erosion

   use fabm_types

   use ersem_shared

   implicit none

   private

   type,extends(type_base_model),public :: type_ersem_benthic_erosion
      type (type_dependency_id)                     :: id_dens
      type (type_horizontal_dependency_id)          :: id_bedstress
      type (type_horizontal_dependency_id)          :: id_porosity
      type (type_horizontal_diagnostic_variable_id) :: id_v_er

      real(rk) :: v_cr
      real(rk) :: M
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self,configunit)
      class (type_ersem_benthic_erosion),intent(inout),target :: self
      integer,                           intent(in)           :: configunit

      call self%get_parameter(self%M,   'M',   'g/s/m4','erosion constant (Puls & Sundermann 1990)',       default=100.0_rk)
      call self%get_parameter(self%v_cr,'v_cr','m/s',   'critical bed shear velocity for sediment erosion',default=0.02_rk)

      call self%register_dependency(self%id_dens,     standard_variables%density)
      call self%register_dependency(self%id_bedstress,standard_variables%bottom_stress)
      call self%register_dependency(self%id_porosity, sediment_porosity)
      call self%register_diagnostic_variable(self%id_v_er,'v_er','m/d','erosion rate',standard_variable=sediment_erosion,source=source_do_bottom)
   end subroutine initialize

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_benthic_erosion),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: tau_bot, rho_wat, rho_sed, er, v_er, porosity
      real(rk),parameter :: rho_grain = 2650._rk  ! Sediment density (kg m-3) - single grain only, not the water or air around them! Quartz: 2650 kg/m3

      _HORIZONTAL_LOOP_BEGIN_
         _GET_(self%id_dens,rho_wat)
         _GET_HORIZONTAL_(self%id_bedstress,tau_bot)
         _GET_HORIZONTAL_(self%id_porosity,porosity)
         rho_sed = (1-porosity)*rho_grain                       ! Sediment density (dry sediment per total volume) - for erosion, it should be defined at the sediment surface.
         er = self%M*max(0.0_rk,tau_bot/rho_wat - self%v_cr**2) ! Erosion rate (g s-1 m-2) for sediment. Note: square of bed shear velocity = bed stress (Pa)/density (kg m-3)
         v_er = er/(1000*rho_sed)*86400                         ! Erosion in m d-1 (note: rho is expressed in kg m-3, hence the multiplication by 1000)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_v_er,v_er)
      _HORIZONTAL_LOOP_END_
   end subroutine do_bottom

end module
