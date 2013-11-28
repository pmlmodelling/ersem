#include "fabm_driver.h"
module pml_ersem_gas_dynamics

   use fabm_types

   implicit none

   private
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman

   type,extends(type_base_model),public :: type_pml_ersem_gas_dynamics
!     Variable identifiers
      type (type_state_variable_id)     :: id_O3c,id_O2o
      type (type_conserved_quantity_id) :: id_totc

      type (type_diagnostic_variable_id) :: id_eO2mO2
   contains
      procedure :: initialize
   end type
      
contains
      
   subroutine initialize(self,configunit)
!
! !INPUT PARAMETERS:
      class (type_pml_ersem_gas_dynamics), intent(inout), target :: self
      integer,                             intent(in)            :: configunit
!
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%register_state_variable(self%id_O3c,'O3c','mmol C/m^3','Carbon Dioxide', 2200._rk,minimum=0._rk)
      call self%register_state_variable(self%id_O2o,'O2o','mmol O/m^3','Oxygen',          300._rk,minimum=0._rk)

      call self%register_conserved_quantity(self%id_totc,standard_variables%total_carbon)
      call self%add_conserved_quantity_component(self%id_totc,self%id_O3c)
   
      call self%register_diagnostic_variable(self%id_eO2mO2,'eO2mO2','1','oxygen saturation', &
         missing_value=1._rk,standard_variable=standard_variables%fractional_saturation_of_oxygen)

   end subroutine

end module