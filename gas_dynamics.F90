#include "fabm_driver.h"
module pml_ersem_gas_dynamics

   use fabm_types
   use fabm_driver

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_pml_ersem_gas_dynamics
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman

   type,extends(type_base_model) :: type_pml_ersem_gas_dynamics
!     Variable identifiers
      type (type_state_variable_id)     :: id_O3c,id_O2o
      type (type_conserved_quantity_id) :: id_totc
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

   end subroutine

end module