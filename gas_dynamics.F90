#include "fabm_driver.h"
module pml_ersem_gas_dynamics

   use fabm_types
   use fabm_driver

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public pml_ersem_gas_dynamics_create, type_pml_ersem_gas_dynamics
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman

   type,extends(type_base_model) :: type_pml_ersem_gas_dynamics
!     Variable identifiers
      type (type_state_variable_id)      :: id_O3c,id_O2o

   end type
      
contains
      
   function pml_ersem_gas_dynamics_create(configunit,name,parent) result(self)
!
! !DESCRIPTION:
!  Here, the npzd namelist is read and the variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   integer,                          intent(in)    :: configunit
   character(len=*),                 intent(in)    :: name
   class (type_base_model),target,   intent(inout) :: parent
   class (type_pml_ersem_gas_dynamics), pointer       :: self
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   allocate(self)
   call self%initialize(name,parent)

   call self%register_state_variable(self%id_O3c,'O3c','mmol C/m^3','Carbon Dioxide', 2200._rk,minimum=0._rk)
   call self%register_state_variable(self%id_O2o,'O2o','mmol O/m^3','Oxygen',          300._rk,minimum=0._rk)

   end function pml_ersem_gas_dynamics_create

end module