#include "fabm_driver.h"
module pml_ersem_dom

   use fabm_types

   implicit none

!  default: all is private.
   private

   type,extends(type_base_model),public :: type_pml_ersem_dom
      type (type_state_variable_id)      :: id_c,id_p,id_n
   contains
      procedure                          :: initialize
   end type
      
contains
      
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_pml_ersem_dom), intent(inout),target :: self
   integer,                    intent(in)           :: configunit
!
! !LOCAL VARIABLES:
   logical :: has_np
!EOP
!-----------------------------------------------------------------------
!BOC
   call self%get_parameter(has_np,'has_np',default=.true.)
   call self%register_state_variable(self%id_c,'c','mg C/m^3','carbon',0._rk,minimum=0._rk)
   if (has_np) then
      call self%register_state_variable(self%id_p,'p','mmol P/m^3','phosphorus',0._rk,minimum=0._rk)
      call self%register_state_variable(self%id_n,'n','mmol N/m^3','nitrogen',  0._rk,minimum=0._rk)
   end if

   end subroutine

end module