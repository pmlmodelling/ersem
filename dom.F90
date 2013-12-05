#include "fabm_driver.h"
module pml_ersem_dom

   use fabm_types
   use pml_ersem_base

   implicit none

!  default: all is private.
   private

   type,extends(type_ersem_pelagic_base_model),public :: type_pml_ersem_dom
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
   if (has_np) then
      call self%initialize_ersem_base(c_ini=0._rk,p_ini=0._rk,n_ini=0._rk)
   else
      call self%initialize_ersem_base(c_ini=0._rk)
   end if

   end subroutine

end module