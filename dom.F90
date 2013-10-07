#include "fabm_driver.h"
module pml_ersem_dom

   use fabm_types
   use fabm_driver

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_pml_ersem_dom
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman

   type,extends(type_base_model) :: type_pml_ersem_dom
!     Variable identifiers
      type (type_state_variable_id)      :: id_R1c,id_R1p,id_R1n
      type (type_state_variable_id)      :: id_R2c

      end type
      
contains
      
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the npzd namelist is read and the variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_pml_ersem_dom), intent(inout),target :: self
   integer,                    intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   call self%register_state_variable(self%id_R1c,'R1c','mg C/m^3',  'DOC',            0._rk,minimum=0._rk)
   call self%register_state_variable(self%id_R1p,'R1p','mmol P/m^3','DOP',            0._rk,minimum=0._rk)
   call self%register_state_variable(self%id_R1n,'R1n','mmol N/m^3','DON',            0._rk,minimum=0._rk)
   call self%register_state_variable(self%id_R2c,'R2c','mg C/m^3',  'Semi-labile DOC',0._rk,minimum=0._rk)

   end subroutine

end module