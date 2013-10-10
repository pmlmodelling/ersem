#include "fabm_driver.h"
#define IRON
module pml_ersem_pom

   use fabm_types
   use fabm_driver

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_pml_ersem_pom
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman

   type,extends(type_base_model) :: type_pml_ersem_pom
!     Variable identifiers
      type (type_state_variable_id)      :: id_R6c,id_R6p,id_R6n,id_R6s
#ifdef IRON
      type (type_state_variable_id)      :: id_R6f
#endif
   contains

      procedure :: initialize

   end type

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the npzd namelist is read and the variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_pml_ersem_pom), intent(inout), target :: self
   integer,                    intent(in)            :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   call self%register_state_variable(self%id_R6c,'R6c','mg C/m^3',  'POC', 0._rk,minimum=0._rk)
   call self%register_state_variable(self%id_R6p,'R6p','mmol P/m^3','POP', 0._rk,minimum=0._rk)
   call self%register_state_variable(self%id_R6n,'R6n','mmol N/m^3','PON', 0._rk,minimum=0._rk)
   call self%register_state_variable(self%id_R6s,'R6s','mmol N/m^3','POS', 0._rk,minimum=0._rk)
#ifdef IRON   
   call self%register_state_variable(self%id_R6f,'R6f','umol N/m^3','POF', 0._rk,minimum=0._rk)
#endif
   end subroutine

end module