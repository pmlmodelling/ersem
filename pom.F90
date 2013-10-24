#include "fabm_driver.h"
#define IRON
module pml_ersem_pom

   use fabm_types

   implicit none

!  default: all is private.
   private

   type,extends(type_base_model),public :: type_pml_ersem_pom
      type (type_state_variable_id)      :: id_c,id_p,id_n,id_s
#ifdef IRON
      type (type_state_variable_id)      :: id_f
#endif
   contains
      procedure :: initialize
   end type

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_pml_ersem_pom), intent(inout), target :: self
   integer,                    intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   logical :: has_s,has_f
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(has_s,'has_s',default=.true.)
      call self%get_parameter(has_f,'has_f',default=.true.)

      call self%register_state_variable(self%id_c,'c','mg C/m^3',  'carbon',     0._rk,minimum=0._rk)
      call self%register_state_variable(self%id_p,'p','mmol P/m^3','phosphorous',0._rk,minimum=0._rk)
      call self%register_state_variable(self%id_n,'n','mmol N/m^3','nitrogen',   0._rk,minimum=0._rk)

      if (has_s) call self%register_state_variable(self%id_s,'s','mmol Si/m^3','silicate',0._rk,minimum=0._rk)
      if (has_f) call self%register_state_variable(self%id_f,'f','umol Fe/m^3','iron',    0._rk,minimum=0._rk)
   end subroutine

end module