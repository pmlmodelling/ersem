#include "fabm_driver.h"

module pml_ersem_nutrient

   use fabm_types

   implicit none

!  default: all is private.
   private

   type,extends(type_base_model),public :: type_pml_ersem_nutrient
      type (type_state_variable_id) :: id_N
   contains
      procedure :: initialize
   end type
      
contains
      
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_pml_ersem_nutrient), intent(inout), target :: self
   integer,                         intent(in)            :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
   character(len=1024) :: name,units
!EOP
!-----------------------------------------------------------------------
!BOC
   call self%get_parameter(name,'name')
   call self%get_parameter(units,'units')
   call self%register_state_variable(self%id_N,name,units,trim(self%long_name_prefix),1._rk,minimum=0._rk)

   end subroutine

end module