#include "fabm_driver.h"

module pml_ersem_nutrient

   use fabm_types
   use pml_ersem_base

   implicit none

!  default: all is private.
   private

   type,extends(type_ersem_base_model),public :: type_pml_ersem_nutrient
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
   select case (name)
      case ('n')
         call self%initialize_ersem_base(n_ini=0._rk)
      case ('p')
         call self%initialize_ersem_base(p_ini=0._rk)
      case ('f')
         call self%initialize_ersem_base(f_ini=0._rk)
      case ('s')
         call self%initialize_ersem_base(s_ini=0._rk)
      case default
         call self%fatal_error('initialize','unknown element '//trim(name))
   end select

   end subroutine

end module