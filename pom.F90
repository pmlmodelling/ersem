#include "fabm_driver.h"
#define IRON
module pml_ersem_pom

   use fabm_types
   use pml_ersem_base

   implicit none

!  default: all is private.
   private

   type,extends(type_ersem_pelagic_base_model),public :: type_pml_ersem_pom
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
   real(rk) :: rRPmX
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(has_s,'has_s',default=.true.)
      call self%get_parameter(has_f,'has_f',default=.true.)
      call self%get_parameter(rRPmX,'rRPmX',default=0.0_rk)

      if (has_s.and.has_f) then
         call self%initialize_ersem_base(c_ini=0._rk,p_ini=0._rk,n_ini=0._rk,s_ini=0._rk,f_ini=0._rk,w=-rRPmX,sedimentation=.true.)
      elseif (has_s) then
         call self%initialize_ersem_base(c_ini=0._rk,p_ini=0._rk,n_ini=0._rk,s_ini=0._rk,w=-rRPmX,sedimentation=.true.)
      elseif (has_f) then
         call self%initialize_ersem_base(c_ini=0._rk,p_ini=0._rk,n_ini=0._rk,f_ini=0._rk,w=-rRPmX,sedimentation=.true.)
      else
         call self%initialize_ersem_base(c_ini=0._rk,p_ini=0._rk,n_ini=0._rk,w=-rRPmX,sedimentation=.true.)
      end if
   end subroutine

end module