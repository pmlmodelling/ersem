#include "fabm_driver.h"

module pml_ersem_pom

   use fabm_types
   use pml_ersem_base

   implicit none

!  default: all is private.
   private

   type,extends(type_ersem_pelagic_base_model),public :: type_pml_ersem_pom
      real(rk) :: EPSR6X
   contains
      procedure :: initialize
      procedure :: get_light_extinction
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
   real(rk) :: rRPmX,c0
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(has_s,'has_s',default=.true.)
      call self%get_parameter(has_f,'has_f',default=.true.)
      call self%get_parameter(rRPmX,'rRPmX',default=0.0_rk)
      call self%get_parameter(self%EPSR6X,'EPSR6X')
      call self%get_parameter(c0,'c0',default=0.0_rk)

      if (has_s.and.has_f) then
         call self%initialize_ersem_base(c_ini=0._rk,p_ini=0._rk,n_ini=0._rk,s_ini=0._rk,f_ini=0._rk,w=-rRPmX,sedimentation=.true.,c0=c0,n0=c0*qnRPIcX,p0=c0*qpRPIcX,s0=c0*qsRPIcX)
      elseif (has_s) then
         call self%initialize_ersem_base(c_ini=0._rk,p_ini=0._rk,n_ini=0._rk,s_ini=0._rk,w=-rRPmX,sedimentation=.true.,c0=c0,n0=c0*qnRPIcX,p0=c0*qpRPIcX,s0=c0*qsRPIcX)
      elseif (has_f) then
         call self%initialize_ersem_base(c_ini=0._rk,p_ini=0._rk,n_ini=0._rk,f_ini=0._rk,w=-rRPmX,sedimentation=.true.,c0=c0,n0=c0*qnRPIcX,p0=c0*qpRPIcX)
      else
         call self%initialize_ersem_base(c_ini=0._rk,p_ini=0._rk,n_ini=0._rk,w=-rRPmX,sedimentation=.true.,c0=c0,n0=c0*qnrpicX,p0=c0*qpRPIcX)
      end if
   end subroutine

   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
      class (type_pml_ersem_pom),intent(in) :: self
      _DECLARE_ARGUMENTS_GET_EXTINCTION_

      real(rk) :: c

      _LOOP_BEGIN_
         _GET_WITH_BACKGROUND_(self%id_c,c)
         _SET_EXTINCTION_(self%EPSR6X*c)
      _LOOP_END_
   end subroutine

end module