#include "fabm_driver.h"
#define IRON
module pml_ersem_benthic_pom

   use fabm_types
   use pml_ersem_base

   implicit none

!  default: all is private.
   private

   type,extends(type_ersem_benthic_base_model),public :: type_pml_ersem_benthic_pom
      type (type_state_variable_id) :: id_resuspenion_c,id_resuspenion_n,id_resuspenion_p,id_resuspenion_s,id_resuspenion_f
   contains
      procedure :: initialize
   end type

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_pml_ersem_benthic_pom), intent(inout), target :: self
   integer,                            intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   logical :: has_s,has_f
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(has_s,'has_s',default=.true.)
      call self%get_parameter(has_f,'has_f',default=.true.)

      if (has_s.and.has_f) then
         call self%initialize_ersem_benthic_base(c_ini=0._rk,p_ini=0._rk,n_ini=0._rk,s_ini=0._rk,f_ini=0._rk)
      elseif (has_s) then
         call self%initialize_ersem_benthic_base(c_ini=0._rk,p_ini=0._rk,n_ini=0._rk,s_ini=0._rk)
      elseif (has_f) then
         call self%initialize_ersem_benthic_base(c_ini=0._rk,p_ini=0._rk,n_ini=0._rk,f_ini=0._rk)
      else
         call self%initialize_ersem_benthic_base(c_ini=0._rk,p_ini=0._rk,n_ini=0._rk)
      end if

      call self%register_state_dependency(self%id_resuspenion_c,'resuspension_target_c','mmol m-3','pelagic variable taking up resuspended carbon',    required=.false.)
      call self%register_state_dependency(self%id_resuspenion_n,'resuspension_target_n','mmol m-3','pelagic variable taking up resuspended nitrogen',  required=.false.)
      call self%register_state_dependency(self%id_resuspenion_p,'resuspension_target_p','mmol m-3','pelagic variable taking up resuspended phosphorus',required=.false.)
      if (has_s) call self%register_state_dependency(self%id_resuspenion_s,'resuspension_target_s','mmol m-3','pelagic variable taking up resuspended silicate',required=.false.)
      if (has_f) call self%register_state_dependency(self%id_resuspenion_f,'resuspension_target_f','umol m-3','pelagic variable taking up resuspended iron',required=.false.)
   end subroutine

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_pml_ersem_benthic_pom), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
   end subroutine
end module