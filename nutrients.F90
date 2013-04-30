#include "fabm_driver.h"
#define IRON
module pml_ersem_nutrients

   use fabm_types
   use fabm_driver

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public pml_ersem_nutrients_create, type_pml_ersem_nutrients
!
! !PRIVATE DATA MEMBERS:
   real(rk), parameter :: secs_pr_day = 86400.0_rk
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman

   type,extends(type_base_model) :: type_pml_ersem_nutrients
!     Variable identifiers
      type (type_state_variable_id)      :: id_N1p,id_N3n,id_N4n,id_N5s
#ifdef IRON   
      type (type_state_variable_id)      :: id_N7f
#endif

      end type
      
contains
      
   function pml_ersem_nutrients_create(configunit,name,parent) result(self)
!
! !DESCRIPTION:
!  Here, the npzd namelist is read and the variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   integer,                          intent(in)    :: configunit
   character(len=*),                 intent(in)    :: name
   class (type_base_model),target,   intent(inout) :: parent
   class (type_pml_ersem_nutrients), pointer       :: self
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   allocate(self)
   call self%initialize(name,parent)

   call self%register_state_variable(self%id_N1p,'N1p','mmol P/m^3','Phosphate', 1._rk,minimum=0._rk)
   call self%register_state_variable(self%id_N3n,'N3n','mmol N/m^3','Nitrate',   1._rk,minimum=0._rk)
   call self%register_state_variable(self%id_N4n,'N4n','mmol N/m^3','Ammonium',  1._rk,minimum=0._rk)
   call self%register_state_variable(self%id_N5s,'N5s','mmol S/m^3','Silicate',  1._rk,minimum=0._rk)
#ifdef IRON   
   call self%register_state_variable(self%id_N7f,'N7f','umol F/m^3','Inorganic Iron', 1._rk, minimum=0._rk)
#endif

   end function pml_ersem_nutrients_create

end module