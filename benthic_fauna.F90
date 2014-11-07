#include "fabm_driver.h"

module ersem_benthic_fauna

   use fabm_types

   use fabm_particle
   use ersem_shared
   use ersem_benthic_base

   implicit none

   private

  type,extends(type_ersem_benthic_base),public ::
     type_ersem_benthic_fauna


  contains
     procedure :: initialize
     procedure :: do_bottom
  end type

contains

  subroutine initialize(self,configunit)
    class (type_ersem_benthic_fauna),intent(inout),target :: self
    integer,                                 intent(in)           ::configunit
  end subroutine

  subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)

     class (type_ersem_benthic_fauna),intent(in) :: self
     _DECLARE_ARGUMENTS_DO_BOTTOM_
     
     _HORIZONTAL_LOOP_BEGIN_
      _HORIZONTAL_LOOP_END_

  end subroutine do_bottom

end module

