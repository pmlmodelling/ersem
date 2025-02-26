#include "fabm_driver.h"

module dvm_conservative_migrator

use fabm_types
use fabm_expressions
!use ecosmo_shared

implicit none

private

type,extends(type_base_model), public  :: type_conservative_migrator
    type (type_state_variable_id)         :: id_c
    contains

!     Model procedures
    procedure :: initialize

end type type_conservative_migrator

contains
subroutine initialize(self,configunit)
    !
    ! !INPUT PARAMETERS:
    class (type_conservative_migrator), intent(inout),target  :: self
    integer,  intent(in) :: configunit
    !
    ! !REVISION HISTORY
    !
    !  Veli Çağlar Yumruktepe:
    !       XXX
    call self%register_state_variable(self%id_c, 'c', 'mgC/m3', 'concentration in carbon units', minimum=0.0_rk)

end subroutine initialize

end module dvm_conservative_migrator

