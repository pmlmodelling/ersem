#include "fabm_driver.h"

module dvm_apply_move

use fabm_types
use fabm_expressions
use, intrinsic :: ieee_arithmetic

implicit none 

private 

public type_apply_move

type, extends(type_base_model) :: type_apply_move

    type (type_dependency_id)          :: id_distributed
    type (type_state_variable_id)      :: id_target
    type (type_horizontal_dependency_id)      :: id_target0
    type (type_horizontal_dependency_id)            :: id_distributed0
    type (type_bottom_dependency_id)   :: id_integral
    contains
        procedure :: initialize
        procedure :: check_state
    end type

contains

    subroutine initialize(self, configunit)

        class (type_apply_move), intent(inout), target :: self
        integer,         intent(in)            :: configunit

        call self%register_state_dependency(self%id_target, 'target', '', 'variable to apply sources and sinks to')
        call self%register_dependency(self%id_target0,vertical_integral(self%id_target))
        call self%register_dependency(self%id_distributed,'migrator_distributed_mass','-','migrators final mass distribution')
        call self%register_dependency(self%id_distributed0,vertical_integral(self%id_distributed))
        call self%register_dependency(self%id_integral,'integral','','depth-integrated target variable')

    end subroutine initialize

    subroutine check_state(self,_ARGUMENTS_CHECK_STATE_)
        class (type_apply_move), intent(in) :: self
        _DECLARE_ARGUMENTS_CHECK_STATE_

        real(rk) :: local, distributed, distributed0, integral, difference,target0

        _LOOP_BEGIN_
            _GET_(self%id_target,local)
            _GET_(self%id_distributed,distributed)
            _GET_HORIZONTAL_(self%id_distributed0,distributed0)
            _GET_BOTTOM_(self%id_integral,integral)
            
            difference = abs( integral - distributed0 )
            _SET_(self%id_target, distributed)

        _LOOP_END_

    end subroutine check_state

end module    
