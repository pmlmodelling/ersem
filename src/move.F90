#include "fabm_driver.h"

module dvm_move

use fabm_types
use fabm_expressions

implicit none 

private 

public type_move

type, extends(type_base_model) :: type_move
    type (type_dependency_id)          :: id_random_weights
    type (type_dependency_id)          :: id_thickness
    type (type_state_variable_id)      :: id_target
    type (type_horizontal_dependency_id)      :: id_target0
    !type (type_bottom_dependency_id)   :: id_integral
    type (type_bottom_dependency_id)   :: id_integral_random_weights
    type (type_diagnostic_variable_id) :: id_distributed
    real(rk) :: m, ratioMig
    contains
        procedure :: initialize
        procedure :: check_state
    end type

contains

    subroutine initialize(self, configunit)

        class (type_move), intent(inout), target :: self
        integer,         intent(in)            :: configunit

        call self%register_state_dependency(self%id_target, 'target', '', 'variable to apply sources and sinks to')
        call self%register_dependency(self%id_target0,vertical_integral(self%id_target))
        call self%register_dependency(self%id_random_weights,'migrator_random_weights','-','migrators distribution random weights')
        !call self%register_dependency(self%id_integral,'integral','','depth-integrated target variable')
        call self%register_dependency(self%id_integral_random_weights,'migrator_integral_random_weights','','migrators distribution integral random weights')
        call self%register_dependency(self%id_thickness, standard_variables%cell_thickness)
        call self%get_parameter( self%ratioMig,  'ratioMig',    '-',      'ratio of moving biomass', default=0.5_rk)

    end subroutine initialize

   subroutine check_state(self,_ARGUMENTS_CHECK_STATE_)
       class (type_move), intent(in) :: self
       _DECLARE_ARGUMENTS_CHECK_STATE_

        real(rk) :: integral_random_weights, random_weights !integral
        real(rk) :: local, distributed, thickness, target0
        !_GET_BOTTOM_(self%id_integral,integral)
        _GET_BOTTOM_(self%id_integral_random_weights,integral_random_weights)
        _GET_HORIZONTAL_(self%id_target0,target0)
        
        _LOOP_BEGIN_
           _GET_(self%id_target,local)
           _GET_(self%id_random_weights,random_weights)
           _GET_(self%id_thickness,thickness)
    
           distributed = random_weights / integral_random_weights ! this ensures that total weight distribution is 1, so mass conserved
           _SET_(self%id_target, max(0.0_rk,local * (1.0_rk - self%ratioMig) + (distributed * target0/max(thickness,1.0E-20_rk)) * self%ratioMig ) )
    
        _LOOP_END_
     end subroutine check_state

end module
