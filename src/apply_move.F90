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
    type (type_horizontal_diagnostic_variable_id)      :: id_target0_out
    type (type_horizontal_dependency_id)            :: id_distributed0
    type (type_horizontal_dependency_id)            :: id_thickness0
    type (type_bottom_dependency_id)   :: id_integral
    type (type_bottom_dependency_id)   :: id_totalthk
    type (type_dependency_id)                 :: id_thickness
    type (type_bottom_dependency_id)          :: id_topo
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
        call self%register_dependency(self%id_totalthk,'totalthk','','totalthk')
        call self%register_dependency(self%id_thickness, standard_variables%cell_thickness)
        call self%register_dependency(self%id_thickness0,vertical_integral(self%id_thickness))
        call self%register_dependency(self%id_topo,standard_variables%bottom_depth )
        call self%register_diagnostic_variable(self%id_target0_out,'migrator_integrated_mass','-','migrators final integrated mass', missing_value=0.0_rk, source=source_check_state)

    end subroutine initialize

    subroutine check_state(self,_ARGUMENTS_CHECK_STATE_)
        class (type_apply_move), intent(in) :: self
        _DECLARE_ARGUMENTS_CHECK_STATE_

        real(rk) :: local, distributed, distributed0, integral, difference,thickness0,thickness,topo,totalthk,target0

        _LOOP_BEGIN_
            _GET_(self%id_target,local)
            _GET_(self%id_distributed,distributed)
            _GET_(self%id_thickness,thickness)
            _GET_HORIZONTAL_(self%id_distributed0,distributed0)
            _GET_HORIZONTAL_(self%id_thickness0,thickness0)
            _GET_HORIZONTAL_(self%id_target0,target0)
            _GET_BOTTOM_(self%id_integral,integral)
            _GET_BOTTOM_(self%id_totalthk,totalthk)
            _GET_BOTTOM_(self%id_topo,topo)
            difference = abs( integral - distributed0 )
            !if (difference > 0.0_rk) then
                !write(*,*)'INTEGRALS_ARE_NOT_EQUAL',difference,thickness0,totalthk,topo
            !    _SET_(self%id_target, (integral - distributed0)/totalthk + distributed ) 
            !else
            !if (distributed < 0.0_rk) write(*,*)'DISTRIBUTED_NEGATIVE', thickness, topo, local
                _SET_(self%id_target, distributed)
            !end if
                !difference = integral - distributed0
            !conc = difference / thickness0
            
!            write(*,*) distributed
!            if (.not.ieee_is_finite( distributed ) ) write(*,*)'DISTRIBUTED_NOT_FINITE',distributed
!            _SET_(self%id_target, max(0.0_rk, distributed + (difference / thickness0) ) )
            !if (abs(difference)/integral < 0.0001_rk) then
!            if (topo >= 50.0_rk) then
!                _SET_(self%id_target, distributed)
!            end if
            !end if
            _LOOP_END_

    end subroutine check_state

end module    