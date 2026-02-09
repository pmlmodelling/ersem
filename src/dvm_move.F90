#include "fabm_driver.h"

! -----------------------------------------------------------------------------
! For diel vertical migration of migrating mesozooplankton. Moves migrator 
! vertically to new location in water column based on distribution calculated 
! in dvm_weight_distribution 
!
! Migration only occurs when mesozooplankton are not overwintering.
!
! Adapted from code written by Caglar Yumruktepe (NERSC), available at:
! https://github.com/nansencenter/nersc
! -----------------------------------------------------------------------------

module dvm_move

use fabm_types
use fabm_expressions

implicit none 

private 

public type_move

type, extends(type_base_model) :: type_move
    type (type_dependency_id)          :: id_weights
    type (type_dependency_id)          :: id_thickness
    type (type_state_variable_id)      :: id_target
    type (type_bottom_dependency_id)   :: id_integral
    type (type_bottom_dependency_id)   :: id_integral_weights
    type (type_bottom_dependency_id)   :: id_overwintering
    real(rk) :: tstep, ratioMig
    contains
        procedure :: initialize
        procedure :: do
    end type

contains

    subroutine initialize(self, configunit)

        class (type_move), intent(inout), target :: self
        integer,         intent(in)            :: configunit

        call self%register_state_dependency(self%id_target, 'target', '', 'variable to apply sources and sinks to')
        call self%register_dependency(self%id_weights,'migrator_weights','-','migrators distribution weights')
        call self%register_dependency(self%id_integral,'integral','','depth-integrated target variable')
        call self%register_dependency(self%id_integral_weights,'migrator_integral_weights','','migrators distribution integral weights')
        call self%register_dependency(self%id_thickness, standard_variables%cell_thickness)
        call self%register_dependency(self%id_overwintering, 'overwintering','','Is mesozooplankton in diapause?')
        call self%get_parameter(self%ratioMig, 'ratioMig', '-', 'ratio of moving biomass', default=0.5_rk)
        call self%get_parameter(self%tstep, 'tstep', 'sec', 'time-step in seconds', default=900.0_rk)

    end subroutine initialize

    subroutine do(self,_ARGUMENTS_DO_)
       class (type_move), intent(in) :: self
       _DECLARE_ARGUMENTS_DO_

        real(rk) :: integral_weights, weights, integral
        real(rk) :: local, distributed, thickness, overwintering
        
        _LOOP_BEGIN_
           _GET_BOTTOM_(self%id_integral,integral)
           _GET_BOTTOM_(self%id_integral_weights,integral_weights)         
           _GET_BOTTOM_(self%id_overwintering,overwintering)
           _GET_(self%id_target,local)
           _GET_(self%id_weights,weights)
           _GET_(self%id_thickness,thickness)

           if (overwintering < 1.0_rk) then
              ! Make sure total weight distribution is 1, so mass conserved
              distributed = weights / integral_weights               
              ! Move migrator to new location in watercolumn
              _ADD_SOURCE_(self%id_target, (distributed * integral - local) * self%ratioMig / self%tstep)
              !_ADD_SOURCE_(self%id_target,  max(-max(local,0.0_rk),(distributed * integral /max(thickness,1.0E-20_rk) - max(local,0.0_rk)) * self%ratioMig) / self%tstep)
           end if
        _LOOP_END_
        
    end subroutine do

end module
