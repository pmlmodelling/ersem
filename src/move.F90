#include "fabm_driver.h"

module dvm_move

use fabm_types
use fabm_expressions
!use, intrinsic :: ieee_arithmetic

implicit none 

private 

public type_move

type, extends(type_base_model) :: type_move
!type, extends(type_base_model), public :: type_move
    type (type_dependency_id)          :: id_random_weights
    type (type_dependency_id)          :: id_thickness
    type (type_state_variable_id)      :: id_target
    type (type_horizontal_dependency_id)      :: id_target0
    type (type_bottom_dependency_id)   :: id_integral
    type (type_bottom_dependency_id)   :: id_integral_random_weights
!    type (type_bottom_diagnostic_variable_id) :: id_integral_random_weights_fabm
    type (type_diagnostic_variable_id) :: id_distributed
!    type (type_bottom_diagnostic_variable_id)      :: id_target0_out
    real(rk) :: tstep, m, ratioMig
    contains
        procedure :: initialize
!        procedure :: do
        procedure :: check_state
    end type

contains

    subroutine initialize(self, configunit)

        class (type_move), intent(inout), target :: self
        integer,         intent(in)            :: configunit

        call self%register_state_dependency(self%id_target, 'target', '', 'variable to apply sources and sinks to')
        call self%register_dependency(self%id_target0,vertical_integral(self%id_target))
        call self%register_dependency(self%id_random_weights,'migrator_random_weights','-','migrators distribution random weights')
        call self%register_dependency(self%id_integral,'integral','','depth-integrated target variable')
        call self%register_dependency(self%id_integral_random_weights,'migrator_integral_random_weights','','migrators distribution integral random weights')
!        call self%register_diagnostic_variable(self%id_distributed,'migrator_distributed_mass','-','migrators final mass distribution', missing_value=0.0_rk, source=source_do)
        call self%register_dependency(self%id_thickness, standard_variables%cell_thickness)
        call self%get_parameter( self%tstep,  'tstep',    'sec',      'time-step in seconds', default=1200.0_rk)
        call self%get_parameter( self%ratioMig,  'ratioMig',    '-',      'ratio of moving biomass', default=0.5_rk)
!        call self%register_diagnostic_variable(self%id_target0_out,'migrator_integrated_mass','-','migrators final integrated mass', missing_value=0.0_rk, source=source_check_state)
!        call self%register_dependency(self%id_integral_random_weights_fabm,vertical_integral(self%id_random_weights))

    end subroutine initialize

   subroutine check_state(self,_ARGUMENTS_CHECK_STATE_)
       class (type_move), intent(in) :: self
       _DECLARE_ARGUMENTS_CHECK_STATE_

    ! subroutine do(self, _ARGUMENTS_DO_)
    !     class (type_move), intent(in) :: self
    !     _DECLARE_ARGUMENTS_DO_
    
        real(rk) :: integral, integral_random_weights, integral_random_weights_fabm, random_weights
        real(rk) :: local, distributed, thickness, target0
        real(rk) :: mortality_switch, local_loss
        !real(rk) :: final_concentration
    
!        _GET_HORIZONTAL_(self%id_integral,integral)
!        _GET_HORIZONTAL_(self%id_integral_random_weights,integral_random_weights)
!        _GET_HORIZONTAL_(self%id_integral_random_weights_fabm,integral_random_weights_fabm)
        !write(*,*)integral
!        _GET_BOTTOM_(self%id_integral,integral)
!        _GET_BOTTOM_(self%id_integral_random_weights,integral_random_weights)
        _LOOP_BEGIN_
           _GET_BOTTOM_(self%id_integral,integral)
           _GET_BOTTOM_(self%id_integral_random_weights,integral_random_weights)
!           _GET_BOTTOM_(self%id_integral_random_weights_fabm,integral_random_weights_fabm)
           _GET_(self%id_target,local)
           _GET_HORIZONTAL_(self%id_target0,target0)
           _GET_(self%id_random_weights,random_weights)
           _GET_(self%id_thickness,thickness)
    
           distributed = random_weights / integral_random_weights ! this ensures that total weight distribution is 1, so mass conserved
           !if (distributed < 0.0_rk) then
           !   write(*,*) 'DVM_DISTRIBUTED',distributed,random_weights,integral_random_weights
           !end if
!           if (integral_random_weights<0.0_rk) write(*,'(A8,1x,E12.4,1x,E12.4,1x,E12.4,1x,E12.4)') 'INT_MOVE', integral_random_weights, integral_random_weights_fabm,random_weights,integral
           !if (integral_random_weights<0.0_rk) 
           !write(*,'(A8,1x,E12.4,1x,E12.4,1x,E12.4)') 'INT_MOVE', integral_random_weights,random_weights,integral
!           _SET_DIAGNOSTIC_(self%id_distributed,distributed)
           !write(*,*)par,mortality_switch,local_loss,local
!           if (local < 0.0_rk) write(*,*) 'LOCAL_BELOW_0.0', local
!           if (distributed < 0.0_rk) write(*,*) 'DITRIBUTED_BELOW_0.0', distributed
!           if (integral < 0.0_rk) write(*,*) 'INTEGRAL_BELOW_0.0', integral
!           if (distributed * integral * self%ratioMig / thickness - local < 0.0_rk) write(*,'(A20,1x,F14.6,1x,F14.6,1x,F14.6,1x,F14.6,1x,F14.6,1x,F14.6)') &
!                'CALCULATED_BELOW_0.0', distributed * integral * self%ratioMig / thickness - local, &
!                distributed, integral, self%ratioMig, thickness, local  
!           _ADD_SOURCE_(self%id_target,  max(-max(local,0.0_rk),(distributed * integral/thickness - max(local,0.0_rk) ) * self%ratioMig) / self%tstep) 
           !final_concentration = local * (1.0_rk - self%ratioMig) + (distributed * integral/thickness) * self%ratioMig
!           write(*,*) (local * (1.0_rk - self%ratioMig) + (distributed * integral/thickness) * self%ratioMig) &
!                     , local, distributed, integral
!           if (.not.ieee_is_finite( local ) ) write(*,*)'LOCAL_NOT_FINITE',local
!           if (.not.ieee_is_finite( distributed ) ) write(*,*)'DISTRIBUTED_NOT_FINITE',distributed
!           if (.not.ieee_is_finite( self%ratioMig ) ) write(*,*)'RATIO_NOT_FINITE',self%ratioMig
!           if (.not.ieee_is_finite( integral ) ) write(*,*)'INTEGRAL_NOT_FINITE',integral
!           if (.not.ieee_is_finite( thickness ) ) write(*,*)'THICKNESS_NOT_FINITE',thickness
!           if (.not.ieee_is_finite( integral/thickness ) ) write(*,*)'INToverTHK_NOT_FINITE',integral/thickness
!        if ( ieee_is_finite(final_concentration) ) then
!           if (thickness>1.0_rk) then 
           !_SET_DIAGNOSTIC_(self%id_distributed,local * (1.0_rk - self%ratioMig) + (distributed * integral/thickness) * self%ratioMig)
! BU BU          _ADD_SOURCE_(self%id_target,  max(-max(local,0.0_rk),(distributed * integral/thickness - max(local,0.0_rk) ) * self%ratioMig) / self%tstep) 
!           else
!            _SET_DIAGNOSTIC_(self%id_distributed,local)
!           end if
           !_SET_(self%id_target, final_concentration )
!           _SET_(self%id_target, local )
!        else
!            _SET_(self%id_target, local )
!        end if 
            !write(*,*) integral,target0
           ! 
!           _ADD_SOURCE_(self%id_target,  max(-max(local,0.0_rk),(distributed * target0/max(thickness,1.0E-20_rk) - max(local,0.0_rk) ) * self%ratioMig) / self%tstep)
           _SET_(self%id_target, max(0.0_rk,local * (1.0_rk - self%ratioMig) + (distributed * target0/max(thickness,1.0E-20_rk)) * self%ratioMig ) )
           !_SET_(self%id_target, local )
           !write(*,*)self%dt
    
        _LOOP_END_
!      _SET_BOTTOM_DIAGNOSTIC_(self%id_target0_out,target0)
!     end subroutine do
     end subroutine check_state

    !  subroutine check_state(self,_ARGUMENTS_CHECK_STATE_)
    !     class (type_move), intent(in) :: self
    !     _DECLARE_ARGUMENTS_CHECK_STATE_

    !     real(rk) :: local, distributed

    !     _LOOP_BEGIN_
    !         _GET_(self%id_target,local)
    !         _GET_(self%id_distributed,distributed)
    !         write(*,*) distributed
    !         !if (.not.ieee_is_finite( distributed ) ) write(*,*)'DISTRIBUTED_NOT_FINITE',distributed
    !         _SET_(self%id_target, local )
    !     _LOOP_END_

    ! end subroutine check_state

end module
