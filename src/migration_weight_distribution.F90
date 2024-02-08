#include "fabm_driver.h"

module ecosmo_migration_weight_distribution

use fabm_types
use fabm_expressions
!use ecosmo_shared

private 

type, extends(type_base_model), public :: type_weights
    type (type_bottom_diagnostic_variable_id) :: id_integral
    type (type_bottom_diagnostic_variable_id) :: id_integral_random_weights
    type (type_diagnostic_variable_id)        :: id_random_weights
    type (type_state_variable_id)             :: id_target
    type (type_dependency_id)                 :: id_present
    type (type_dependency_id)                 :: id_thickness
    type (type_horizontal_dependency_id) :: id_par0, id_parmean0

    contains
        procedure :: initialize => depth_integral_initialize
        procedure :: do_column  => depth_integral_do_column
end type

contains

subroutine depth_integral_initialize(self, configunit)
    class (type_weights), intent(inout), target :: self
    integer,                     intent(in)            :: configunit
    call self%register_diagnostic_variable(self%id_integral_random_weights,'migrator_integral_random_weights','-','migrators distribution integral random weights', missing_value=0.0_rk, &
        act_as_state_variable=.true., source=source_do_column)
    call self%register_diagnostic_variable(self%id_random_weights,'migrator_random_weights','-','migrators distribution random weights', missing_value=0.0_rk, &
        act_as_state_variable=.true., source=source_do_column)
    call self%register_state_dependency(self%id_target, 'target', '', 'variable to depth-integrate')   
    call self%register_dependency(self%id_present,'present', '-', 'migrators are present here')
    call self%register_diagnostic_variable(self%id_integral, 'integral', '','integral', missing_value=0.0_rk, &
        act_as_state_variable=.true., source=source_do_column)
    call self%register_dependency(self%id_thickness, standard_variables%cell_thickness)
    call self%register_dependency(self%id_par0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)
    call self%register_dependency(self%id_parmean0,temporal_mean(self%id_par0,period=86400._rk,resolution=3600._rk))

end subroutine depth_integral_initialize

subroutine depth_integral_do_column(self, _ARGUMENTS_DO_COLUMN_)
    class (type_weights), intent(in) :: self
    _DECLARE_ARGUMENTS_DO_COLUMN_

    real(rk) :: local, weight, thickness, integral, integral_random, present, local_random, minimum_random, minimum_value
    real(rk) :: totaldepth
    real(rk) :: parmean0

    integral = 0.0_rk
    integral_random = 0.0_rk
    totaldepth = 0.0_rk
    _GET_SURFACE_(self%id_parmean0,parmean0)
    call random_seed()

    _VERTICAL_LOOP_BEGIN_

       _GET_(self%id_target,local)
       _GET_(self%id_present,present)
       _GET_(self%id_thickness,thickness)
       integral = integral + local*thickness
       totaldepth = totaldepth + thickness

       call random_number(local_random)
       call random_number(minimum_random)
       !local_random = 0.1_rk + 0.9_rk * local_random * present
       !local_random = local_random * present
       minimum_value = 0.02 + (0.1 - 0.02) * minimum_random
       !write(*,*)parmean0
       if (parmean0 < 1E-4_rk) then
            if (totaldepth <= 400.0_rk) then
                local_random = thickness * (minimum_value + (1.0_rk - minimum_value) * local_random)
            else
                local_random = thickness * minimum_value
            end if 
       else
            local_random = thickness * (minimum_value + (1.0_rk - minimum_value) * local_random * present)
       end if
       !write(*,*)local_random
       integral_random = integral_random + local_random

       _SET_DIAGNOSTIC_(self%id_random_weights,local_random)

    _VERTICAL_LOOP_END_
    _SET_BOTTOM_DIAGNOSTIC_(self%id_integral,integral)
    _SET_BOTTOM_DIAGNOSTIC_(self%id_integral_random_weights,integral_random)

end subroutine depth_integral_do_column

end module
