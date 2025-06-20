#include "fabm_driver.h"

module dvm_weight_distribution

use fabm_types
use fabm_expressions

implicit none

private 

type, extends(type_base_model), public :: type_weight_distribution
    type (type_bottom_diagnostic_variable_id) :: id_integral
    type (type_bottom_diagnostic_variable_id) :: id_integral_weights
    type (type_diagnostic_variable_id)        :: id_weights
    type (type_state_variable_id)             :: id_target
    type (type_dependency_id)                 :: id_present
    type (type_dependency_id)                 :: id_migrator_food
    type (type_dependency_id)                 :: id_thickness
    type (type_surface_dependency_id)         :: id_par0
    type (type_bottom_dependency_id)          :: id_topo
    type (type_dependency_id)                 :: id_depth    

    contains
        procedure :: initialize
        procedure :: do_column
end type

contains

    subroutine initialize(self, configunit)
        class (type_weight_distribution), intent(inout), target :: self
        integer,                     intent(in)            :: configunit
        call self%register_diagnostic_variable(self%id_integral_weights,'migrator_integral_weights','-','migrators distribution integral weights', missing_value=0.0_rk, &
            act_as_state_variable=.true., source=source_do_column)
        call self%register_diagnostic_variable(self%id_weights,'migrator_weights','-','migrators distribution weights', missing_value=0.0_rk, &
            act_as_state_variable=.true., source=source_do_column)
        call self%register_state_dependency(self%id_target, 'target', '', 'variable to depth-integrate')   
        call self%register_dependency(self%id_present,'present', '-', 'migrators are present here')
        call self%register_dependency(self%id_migrator_food,'migrator_food', 'mgC/m3', 'food availability for the migrators')
        !call self%register_dependency(self%id_migrator_food0,vertical_integral(self%id_migrator_food))
        call self%register_diagnostic_variable(self%id_integral, 'integral', '','integral', missing_value=0.0_rk, &
            act_as_state_variable=.true., source=source_do_column)
        !call self%register_diagnostic_variable(self%id_totalthk, 'totalthk', '','totalthk', missing_value=0.0_rk, &
        !    source=source_do_column)
        call self%register_dependency(self%id_thickness, standard_variables%cell_thickness)
        call self%register_dependency(self%id_par0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)
        call self%register_dependency(self%id_depth,standard_variables%depth)
        call self%register_dependency(self%id_topo,standard_variables%bottom_depth_below_geoid )

    end subroutine initialize

    subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
        class (type_weight_distribution), intent(in) :: self
        _DECLARE_ARGUMENTS_DO_COLUMN_

        real(rk) :: thickness, integral_weights, present, weights, minimum_value
        real(rk) :: depth,topo, depth_threshold
        real(rk) :: par0, local
        real(rk) :: food 
        real(rk) :: integral  
        real(rk) :: search_food

        integral = 0.0_rk
        integral_weights = 0.0_rk
        depth_threshold = 600.0_rk


        _VERTICAL_LOOP_BEGIN_

            _GET_(self%id_target,local)
            _GET_(self%id_present,present)
            _GET_(self%id_thickness,thickness)
            _GET_(self%id_migrator_food,food)
            _GET_(self%id_depth,depth)
            _GET_SURFACE_(self%id_par0,par0)
            _GET_BOTTOM_(self%id_topo,topo)

            integral = integral + local*thickness

            minimum_value = 0.002 

            ! finalise this
            if (par0 <= 1.0_rk) then
                ! night
                search_food = food
            else
                ! day    
                search_food = food !1.0_rk ! food has no effect
            end if

            thickness = max(thickness, 1.0E-20_rk)

            if (present > 0.5_rk) then
                 weight = (minimum_value + (1.0_rk - minimum_value) * search_food ) 
            else
                if (depth <= min(depth_threshold,topo)) then
                     ! add comment here   
                     weight = (minimum_value  + (0.01_rk - minimum_value) * search_food )  
                else
                     weight = (minimum_value + (0.01_rk - minimum_value) &
                         * exp(-0.025_rk * (depth - min(depth_threshold, topo))) * search_food ) 
                end if 
            end if

            integral_weights = integral_weights + weight*thickness
            
            _SET_DIAGNOSTIC_(self%id_weights,weight)

        _VERTICAL_LOOP_END_

        _SET_BOTTOM_DIAGNOSTIC_(self%id_integral,integral)
        _SET_BOTTOM_DIAGNOSTIC_(self%id_integral_weights,integral_weights)

    end subroutine do_column

end module
