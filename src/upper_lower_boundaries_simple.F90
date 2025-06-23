#include "fabm_driver.h"

module dvm_upper_lower_boundaries_simple

use fabm_types
use fabm_expressions

implicit none

private 

type, extends(type_base_model), public :: type_upper_lower_boundaries_simple

    type (type_dependency_id)                       :: id_par, id_parmean, id_depth 
    type (type_surface_dependency_id)               :: id_par0
    type (type_diagnostic_variable_id)              :: id_present
    type (type_bottom_dependency_id)                :: id_topo
    real(rk) :: day_light, upper_light, lower_light
    contains
        procedure :: initialize
        procedure :: do

end type

contains

    subroutine initialize(self, configunit)
        class (type_upper_lower_boundaries_simple), intent(inout), target :: self
        integer, intent(in)                                  :: configunit
        
        call self%register_diagnostic_variable(self%id_present,'migrator_presence','-','migrators are present here')
        call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
        call self%register_dependency(self%id_par0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)
        call self%register_dependency(self%id_parmean,temporal_mean(self%id_par,period=86400._rk,resolution=3600._rk,missing_value=1.0_rk))
        call self%register_dependency(self%id_depth,standard_variables%depth)
        call self%register_dependency(self%id_topo,standard_variables%bottom_depth_below_geoid )
        call self%get_parameter( self%day_light,'day_light','log Wm-2','minimum light level for day', default=0._rk)
        call self%get_parameter( self%upper_light,'upper_light','log Wm-2','light level for upper isolume', default=-6.5_rk)
        call self%get_parameter( self%lower_light,'lower_light','log Wm-2','light level for lower isolume',default=-15._rk)

    end subroutine initialize


    subroutine do(self, _ARGUMENTS_DO_)

        class (type_upper_lower_boundaries_simple), intent(in) :: self
        _DECLARE_ARGUMENTS_DO_
    
        real(rk) :: par, par0, parmean
        real(rk) :: parlog, par0log, parmeanlog
        real(rk) :: depth
        real(rk) :: upper_presence, lower_presence
        real(rk) :: topo


        _LOOP_BEGIN_
            _GET_SURFACE_(self%id_par0,par0)
            _GET_BOTTOM_(self%id_topo,topo)

            par0log = max(-20.0_rk, log10(par0))

            _GET_(self%id_par,par)
            _GET_(self%id_parmean,parmean)
            _GET_(self%id_depth,depth)

            parlog = max(-20.0_rk, log10(par))
            parmeanlog = max(-20.0_rk, log10(parmean))

            ! SPECIFY THE POSSIBLE LOCATIONS OF HIGH MIGRATOR CONCENTRATION !

            ! Initialize presence variables
            upper_presence = 0.0_rk
            lower_presence = 0.0_rk

            ! First calculate possibilities above the lower boundary
            ! Lowerlight Rules
            if (parmeanlog >  self%lower_light) then
                upper_presence = 1.0_rk
            else
                upper_presence = 0.0_rk
            end if

            if (par0log > self%day_light) then
                ! Day time
                ! there is an upper and a lower light boundary
                                
                ! Upperlight Rules
                if (parlog < self%upper_light) then
                    lower_presence = 1.0_rk
                else
                    lower_presence = 0.0_rk
                end if
                                
                ! Set diagnostic based on presence
                if (upper_presence + lower_presence > 1.0_rk) then
                    _SET_DIAGNOSTIC_(self%id_present, 1.0_rk)
                else
                    if (upper_presence > 0.9_rk .and. depth >= max(topo - 20.0_rk, 0.0_rk) ) then 
                        ! Within 20m of the bottom - prevents mesozoo being squeezed 
                        ! into a very thin layer at the bottom
                        _SET_DIAGNOSTIC_(self%id_present,1.0_rk)
                    else 
                        _SET_DIAGNOSTIC_(self%id_present, 0.0_rk)
                    end if
                end if
                
            else
                ! Night time
                ! Lower light boundary only
                
                ! Set diagnostic based on presence
                if (upper_presence + lower_presence > 0.9_rk) then
                    _SET_DIAGNOSTIC_(self%id_present, 1.0_rk)
                else
                    _SET_DIAGNOSTIC_(self%id_present, 0.0_rk)
                end if
                        
            end if

        _LOOP_END_

    end subroutine do

end module
