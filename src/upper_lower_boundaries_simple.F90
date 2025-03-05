#include "fabm_driver.h"

module dvm_upper_lower_boundaries_simple

use fabm_types
use fabm_expressions

implicit none

private 

type, extends(type_base_model), public :: type_upper_lower_boundaries_simple

    type (type_dependency_id)                       :: id_par, id_parmean, id_migrator_food, id_depth 
    type (type_horizontal_dependency_id)            :: id_parmean0 
    type (type_horizontal_dependency_id)            :: id_migrator_food0 
    type (type_surface_dependency_id)               :: id_par0
    type (type_diagnostic_variable_id)              :: id_present
    type (type_bottom_dependency_id)                :: id_topo
    type (type_horizontal_dependency_id)            :: id_daylength
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
        call self%register_dependency(self%id_parmean0,temporal_mean(self%id_par0,period=86400._rk,resolution=3600._rk,missing_value=50.0_rk))
        call self%register_dependency(self%id_parmean,temporal_mean(self%id_par,period=86400._rk,resolution=3600._rk,missing_value=1.0_rk))
        call self%register_dependency(self%id_migrator_food,'migrator_food','mgC/m3','food availability for the migrators')
        call self%register_dependency(self%id_migrator_food0,vertical_integral(self%id_migrator_food))
        call self%register_dependency(self%id_depth,standard_variables%pressure)
        call self%register_dependency(self%id_topo,standard_variables%bottom_depth )
        call self%register_dependency(self%id_daylength,'daylength','hours','number of hours light is available at the surface')
        call self%get_parameter( self%day_light,'day_light','log Wm-2','minimum light level for day', default=0._rk)
        call self%get_parameter( self%upper_light,'upper_light','log Wm-2','light level for upper isolume', default=-6.5_rk)
        call self%get_parameter( self%lower_light,'lower_light','log Wm-2','light level for lower isolume',default=-15._rk)

    end subroutine initialize


    subroutine do(self, _ARGUMENTS_DO_)

        class (type_upper_lower_boundaries_simple), intent(in) :: self
        _DECLARE_ARGUMENTS_DO_
    
        real(rk) :: par, par0, parmean, parmean0, nhours, food
        real(rk) :: parlog, par0log, parmeanlog, parmean0log
        real(rk) :: depth
        real(rk) :: upper_presence, lower_presence
        real(rk) :: topo

        _LOOP_BEGIN_

            _GET_SURFACE_(self%id_parmean0,parmean0)
            _GET_SURFACE_(self%id_par0,par0)
            _GET_SURFACE_(self%id_daylength,nhours)
            _GET_SURFACE_(self%id_migrator_food0,food)
            _GET_BOTTOM_(self%id_topo,topo)

            !nhours = min(24.0_rk, max(0.0_rk,nhours * 86400.0_rk))
            par0log = max(-20.0_rk, log10(par0))
            parmean0log = max(-20.0_rk, log10(parmean0))

            _GET_(self%id_par,par)
            _GET_(self%id_parmean,parmean)
            _GET_(self%id_depth,depth)

            parlog = max(-20.0_rk, log10(par))
            parmeanlog = max(-20.0_rk, log10(parmean))

            ! SPECIFY THE POSSIBLE LOCATIONS OF HIGH MIGRATOR CONCENTRATION !

            ! There are 4 cases
            ! 1. Winter Arctic night (surface parmean < 1E-10)
            ! 2. Summer Arctic night (number of daylight hours > 23.9)
            ! 3. Normal day cycle day time
            ! 4. Normal day cycle night time

            ! CASE 1
            upper_presence = 0.0_rk
            lower_presence = 0.0_rk

            if (nhours < 0.1_rk) then
                upper_presence = 0.0_rk
                lower_presence = 0.0_rk
                
                    if (depth < 240.0_rk) then
                        upper_presence = 1.0_rk
                    else
                        upper_presence = 0.0_rk
                    end if
!                end if
                
                ! Set diagnostic based on presence
                if (upper_presence + lower_presence > 0.9_rk) then
                    _SET_DIAGNOSTIC_(self%id_present, 1.0_rk)
                else
                    _SET_DIAGNOSTIC_(self%id_present, 0.0_rk)
                end if
                
            else

                ! CASE 2
                if (nhours > 23.9_rk) then
                    ! there is an upper and a lower light boundary
                    ! first calculate possibilities above the lower boundary

                    ! Initialize presence variables
                    upper_presence = 0.0_rk
                    lower_presence = 0.0_rk
                    
                    ! Lowerlight Rules
                            if (parmeanlog > -15.0_rk) then
                                upper_presence = 1.0_rk
                            else
                                upper_presence = 0.0_rk
                            end if
                    
                    ! Upperlight Rules
                    if (parmean0log <= 0.65_rk) then
                            if (parlog < -8.0_rk) then
                                lower_presence = 1.0_rk
                            else
                                lower_presence = 0.0_rk
                            end if
                    else
                            if (parlog < -4.0_rk) then
                                lower_presence = 1.0_rk
                            else
                                lower_presence = 0.0_rk
                            end if
                    end if
                                        
                    ! Set diagnostic based on presence
                    if (upper_presence + lower_presence > 1.0_rk) then
                        _SET_DIAGNOSTIC_(self%id_present, 1.0_rk)
                    else
                        if (upper_presence > 0.9_rk .and. depth >= max(topo - 20.0_rk, 0.0_rk) ) then 
                            _SET_DIAGNOSTIC_(self%id_present,1.0_rk)
                        else 
                            _SET_DIAGNOSTIC_(self%id_present, 0.0_rk)
                        end if
                    end if
                else

                    ! CASE 3
                    if (par0log > self%day_light) then
                        ! there is an upper and a lower light boundary
                        ! first calculate possibilities above the lower boundary
                        
                        ! Initialize presence variables
                        upper_presence = 0.0_rk
                        lower_presence = 0.0_rk
                        
                        ! Lowerlight Rules
                        if (parmeanlog >  self%lower_light) then
                            upper_presence = 1.0_rk
                        else
                            upper_presence = 0.0_rk
                        end if
                        
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
                                _SET_DIAGNOSTIC_(self%id_present,1.0_rk)
                            else 
                                _SET_DIAGNOSTIC_(self%id_present, 0.0_rk)
                            end if
                        end if
                        
                    else
                        ! CASE 4
                        ! Initialize presence variables
                        upper_presence = 0.0_rk
                        lower_presence = 0.0_rk
                        
                        ! Calculate possibilities above the lower boundary
                        if (parmeanlog > self%lower_light) then
                            upper_presence = 1.0_rk
                        else
                            upper_presence = 0.0_rk
                        end if
                        
                        ! Set diagnostic based on presence
                        if (upper_presence + lower_presence > 0.9_rk) then
                            _SET_DIAGNOSTIC_(self%id_present, 1.0_rk)
                        else
                            _SET_DIAGNOSTIC_(self%id_present, 0.0_rk)
                        end if
                        
                    end if
                end if

            end if

            ! This should ensure that each point at least receives the 0.0_rk value
            if (upper_presence + lower_presence < 0.9_rk) then
                _SET_DIAGNOSTIC_(self%id_present,0.0_rk)     
            end if 
            ! 


            ! ------------------------------------------------------------- !


        _LOOP_END_

    end subroutine do

end module
