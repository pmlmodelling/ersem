#include "fabm_driver.h"

! -----------------------------------------------------------------------------
! For diel vertical migration of migrating plankton. Calculates an upper and 
! lower boundary in the water column, between which the migrator is able to 
! be present. Boundaries depend on light levels which can be set by the user
! or defaul to a maximum light level of -6.5 log W m-2 and a minimum light 
! level of -15 log W m-2
!
! Adapted from code written by Caglar Yumruktepe (NERSC), available at: 
! https://github.com/nansencenter/nersc
! -----------------------------------------------------------------------------

module dvm_upper_lower_boundaries

use fabm_types
use fabm_expressions

implicit none

private 

type, extends(type_base_model), public :: type_upper_lower_boundaries

    type (type_dependency_id)                       :: id_par, id_parmean, id_depth 
    type (type_surface_dependency_id)               :: id_par0
    type (type_bottom_dependency_id)                :: id_topo
    type (type_state_variable_id),  allocatable,dimension(:) :: id_prey !state variable here - dependency in ersem mesozoo
    type (type_diagnostic_variable_id)              :: id_present
    type (type_diagnostic_variable_id)              :: id_migrator_food

    integer  :: nprey
    real(rk) :: day_light, upper_light, lower_light
    real(rk) :: divide_food_by
    contains
        procedure :: initialize
        procedure :: do

end type

contains

    subroutine initialize(self, configunit)
        class (type_upper_lower_boundaries), intent(inout), target :: self
        integer, intent(in)                                               :: configunit
        integer                                                           :: iprey
        character(len=16)                                                 :: index
        call self%register_diagnostic_variable(self%id_present,'migrator_presence','-','migrators are present here')
        call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
        call self%register_dependency(self%id_par0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)
        call self%register_dependency(self%id_parmean,temporal_mean(self%id_par,period=86400._rk,resolution=3600._rk,missing_value=1.0_rk))
        call self%register_dependency(self%id_depth,standard_variables%depth)
        call self%register_dependency(self%id_topo,standard_variables%bottom_depth_below_geoid )
        call self%get_parameter( self%day_light,'day_light','log W m-2','minimum light level for day', default=0._rk)
        call self%get_parameter( self%upper_light,'upper_light','log W m-2','light level for upper isolume', default=-6.5_rk)
        call self%get_parameter( self%lower_light,'lower_light','log W m-2','light level for lower isolume',default=-15._rk)

        call self%get_parameter(self%nprey,'nprey','','number of prey types',default=1)
        call self%get_parameter(self%divide_food_by,'divide_food_by','','a likely concentration (e.g. half saturation constant) to normalize food',default=40.0_rk)
        call self%register_diagnostic_variable(self%id_migrator_food,'migrator_food','mmol C/m^3','food availability for the migrators', act_as_state_variable=.true., missing_value=0.0_rk, source=source_do)
        ! Get prey-specific coupling links.
        allocate(self%id_prey(self%nprey))
        do iprey=1,self%nprey
            write (index,'(i0)') iprey
            call self%register_state_dependency(self%id_prey(iprey),'prey'//trim(index)//'','mmol C/m^3', 'prey '//trim(index)//' carbon')
        end do

    end subroutine initialize


    subroutine do(self, _ARGUMENTS_DO_)

        class (type_upper_lower_boundaries), intent(in) :: self
        _DECLARE_ARGUMENTS_DO_
    
        real(rk) :: par, par0, parmean
        real(rk) :: parlog, par0log, parmeanlog
        real(rk) :: depth
        real(rk) :: upper_presence, lower_presence
        real(rk) :: topo
        integer  :: iprey
        real(rk),dimension(self%nprey) :: prey


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

            ! Lower boundary exists at all times of the day, so first calculate
            ! possibilities above the lower boundary, which depends on the 
            ! minimum light level (lower_light)
            if (parmeanlog >  self%lower_light) then
                upper_presence = 1.0_rk
            else
                upper_presence = 0.0_rk
            end if

            if (par0log > self%day_light) then
                ! Day time
                ! Migrators are not at the surface and therefore there is an
                ! upper and a lower light boundary
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

            ! Migrators will be distributed in water column according to prey
            ! availability
            ! Get prey concentrations to calculate prey availability
            do iprey=1,self%nprey
                _GET_(self%id_prey(iprey), prey(iprey))
            end do
            _SET_DIAGNOSTIC_(self%id_migrator_food, sum(prey)/self%divide_food_by)

        _LOOP_END_
    
    end subroutine do

end module
