#include "fabm_driver.h"

module dvm_get_dependencies

use fabm_types
use fabm_expressions

implicit none

private

type, extends(type_base_model), public :: type_get_dependencies

    type (type_diagnostic_variable_id) :: id_migrator_food
    integer  :: nprey
    real(rk) :: divide_food_by
    type (type_state_variable_id),  allocatable,dimension(:) :: id_prey

    contains
        procedure :: initialize
        procedure :: do

end type

contains

    subroutine initialize(self, configunit)
        class (type_get_dependencies), intent(inout), target :: self
        integer, intent(in)                                  :: configunit
        integer                                              :: iprey
        character(len=16)                                    :: index

        call self%get_parameter(self%nprey,'nprey','','number of prey types',default=0)
        call self%get_parameter(self%divide_food_by,'divide_food_by','','a likely concentration (e.g. half saturation constant) to normalize food',default=40.0_rk)
        call self%register_diagnostic_variable(self%id_migrator_food,'migrator_food','mgC/m3','food availability for the migrators', act_as_state_variable=.true., missing_value=0.0_rk, source=source_do)
        ! Get prey-specific coupling links.
        allocate(self%id_prey(self%nprey))
        do iprey=1,self%nprey
            write (index,'(i0)') iprey
            call self%register_state_dependency(self%id_prey(iprey),'prey'//trim(index)//'','mgC/m3', 'prey '//trim(index)//' carbon concentration')
        end do

    end subroutine initialize

    subroutine do(self, _ARGUMENTS_DO_)
        class (type_get_dependencies), intent(in) :: self
        _DECLARE_ARGUMENTS_DO_

        integer  :: iprey
        real(rk),dimension(self%nprey) :: prey

        _LOOP_BEGIN_

            ! Get prey concentrations
            do iprey=1,self%nprey
                _GET_(self%id_prey(iprey), prey(iprey))
            end do
            _SET_DIAGNOSTIC_(self%id_migrator_food, sum(prey)/self%divide_food_by)

        _LOOP_END_

    end subroutine do

end module
