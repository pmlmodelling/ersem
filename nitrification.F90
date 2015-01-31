#include "fabm_driver.h"

module ersem_nitrification

   use fabm_types
   use fabm_particle

   use ersem_shared
   use ersem_pelagic_base

   implicit none

   private

   type,extends(type_ersem_pelagic_base),public :: type_ersem_nitrification
      ! Variables
      type (type_state_variable_id) :: id_O2o
      type (type_state_variable_id) :: id_N3n,id_N4n
      type (type_dependency_id)     :: id_ETW,id_phx

      ! Parameters
      real(rk) :: q10
      real(rk) :: sN4N3X,chN3oX
      integer  :: ISWphx
   contains
!     Model procedures
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_ersem_nitrification),intent(inout),target :: self
      integer,                         intent(in)           :: configunit
!
!EOP
!-----------------------------------------------------------------------
!BOC
      ! Initialize pelagic base model (this also sets the time unit to per day, instead of the default per second)
      call self%initialize_ersem_base(sedimentation=.false.)

      ! Retrieve parameter values
      call self%get_parameter(self%q10,   'q10',  '-',               'Q_10 temperature coefficient')
      call self%get_parameter(self%ISWphx,'ISWph','',                'pH impact on nitrification (0: off, 1: on)')
      call self%get_parameter(self%sN4N3X,'sN4N3','1/d',             'specific nitrification rate')
      call self%get_parameter(self%chN3oX,'chN3o','(mmol O_2/m^3)^3','cubic Michaelis-Menten constant for oxygen dependence of nitrification')

      ! Register links to nutrient and oxygen pools.
      call self%register_state_dependency(self%id_N3n,'N3n','mmol N/m^3',  'nitrate')
      call self%register_state_dependency(self%id_N4n,'N4n','mmol N/m^3',  'ammonium')
      call self%register_state_dependency(self%id_O2o,'O2o','mmol O_2/m^3','oxygen')

      ! Register environmental dependencies (temperature, pH)
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      if (self%ISWphx==1) call self%register_dependency(self%id_phx,standard_variables%ph_reported_on_total_scale)

   end subroutine initialize
   
   subroutine do(self,_ARGUMENTS_DO_)

      class (type_ersem_nitrification),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: ETW,phx,O2o,N4nP
      real(rk) :: etB1,o2state,Fph,fN4N3n

      ! Leave spatial loops (if any)
      _LOOP_BEGIN_

         ! Temperature dependence
         _GET_(self%id_ETW,ETW)
         etB1 = self%q10**((ETW-10._rk)/10._rk) - self%q10**((ETW-32._rk)/3._rk)

         ! Oxygen state for nitrogen species transformation
         _GET_(self%id_O2o,O2o)
         o2state = O2o**3/(O2o**3 + self%chN3oX) ! half saturation    

         !..Nitrification..
         _GET_(self%id_N4n,N4nP)
         fN4N3n = self%sN4N3X  * N4nP * etB1 * o2state

         ! Ph influence on nitrification - empirical equation. Use (1) or not (2)
         if (self%ISWphx==1) then
            _GET_(self%id_phx,phx)
            Fph = min(2._rk,max(0._rk,0.6111_rk*phx-3.8889_rk))
            fN4N3n = fN4N3n * Fph
         end if

         _SET_ODE_(self%id_N3n, + fN4N3n)
         _SET_ODE_(self%id_N4n, - fN4N3n)

      ! Leave spatial loops (if any)
      _LOOP_END_

   end subroutine do
   
end module