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
      type (type_state_variable_id) :: id_O2o,id_TA
      type (type_state_variable_id) :: id_N3n,id_N4n,id_O5n
      type (type_dependency_id)     :: id_ETW,id_phx
      type (type_diagnostic_variable_id)     :: id_nitrification, id_fN4O5n, id_fN4N3n

      ! Parameters
      real(rk) :: q10
      real(rk) :: sN4N3X,chN3oX,chN4nX, N4O5minX
      integer  :: ISWphx
      logical  :: ISWn2o
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
      call self%get_parameter(self%chN3oX,'chN3o','(mmol O_2/m^3)^3','Michaelis-Menten constant for cubic oxygen dependence of nitrification')
      call self%get_parameter(self%ISWn2o,'ISWn2o','',               'activate n2o production', default = .false.)
      if (self%ISWn2o) call self%get_parameter(self%N4O5minX,'N4O5minX','-','minimal fraction of N2O production')
      call self%get_parameter(self%chN4nX,'chN4n','(mmol N/m^3)^3','Michaelis-Menten constant for cubic ammonium dependence of nitrification', default=0.0_rk)

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_nitrification,"rate","mmol/m3/d","rate")
      call self%register_diagnostic_variable(self%id_fN4N3n,'fN4N3n','mmol N/m^3/d','nitrification',output=output_time_step_averaged)
      if (self%ISWn2o) call self%register_diagnostic_variable(self%id_fN4O5n,'fN4O5n','mmol N/m^3/d','N2O production',output=output_time_step_averaged)

      ! Register links to nutrient and oxygen pools.
      call self%register_state_dependency(self%id_N3n,'N3n','mmol N/m^3',  'nitrate')
      call self%register_state_dependency(self%id_N4n,'N4n','mmol N/m^3',  'ammonium')
      call self%register_state_dependency(self%id_O2o,'O2o','mmol O_2/m^3','oxygen')
      if (self%ISWn2o) call self%register_state_dependency(self%id_O5n,'O5n','mmol N/m^3','nitrous oxide')
      call self%register_state_dependency(self%id_TA,standard_variables%alkalinity_expressed_as_mole_equivalent)

      ! Register environmental dependencies (temperature, pH)
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      if (self%ISWphx==1) call self%register_dependency(self%id_phx,standard_variables%ph_reported_on_total_scale)

   end subroutine initialize

   subroutine do(self,_ARGUMENTS_DO_)

      class (type_ersem_nitrification),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: ETW,phx,O2o,N4n,N4nP
      real(rk) :: etB1,o2state,n4state,Fph,fN4N3n,fN4O5n

      ! Leave spatial loops (if any)
      _LOOP_BEGIN_

         ! Temperature dependence
         _GET_(self%id_ETW,ETW)
         etB1 = max(0.0_rk,self%q10**((ETW-10._rk)/10._rk) - self%q10**((ETW-32._rk)/3._rk))

         ! Oxygen state for nitrogen species transformation
         _GET_(self%id_O2o,O2o)
         O2o = max(0.0_rk,O2o)
         o2state = O2o**3/(O2o**3 + self%chN3oX) ! half saturation

         ! Ammonium state for nitrogen species transformation
         _GET_WITH_BACKGROUND_(self%id_N4n,N4n)
         n4state = N4n**3/(N4n**3 + self%chN4nX) ! half saturation

         !..Nitrification..
         _GET_(self%id_N4n,N4nP)
         fN4N3n = self%sN4N3X * N4nP * etB1 * o2state * n4state

         ! Ph influence on nitrification - empirical equation. Use (1) or not (2)
         if (self%ISWphx==1) then
            _GET_(self%id_phx,phx)
            Fph = min(2._rk,max(0._rk,0.6111_rk*phx-3.8889_rk))
            fN4N3n = fN4N3n * Fph
         end if
         ! N2O production
         ! A fixed quota of N is assumed to go to N2O during nitrification
         ! Depending on O2 limitation this quota may increase up to 20 time
         ! (Cadispoti, 2010). Note that N2O is given as mmol of N m-3 (i.e. 1/2
         ! mmole N2O m-2)!  (Luca, July 2016)

         if (self%ISWn2o) then
          fN4O5n=self%N4O5minX*fN4N3n*min(20._rk,(1._rk/o2state))
         _SET_ODE_(self%id_O5n,fN4O5n)
         _SET_ODE_(self%id_N3n, -fN4O5n)
         _SET_DIAGNOSTIC_(self%id_fN4O5n,fN4O5n)
         end if

         _SET_DIAGNOSTIC_(self%id_fN4N3n,fN4N3n)



         _SET_ODE_(self%id_N3n, + fN4N3n)
         _SET_ODE_(self%id_N4n, - fN4N3n)
         _SET_ODE_(self%id_TA, -2*fN4N3n)  ! Alkalinity contributions: +1 for NH4, -1 for nitrate

         ! Legacy ERSEM did not account for oxygen removal by nitrification
         if (.not.legacy_ersem_compatibility) _SET_ODE_(self%id_O2o,-2*fN4N3n)

         _SET_DIAGNOSTIC_(self%id_nitrification,fN4N3n)

      ! Leave spatial loops (if any)
      _LOOP_END_

   end subroutine do

end module
