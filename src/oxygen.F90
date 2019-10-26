#include "fabm_driver.h"
module ersem_oxygen

! Computes oxygen saturation, apparent oxygen utilization and handles
! exchange of oxygen across the water surface.

! Note: negative oxygen concentrations are permitted.
! These reflect an oygen debt (e.g., presence of H2S)
! In this case, oxygen saturation will be zero (not negative!),
! while apparent oxygen utilization will still be the difference
! between saturation concentration and [negative] oxygen concentration.
! Thus, the oxygen debt is included as part of utilization.

   use fabm_types
   use ersem_shared

   implicit none

   private

   type,extends(type_base_model),public :: type_ersem_oxygen
!     Variable identifiers
      type (type_state_variable_id)     :: id_O2o
      type (type_dependency_id)         :: id_ETW, id_X1X
      type (type_horizontal_dependency_id) :: id_wnd

      type (type_diagnostic_variable_id) :: id_eO2mO2,id_osat,id_aou
      type (type_horizontal_diagnostic_variable_id) ::  id_fair

      integer :: iswO2X, iswASFLUX
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
   end type

contains

   subroutine initialize(self,configunit)
!
! !INPUT PARAMETERS:
      class (type_ersem_oxygen), intent(inout), target :: self
      integer,                       intent(in)            :: configunit
!
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%iswO2X,'iswO2','','saturation formulation (1: legacy ERSEM, 2: Weiss 1970, 3: Nightingale et al. 2000, 4: Wanninkhof 2014)')
      call self%get_parameter(self%iswASFLUX,'iswASFLUX','','air-sea O2 exchange (0: none, 1: Nightingale et al. 2000, 2: Wanninkhof 1992 without chemical enhancement, 3: Wanninkhof 1992 with chemical enhancement, 4: Wanninkhof and McGillis 1999, 5: Wanninkhof 1992 switching to Wanninkhof and McGillis 1999, 6: Wanninkhof 2014)',default=6)

      call self%register_state_variable(self%id_O2o,'o','mmol O_2/m^3','oxygen',300._rk)

      call self%register_diagnostic_variable(self%id_eO2mO2,'eO2mO2','1','relative saturation', &
         standard_variable=standard_variables%fractional_saturation_of_oxygen)
      call self%register_diagnostic_variable(self%id_osat,'osat','mmol O_2/m^3','saturation concentration')
      call self%register_diagnostic_variable(self%id_aou,'AOU','mmol O_2/m^3','apparent utilisation')
      call self%register_diagnostic_variable(self%id_fair,'fair','mmol O_2/m^2/d','air-sea flux', source=source_do_surface)

      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      call self%register_dependency(self%id_X1X,standard_variables%practical_salinity)
      call self%register_dependency(self%id_wnd,standard_variables%wind_speed)

      self%dt = 3600._rk*24._rk

   end subroutine

   subroutine do(self,_ARGUMENTS_DO_)
      class (type_ersem_oxygen), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: O2o,ETW,X1X,OSAT

      _LOOP_BEGIN_
         _GET_(self%id_O2o,O2o)
         _GET_(self%id_ETW,ETW)
         _GET_(self%id_X1X,X1X)
         OSAT = oxygen_saturation_concentration(self,ETW,X1X)
         _SET_DIAGNOSTIC_(self%id_osat,OSAT)
         _SET_DIAGNOSTIC_(self%id_eO2mO2,max(0.0_rk,O2o/OSAT))
         _SET_DIAGNOSTIC_(self%id_aou,OSAT-O2o)
      _LOOP_END_
   end subroutine

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_ersem_oxygen), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: O2o,ETW,T,X1X,wnd
      real(rk) :: OSAT,sc,ko2o,FAIRO2

      _HORIZONTAL_LOOP_BEGIN_
         _GET_(self%id_O2o,O2o)
         _GET_(self%id_ETW,ETW)
         _GET_(self%id_X1X,X1X)
         _GET_HORIZONTAL_(self%id_wnd,wnd)

         wnd = max(wnd, 0.0_rk)

         OSAT = oxygen_saturation_concentration(self,ETW,X1X)

         if (legacy_ersem_compatibility) then
           ! Old formulation for the Schmidt number,
           ! left for documentation and back compatibility
           sc = max(0.0_rk, 1953.4_rk - 128._rk*ETW + 3.9918_rk*ETW**2 - 0.050091_rk*ETW**3)
         else
           ! New formulation for the Schmidt number for O2 following Wanninkhof 2014
           T = max(min(ETW,40.0_rk), -2.0_rk)
           sc = 1920.4_rk - 135.6_rk*T + 5.2122_rk*T**2._rk - 0.10939_rk*T**3 + 0.00093777_rk*T**4._rk
         endif

         if     (self%iswASFLUX==1) then      ! Nightingale et al. 2000
               ko2o = (0.222_rk*wnd**2.0_rk + 0.333_rk*wnd)*(sc/600._rk)**(-0.5_rk)
         elseif (self%iswASFLUX==2) then      ! Wanninkhof 1992 without chemical enhancement
               ko2o = 0.31_rk*wnd**2.0_rk*(sc/660._rk)**(-0.5_rk)
         elseif (self%iswASFLUX==3) then      ! Wanninkhof 1992 with chemical enhancement
               ko2o = (2.5_rk*(0.5246_rk+0.016256_rk*T+0.00049946_rk*T**2.0_rk) + 0.3_rk*wnd**2.0_rk)*(sc/660._rk)**(-0.5_rk)
         elseif (self%iswASFLUX==4) then      ! Wanninkhof and McGillis 1999
               ko2o = 0.0283_rk*wnd**3.0_rk*(sc/660._rk)**(-0.5_rk)
         elseif (self%iswASFLUX==5) then      ! Wanninkhof 1992 switching to Wanninkhof and McGillis 1999
            if (wnd.gt.11._rk) then
               ko2o = 0.0283_rk*wnd**3.0_rk*(sc/660._rk)**(-0.5_rk)  ! Wanninkhof & McGillis 1999
               ! Note: Old ERSEM had a typo of 0.02383 instead of 0.0283
            else
               ko2o = 0.31_rk*wnd**2.0_rk*(sc/660._rk)**(-0.5_rk)    ! Wanninkhof 1992
            endif
         elseif (self%iswASFLUX==6) then      ! Wanninkhof 2014
               ko2o = 0.251_rk*wnd**2.0_rk*(sc/660._rk)**(-0.5_rk)
         endif

         ! units of ko2 converted from cm/hr to m/day
         ko2o = ko2o*(24._rk/100._rk)

         FAIRO2 = ko2o*(OSAT-O2o)

         _SET_SURFACE_EXCHANGE_(self%id_O2o,FAIRO2)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fair,FAIRO2)
      _HORIZONTAL_LOOP_END_
   end subroutine

   function oxygen_saturation_concentration(self,ETW,X1X) result(OSAT)
      class (type_ersem_oxygen), intent(in) :: self
      real(rk),                      intent(in) :: ETW,X1X
      real(rk)                                  :: OSAT

      real(rk),parameter :: A1 = -173.4292_rk
      real(rk),parameter :: A2 = 249.6339_rk
      real(rk),parameter :: A3 = 143.3483_rk
      real(rk),parameter :: A4 = -21.8492_rk
      real(rk),parameter :: B1 = -0.033096_rk
      real(rk),parameter :: B2 = 0.014259_rk
      real(rk),parameter :: B3 = -0.0017_rk
      real(rk),parameter :: R = 8.3145_rk
      real(rk),parameter :: P = 101325_rk
      real(rk),parameter :: T = 273.15_rk

      ! volume of an ideal gas at standard temp (0C) and pressure (1 atm)
      real(rk),parameter :: VIDEAL = (R * 273.15_rk / P) *1000._rk

      real(rk)           :: ABT

      ! calc absolute temperature
      ABT = ETW + T

      if (self%ISWO2X==1) then
         ! old ERSEM saturation Formulation
         OSAT = 31.25_rk * (475._rk-2.65_rk*X1X) / (33.5_rk+ETW)
      else
         ! calc theoretical oxygen saturation for temp + salinity
         ! From WEISS 1970 DEEP SEA RES 17, 721-735.
         ! units of ln(ml(STP)/l)
         OSAT = A1 + A2 * (100._rk/ABT) + A3 * log(ABT/100._rk) &
                  + A4 * (ABT/100._rk) &
                  + X1X * ( B1 + B2 * (ABT/100._rk) + B3 * ((ABT/100._rk)**2))

         ! convert units to ml(STP)/l then to mMol/m3
         OSAT = exp( OSAT )
         OSAT = OSAT * 1000._rk / VIDEAL
      end if
   end function

end module
