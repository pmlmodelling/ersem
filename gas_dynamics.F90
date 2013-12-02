#include "fabm_driver.h"
module pml_ersem_gas_dynamics

   use fabm_types

   implicit none

   private
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman

   type,extends(type_base_model),public :: type_pml_ersem_gas_dynamics
!     Variable identifiers
      type (type_state_variable_id)     :: id_O3c,id_O2o
      type (type_conserved_quantity_id) :: id_totc
      type (type_dependency_id)         :: id_ETW, id_X1X
      type (type_horizontal_dependency_id) :: id_wnd

      type (type_diagnostic_variable_id) :: id_eO2mO2,id_osat,id_aou

      integer :: ISWO2X
   contains
      procedure :: initialize
      procedure :: do
   end type
      
contains
      
   subroutine initialize(self,configunit)
!
! !INPUT PARAMETERS:
      class (type_pml_ersem_gas_dynamics), intent(inout), target :: self
      integer,                             intent(in)            :: configunit
!
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%ISWO2X,'ISWO2X')

      call self%register_state_variable(self%id_O3c,'O3c','mmol C/m^3','Carbon Dioxide', 2200._rk,minimum=0._rk)
      call self%register_state_variable(self%id_O2o,'O2o','mmol O/m^3','Oxygen',          300._rk,minimum=0._rk)

      call self%register_conserved_quantity(self%id_totc,standard_variables%total_carbon)
      call self%add_conserved_quantity_component(self%id_totc,self%id_O3c)
   
      call self%register_diagnostic_variable(self%id_eO2mO2,'eO2mO2','1','relative oxygen saturation', &
         standard_variable=standard_variables%fractional_saturation_of_oxygen)
      call self%register_diagnostic_variable(self%id_osat,'osat','mmol O/m^3','oxygen saturation concentration')
      call self%register_diagnostic_variable(self%id_aou,'AOU','mmol O/m^3','apparent oxygen utilisation')

      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      call self%register_dependency(self%id_X1X,standard_variables%practical_salinity)
      call self%register_dependency(self%id_wnd,standard_variables%wind_speed)

   end subroutine

   subroutine do(self,_ARGUMENTS_DO_)
      class (type_pml_ersem_gas_dynamics), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: O2o,ETW,X1X,OSAT
      
      _LOOP_BEGIN_
         _GET_(self%id_O2o,O2o)
         _GET_(self%id_ETW,ETW)
         _GET_(self%id_X1X,X1X)
         OSAT = oxygen_saturation_concentration(self,ETW,X1X)
        _SET_DIAGNOSTIC_(self%id_osat,OSAT)
        _SET_DIAGNOSTIC_(self%id_eO2mO2,O2o/OSAT)
        _SET_DIAGNOSTIC_(self%id_aou,OSAT-O2o)
      _LOOP_END_
   end subroutine

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_pml_ersem_gas_dynamics), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: O2o,ETW,X1X,wnd
      real(rk) :: OSAT,ko2o
      
      _HORIZONTAL_LOOP_BEGIN_
         _GET_(self%id_O2o,O2o)
         _GET_(self%id_ETW,ETW)
         _GET_(self%id_X1X,X1X)
         _GET_HORIZONTAL_(self%id_wnd,wnd)
         OSAT = oxygen_saturation_concentration(self,ETW,X1X)
         if (wnd.gt.11._rk) then
            ko2o = sqrt((1953.4_rk-128._rk*etw+3.9918_rk*etw**2-  &
               0.050091_rk*etw**3)/660._rk) * (0.02383_rk * wnd**3._rk)
         else
            ko2o = sqrt((1953.4_rk-128._rk*etw+3.9918_rk*etw**2- &
               0.050091_rk*etw**3)/660._rk) * (0.31_rk * wnd**2._rk)
         endif

! units of ko2 converted from cm/hr to m/day
         ko2o = ko2o*(24._rk/100._rk)

         _SET_SURFACE_EXCHANGE_(self%id_O2o,ko2o*(OSAT-O2o))
      _HORIZONTAL_LOOP_END_
   end subroutine

   function oxygen_saturation_concentration(self,ETW,X1X) result(OSAT)
      class (type_pml_ersem_gas_dynamics), intent(in) :: self
      real(rk),                            intent(in) :: ETW,X1X
      real(rk)                                        :: OSAT

      real(rk) :: A1,A2,A3,A4,B1,B2,B3
      real(rk) :: R,P,T

      DATA A1/-173.4292_rk/,A2/249.6339_rk/,A3/143.3483_rk/,A4/-21.8492_rk/
      DATA B1/-0.033096_rk/,B2/0.014259_rk/,B3/-0.0017_rk/
      DATA R/8.3145_rk/,P/101325_rk/,T/273.15_rk/

      real(rk) :: VIDEAL,ABT

      VIDEAL = (R * 298.15_rk / P) *1000._rk

!  calc absolute temperature
      ABT = ETW + T

      IF (self%ISWO2X.EQ.1) THEN

   !  old ERSEM saturation Formulation

         OSAT = 31.25_rk * (475._rk-2.65_rk*X1X) / (33.5_rk+ETW)

      ELSE IF (self%iswO2X .EQ. 2 ) THEN

   !  calc theoretical oxygen saturation for temp + salinity
   !  From WEISS 1970 DEEP SEA RES 17, 721-735.
   !  units of ln(ml(STP)/l)

         OSAT = A1 + A2 * (100._rk/ABT) + A3 * log(ABT/100._rk) &
         + A4 * (ABT/100._rk) &
         + X1X * ( B1 + B2 * (ABT/100._rk) + B3 * ((ABT/100._rk)**2))

   !  convert units to ml(STP)/l then to mMol/m3

         OSAT = exp( OSAT )
         OSAT = OSAT * 1000._rk / VIDEAL
      ENDIF
   end function
end module