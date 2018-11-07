#include "fabm_driver.h"

module ersem_benthic_carbonate

   use fabm_types
   use ersem_carbonate

   implicit none

   private

   type,extends(type_base_model),public :: type_ersem_benthic_carbonate
!     Variable identifiers
      type (type_horizontal_dependency_id)     :: id_G3c,id_benTA
      type (type_dependency_id)         :: id_ETW, id_X1X, id_dens, id_pres
      type (type_horizontal_dependency_id)         :: id_Carb_in,id_pco2_in

      type (type_horizontal_diagnostic_variable_id) :: id_ph,id_pco2,id_CarbA, id_BiCarb, id_Carb
      type (type_horizontal_diagnostic_variable_id) :: id_Om_cal,id_Om_arg
      integer :: phscale

   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self,configunit)

      class (type_ersem_benthic_carbonate), intent(inout), target :: self
      integer,                      intent(in)            :: configunit

      call self%get_parameter(self%phscale,'pHscale','','pH scale (0: SWS, 1: total)',default=1,minimum=-1,maximum=1)
      call self%register_horizontal_dependency(self%id_G3c,'G3c','mmol C/m^2','carbon dioxide')
      call self%register_horizontal_dependency(self%id_benTA,'benTA','mmol eq/m^2','benthic alkalinity')
      if (self%phscale==1) then
             call self%register_diagnostic_variable(self%id_ph, 'pH', '-', 'pH on total scale',missing_value=0._rk,domain=domain_bottom,source=source_do_bottom)
      elseif (self%phscale==0) then
             call self%register_diagnostic_variable(self%id_ph, 'pH', '-', 'pH on seawater scale',missing_value=0._rk,domain=domain_bottom,source=source_do_bottom)
      elseif (self%phscale==-1) then
             call self%register_diagnostic_variable(self%id_ph, 'pH', '-', 'pH on seawater scale',missing_value=0._rk,domain=domain_bottom,source=source_do_bottom)
      end if
      call self%register_diagnostic_variable(self%id_pco2,  'pCO2',  '1e-6',    'partial pressure of CO2',source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_CarbA, 'CarbA', 'mmol/m^3','carbonic acid concentration',source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_BiCarb,'BiCarb','mmol/m^3','bicarbonate concentration',source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_Carb,  'Carb',  'mmol/m^3','carbonate concentration',source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_Om_cal,'Om_cal','-','calcite saturation',source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_Om_arg,'Om_arg','-','aragonite saturation',source=source_do_bottom)
      call self%register_dependency(self%id_ETW, standard_variables%temperature)
      call self%register_dependency(self%id_X1X, standard_variables%practical_salinity)
      call self%register_dependency(self%id_dens,standard_variables%density)
      call self%register_dependency(self%id_pres,standard_variables%pressure)
      call self%register_dependency(self%id_pco2_in,'pCO2','1e-6','previous pCO2')
      call self%register_dependency(self%id_Carb_in,'Carb','mmol/m^3','previous carbonate concentration')

   end subroutine

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_benthic_carbonate), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: O3c,ETW,X1X,density,pres
      real(rk) :: TA,Ctot
      real(rk) :: pH,PCO2,H2CO3,HCO3,CO3,k0co2
      real(rk) :: Om_cal,Om_arg
      logical  :: success

      _HORIZONTAL_LOOP_BEGIN_

         _GET_HORIZONTAL_(self%id_G3C,O3C)
         _GET_(self%id_ETW,ETW)
         _GET_(self%id_X1X,X1X)
         _GET_(self%id_dens,density)
         _GET_(self%id_pres,pres)
         _GET_HORIZONTAL_(self%id_benTA,TA)

         TA = TA/1.0e6_rk                ! from umol kg-1 to mol kg-1
         Ctot  = O3C / 1.e3_rk / density ! from mmol m-3 to mol kg-1
         call co2dyn (ETW,X1X,pres*0.1_rk,ctot,TA,pH,PCO2,H2CO3,HCO3,CO3,k0co2,success,self%phscale)   ! NB pressure from dbar to bar
         if (success) then
            ! Carbonate system iterative scheme converged -  save associated diagnostics.
            ! Convert outputs from fraction to ppm (pCO2) and from mol kg-1 to mmol m-3 (concentrations).
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ph,pH)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_pco2,PCO2*1.e6_rk)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_CarbA, H2CO3*1.e3_rk*density)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Bicarb,HCO3*1.e3_rk*density)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Carb,  CO3*1.e3_rk*density)
         else
            ! Carbonate system iterative scheme did not converge.
            ! All diagnostics retain their previous value.
            ! Use previous carbonate concentration (but current environment) for carbonate saturation states.
            _GET_HORIZONTAL_(self%id_Carb_in,CO3)
            CO3 = CO3/1.e3_rk/density  ! from mmol/m3 to mol/kg
         end if

         ! Call carbonate saturation state subroutine
         call CaCO3_Saturation (ETW, X1X, pres*1.e4_rk, CO3, Om_cal, Om_arg)  ! NB pressure from dbar to Pa
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Om_cal,Om_cal)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Om_arg,Om_arg)

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module
