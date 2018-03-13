#include "fabm_driver.h"

module ersem_benthic_nitrogen_cycle

   use fabm_types
   use ersem_shared

   implicit none

   private

   ! Model for nitrogen cycle
   type,extends(type_base_model),public :: type_ersem_benthic_nitrogen_cycle
      type (type_bottom_state_variable_id) :: id_K3n,id_K4n,id_G2o,id_K3n2,id_K4n2,id_G2o2,id_G4n,id_benTA,id_benTA2
      type (type_state_variable_id)        :: id_N4n
      type (type_dependency_id)            :: id_ETW,id_ph
      type (type_horizontal_dependency_id) :: id_D1m,id_K6_sms,id_layer2_thickness
      type (type_horizontal_diagnostic_variable_id) :: id_jM3M4n,id_jM3G4n
      type (type_horizontal_diagnostic_variable_id) :: id_nrate
      real(rk) :: q10nit,hM4M3,sM4M3,xno3
      real(rk) :: pammon,pdenit,xn2,hM3G4
      integer :: ISWph
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

   type,extends(type_base_model),public :: type_ersem_K6_calculator
      type (type_horizontal_diagnostic_variable_id) :: id_K6
   end type

contains

   subroutine initialize(self,configunit)
      class (type_ersem_benthic_nitrogen_cycle),intent(inout),target :: self
      integer,                                     intent(in)           :: configunit
      class (type_ersem_K6_calculator),pointer :: child

      ! Set time unit to d-1. This implies that all rates (sink/source terms) are given in d-1.
      self%dt = 86400._rk

      call self%get_parameter(self%q10nit,'q10nit','-',            'Q_10 temperature coefficient')
      call self%get_parameter(self%hM4M3, 'hM4M3', 'mmol/m^3',     'Michaelis-Menten constant for nitrate limitation')
      call self%get_parameter(self%ISWph, 'ISWph', '',             'pH impact on nitrification (0: off, 1: on)',default=0)
      call self%get_parameter(self%sM4M3, 'sM4M3', '1/d',          'maximum nitrification rate at 10 degrees Celsius')
      call self%get_parameter(self%xno3,  'xno3',  'mol O_2/mol N','oxygen consumed per nitrate produced')

      call self%register_diagnostic_variable(self%id_nrate,'nrate','mmol N/m^2/d','benthic nitrification rate',domain=domain_bottom,source=source_do_bottom)
      call self%register_state_dependency(self%id_K3n,'K3n','mmol N/m^2','benthic nitrate in 1st layer')
      call self%register_state_dependency(self%id_K4n,'K4n','mmol N/m^2','benthic ammonium in 1st layer')
      call self%register_state_dependency(self%id_G2o,'G2o','mmol O_2/m^2','benthic oxygen in 1st layer')
      call self%register_state_dependency(self%id_N4n,'N4n','mmol N/m^3','pelagic ammonium')
      if (.not.legacy_ersem_compatibility) then
         call self%register_state_dependency(self%id_benTA,'benTA','mEq/m^2','benthic alkalinity in aerobic layer')
         call self%register_state_dependency(self%id_benTA2,'benTA2','mEq/m^2','benthic alkalinity in anaerobic layer')
      end if
      call self%register_dependency(self%id_D1m,depth_of_bottom_interface_of_layer_1)

      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      if (self%ISWph==1) call self%register_dependency(self%id_ph,standard_variables%ph_reported_on_total_scale)

      ! Denitrification
      call self%get_parameter(self%pammon,'pammon','-','fraction of oxygen demand fulfilled by denitrification under anaerobic conditions')
      call self%get_parameter(self%pdenit,'pdenit','-','fraction of denitrification producing dinitrogen gas (remainder produces ammonium)')
      call self%get_parameter(self%xn2,   'xn2','mol O_2/mol N','oxygen demand fulfilled by reduction of nitrate to dinitrogen gas')
      call self%get_parameter(self%hM3G4,'hM3G4','mmol N/m^3','Michaelis-Menten constant for nitrate limitation of denitrification')

      ! Create our own state avriable for dinitrogen gas
      ! (only to track its total production, which can then be considered in nitrogen mass balance)
      call self%register_state_variable(self%id_G4n,'G4n','mmol N/m^2','dinitrogen gas')
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_G4n)

      call self%register_state_dependency(self%id_K3n2,'K3n2','mmol N/m^2','benthic nitrate in 2nd layer')
      call self%register_state_dependency(self%id_K4n2,'K4n2','mmol N/m^2','benthic ammonium in 2nd layer')
      call self%register_state_dependency(self%id_G2o2,'G2o2','mmol O_2/m^2',  'oxygen in 2nd layer')
      call self%register_dependency(self%id_layer2_thickness,'layer2_thickness','m','thickness of 2nd layer')

! diagnostic fluxes
      call self%register_diagnostic_variable(self%id_jM3M4n,'jM3M4n','mmol N/m^2/d','layer 2 ammonification flux',  source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_jM3G4n,'jM3G4n','mmol N/m^2/d','layer 2 de-nitrification flux',source=source_do_bottom)

      ! Create a child model that provides a K6 diagnostic. Other models (e.g., anaerobic bacteria) can attach to that to provide it with sink/source terms.
      ! In turn, these are then picked up by this model (type_ersem_benthic_nitrogen_cycle) and translated into chnages in NO3 and O2.
      allocate(child)
      call self%add_child(child,'K6_calculator',configunit=configunit)
      call child%register_diagnostic_variable(child%id_K6,'K6','mmol O_2/m^2','oxygen debt due to anaerobic respiration',act_as_state_variable=.true.,output=output_none,domain=domain_bottom,source=source_do_bottom)
      call self%register_dependency(self%id_K6_sms,'K6_sms','mmol O_2/m^2/s','sources-sinks of oxygen debt')
      call self%request_coupling('K6_sms','K6_calculator/K6_sms_tot')
   end subroutine initialize

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_benthic_nitrogen_cycle),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: K3nP,K4nP,N4n
      real(rk) :: ETW,ph
      real(rk) :: Mu_m,eT,eN,Fph,jM4M3n,D1m
      real(rk) :: K6_sms,layer2_thickness,K3n2
      real(rk) :: MU_m2,eN2,jMIno3,jM3M4n,jM3G4n

      _HORIZONTAL_LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_K3n,K3nP)
         _GET_HORIZONTAL_(self%id_K4n,K4nP)
         _GET_(self%id_N4n,N4n)

         _GET_HORIZONTAL_(self%id_D1m,D1m)

         _GET_(self%id_ETW,ETW)

         ! Nitrification (oxygenated layer only):
         ! Reaction: (NH4+,OH-) + 2*O2 -> (NO3-,H+) + 2*H2O
         ! Treated as first order reaction for NH4 in the oxygenated layer, but dampened by presence of nitrate.

         ! Average nitrate concentration in first (oxygenated) sediment layer
         ! (mmol/m3 - per sediment volume, not pore water volume, as porosity is not considered!)
         MU_m = max(0.0_rk,K3nP)/D1m

         ! Limitation factor for temperature (Q_10)
         eT = self%q10nit**((ETW-10._rk)/10._rk)

         ! Limitation factor for nitrate (hyperbolic, 1 at zero nitrate, dropping to 0 at high nitrate).
         eN = self%hM4M3/(self%hM4M3+MU_m)

         ! Ph influence on nitrification - empirical equation
         if(self%ISWph==1) then
            _GET_(self%id_ph,ph)
           Fph = min(2._rk,max(0._rk,0.6111_rk*ph-3.8889_rk))
           jM4M3n = Fph * self%sM4M3 * K4nP * eT * eN
         else
           jM4M3n = self%sM4M3*K4nP*eT*eN
         end if
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_nrate,jM4M3n)

         ! Impose nitrification source-sink terms for benthic ammonium, nitrate, oxygen.
         _SET_BOTTOM_ODE_(self%id_K3n,jM4M3n)
         _SET_BOTTOM_ODE_(self%id_K4n,-jM4M3n)
         _SET_BOTTOM_ODE_(self%id_G2o,-self%xno3*jM4M3n)

         ! Alkalinity contributions: +1 for NH4, -1 for nitrate
         if (.not.legacy_ersem_compatibility) _SET_BOTTOM_ODE_(self%id_benTA,-2*jM4M3n)

         ! Retrieve oxygen debt to to anaerobic respiration, depth-integrated nitrate in oxidized layer, thickness of oxidized layer.
         _GET_HORIZONTAL_(self%id_K6_sms,K6_sms)
         _GET_HORIZONTAL_(self%id_K3n2,K3n2)
         _GET_HORIZONTAL_(self%id_layer2_thickness,layer2_thickness)

         ! FABM provides sources-sinks of K6 in per second - convert to our internal time unit (per day).
         K6_sms = K6_sms * self%dt

         ! Average nitrate concentration in second (denitrification) layer
         ! (mmol/m3 - per sediment volume, not pore water volume, as porosity is not considered!)
         MU_m2 = K3n2/layer2_thickness

         ! Denitrification is a Michaelis-Menten function of nitrate availability
         eN2 = MU_m2/(MU_m2+self%hM3G4)

         ! Convert oxygen demand (mmol O2/m2/d) into nitrate reduction (mmol NO3/m2/d),
         ! accounting for the fractions of nitrate converted to NH4 and to N2.
         jMIno3 = -K6_sms/(self%xno3 *(1._rk - self%pdenit) + self%xn2*self%pdenit)

         ! Partition nitrate reduction flux into NH4 and N2 production fluxes.
         jM3M4n = jMIno3*eN2*self%pammon*(1._rk-self%pdenit)
         jM3G4n = jMIno3*eN2*self%pammon*       self%pdenit

         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_jM3M4n,jM3M4n)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_jM3G4n,jM3G4n)

         _SET_BOTTOM_ODE_(self%id_K4n2, jM3M4n)
         _SET_BOTTOM_ODE_(self%id_K3n2, -jM3M4n -jM3G4n)
         _SET_BOTTOM_ODE_(self%id_G4n, jM3G4n)

         ! Alkalinity contributions: +1 for NH4, -1 for nitrate
         if (.not.legacy_ersem_compatibility) _SET_BOTTOM_ODE_(self%id_benTA2, 2*jM3M4n + jM3G4n)

         ! Oxygen dynamics: after denitrification is taken into account, use actual oxygen to pay off remaining oxygen debt.
         _SET_BOTTOM_ODE_(self%id_G2o2, K6_sms + self%xno3*jM3M4n + self%xn2*jM3G4n)

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module
