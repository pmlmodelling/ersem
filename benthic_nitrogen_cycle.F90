#include "fabm_driver.h"

module ersem_benthic_nitrogen_cycle

   use fabm_types
   use ersem_shared

   implicit none

   private

   ! Model for nitrogen cycle
   type,extends(type_base_model),public :: type_ersem_benthic_nitrogen_cycle
      type (type_bottom_state_variable_id) :: id_K3n,id_K4n,id_G2o,id_K3n2,id_K4n2,id_G2o2,id_G4n
      type (type_state_variable_id)        :: id_N4n
      type (type_dependency_id)            :: id_ETW,id_phx
      type (type_horizontal_dependency_id) :: id_D1m,id_K6_sms,id_layer2_thickness
      real(rk) :: q10nitX,hM4M3X,sM4M3X,xno3X
      real(rk) :: pammonX,pdenitX,xn2x,hM3G4X
      integer :: ISWphx
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

      call self%get_parameter(self%q10nitX,'q10nitX','-',            'Q_10 temperature coefficient')
      call self%get_parameter(self%hM4M3X, 'hM4M3X', 'mmol/m^3',     'Michaelis-Menten constant for nitrate limitation')
      call self%get_parameter(self%ISWphx, 'ISWphx', '',             'pH impact on nitrification (0: off, 1: on)',default=0)
      call self%get_parameter(self%sM4M3X, 'sM4M3X', '1/d',          'maximum nitrification rate at 10 degrees Celsius')
      call self%get_parameter(self%xno3X,  'xno3X',  'mol O_2/mol N','oxygen consumed per nitrate produced')

      call self%register_state_dependency(self%id_K3n,'K3n','mmol/m^2','nitrate')
      call self%register_state_dependency(self%id_K4n,'K4n','mmol/m^2','ammonium')
      call self%register_state_dependency(self%id_G2o,'G2o','mmol/m^2','oxygen')
      call self%register_state_dependency(self%id_N4n,'N4n','mmol N/m^3','pelagic ammonium')
      call self%register_dependency(self%id_D1m,depth_of_bottom_interface_of_layer_1)

      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      if (self%ISWphx==1) call self%register_dependency(self%id_phx,standard_variables%ph_reported_on_total_scale)

      ! Denitrification
      call self%get_parameter(self%pammonx,'pammonx','-','fraction of oxygen demand fulfilled by denitrification under anaerobic conditions')
      call self%get_parameter(self%pdenitX,'pdenitX','-','fraction of denitrification producing dinitrogen gas (remainder produces ammonium)')
      call self%get_parameter(self%xn2x,   'xn2x','-','oxygen produced per N2 produced')
      call self%get_parameter(self%hM3G4X,'hM3G4X','mmol N/m^3','Michaelis-Menten constant for nitrate limitation of denitrification')

      ! Create our own state avriable for dinitrogen gas
      ! (only to track its total production, which can then be considered in nitrogen mass balance)
      call self%register_state_variable(self%id_G4n,'G4n','mmol N/m^2','dinitrogen gas')
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_G4n)

      call self%register_state_dependency(self%id_K3n2,'K3n2','mmol N/m^2','benthic nitrate in 2nd layer')
      call self%register_state_dependency(self%id_K4n2,'K4n2','mmol N/m^2','benthic ammonium in 2nd layer')
      call self%register_state_dependency(self%id_G2o2,'G2o2','mmol O2/m^2',  'oxygen in 2nd layer')
      call self%register_dependency(self%id_layer2_thickness,'layer2_thickness','m','thickness of 2nd layer')

      ! Create a child model that provides a K6 diagnostic. Other models (e.g., anaerobic bacteria) can attach to that to provide it with sink/source terms.
      ! In turn, these are then picked up by this model (type_ersem_benthic_nitrogen_cycle) and translated into chnages in NO3 and O2.
      allocate(child)
      call self%add_child(child,'K6_calculator',configunit=configunit)
      call child%register_diagnostic_variable(child%id_K6,'K6','mmol O2/m^2','oxygen debt due to anaerobic respiration', act_as_state_variable=.true., output=output_none,domain=domain_bottom)
      call self%register_dependency(self%id_K6_sms,'K6_sms','mmol O2/m^2/s','sources-sinks of oxygen debt')
      call self%request_coupling('K6_sms','K6_calculator/K6_sms_tot')
   end subroutine initialize

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_benthic_nitrogen_cycle),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: K3n,K3nP,K4nP,N4n
      real(rk) :: ETW,phx
      real(rk) :: Mu_m,eT,eN,Fph,jM4M3n,D1m
      real(rk) :: K6_sms,layer2_thickness,K3n2
      real(rk) :: MU_m2,eN2,jMIno3,jM3M4n,jM3G4n

      _HORIZONTAL_LOOP_BEGIN_
         _GET_HORIZONTAL_WITH_BACKGROUND_(self%id_K3n,K3n)
         _GET_HORIZONTAL_(self%id_K3n,K3nP)
         _GET_HORIZONTAL_(self%id_K4n,K4nP)
         _GET_(self%id_N4n,N4n)

         _GET_HORIZONTAL_(self%id_D1m,D1m)
      
         _GET_(self%id_ETW,ETW)

         ! Nitrification (oxygenated layer only):
         ! Reaction: (NH4+,OH-) + 2*O2 -> (NO3-,H+) + 2*H2O
         ! treated as first order reaction for NH4 in the oxygenated layer
         ! Average nitrate density (mmol/m3): MU_m

         ! Average nitrate density in first/oxygenated layer
         MU_m = max(0.0_rk,K3n)/D1m

         ! Limitation factor for temperature (Q_10)
         eT = self%q10nitX**((ETW-10._rk)/10._rk)

         ! Limitation factor for nitrate (hyperbolic, 1 at zero nitrate, dropping to 0 at high nitrate).
         eN = self%hM4M3X/(self%hM4M3X+MU_m)

         ! Ph influence on nitrification - empirical equation
         if(self%ISWphx==1) then
            _GET_(self%id_phx,phx)
           Fph = min(2._rk,max(0._rk,0.6111_rk*phx-3.8889_rk))
           jM4M3n = Fph * self%sM4M3X * K4nP * eT * eN
         else
           jM4M3n = self%sM4M3X*K4nP*eT*eN
         end if

         ! Impose nitrification source-sink terms for benthic ammonium, nitrate, oxygen.
         _SET_BOTTOM_ODE_(self%id_K3n,jM4M3n)
         _SET_BOTTOM_ODE_(self%id_K4n,-jM4M3n)
         _SET_BOTTOM_ODE_(self%id_G2o,-self%xno3X*jM4M3n)

         ! Retrieve oxygen debt to to anaerobic respiration, depth-integrated nitrate in oxidized layer, thickness of oxidized layer.
         _GET_HORIZONTAL_(self%id_K6_sms,K6_sms) 
         _GET_HORIZONTAL_(self%id_K3n2,K3n2)
         _GET_HORIZONTAL_(self%id_layer2_thickness,layer2_thickness)

         ! FABM provides sources-sinks of K6 in per second - convert to our internal time unit (per day).
          K6_sms = K6_sms * self%dt

         ! From nitrate per m2 in layer 2 to nitrate concentration per unit sediment
         ! (not pore water concentration as it does not consider porosity!)
         MU_m2 = K3n2/layer2_thickness
         eN2 = MU_m2/(MU_m2+self%hM3G4X)

         ! Oxygen debt in anaerobic layer:
         ! expression in nitrate reduction (100%):
        jMIno3 = -K6_sms/(self%xno3X *(1._rk - self%pdenitX) + self%xn2X*self%pdenitX)
        jM3M4n = jMIno3*eN2*self%pammonX*(1._rk-self%pdenitX)

        jM3G4n = jMIno3*eN2*self%pammonX*       self%pdenitX

        _SET_BOTTOM_ODE_(self%id_K4n2, jM3M4n)
        _SET_BOTTOM_ODE_(self%id_K3n2, -jM3M4n -jM3G4n)
        _SET_BOTTOM_ODE_(self%id_G4n, jM3G4n)

        ! Oxygen dynamics: after denitrification is taken into account, use actual oxygen to pay off remaining oxygen debt.
        _SET_BOTTOM_ODE_(self%id_G2o2, K6_sms + self%xno3X*jM3M4n + self%xn2X*jM3G4n)

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module
