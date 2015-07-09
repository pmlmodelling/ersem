#include "fabm_driver.h"

module ersem_benthic_fauna

   use fabm_types

   use fabm_particle
   use ersem_shared
   use ersem_benthic_base

   implicit none

   private

   type,extends(type_ersem_benthic_base),public :: type_ersem_benthic_fauna
      type (type_state_variable_id)   :: id_O2o
      type (type_bottom_state_variable_id) :: id_Q6c,id_Q6n,id_Q6p,id_Q6s,id_benTA,id_benTA2
      type (type_bottom_state_variable_id) :: id_G3c,id_G2o,id_K4n,id_K1p,id_K4n2,id_K1p2
      type (type_horizontal_dependency_id), allocatable,dimension(:) :: id_foodc,id_foodn,id_foodp,id_foods
      type (type_horizontal_dependency_id), allocatable,dimension(:) :: id_foodc_an
      type (type_dependency_id), allocatable,dimension(:) :: id_foodpelc,id_foodpeln,id_foodpelp,id_foodpels
      type (type_dependency_id) :: id_ETW
      type (type_horizontal_dependency_id) :: id_Dm
      type (type_horizontal_diagnostic_variable_id) :: id_bioirr,id_biotur,id_fYG3c
      type (type_model_id),allocatable,dimension(:) :: id_food

      ! To achieve compatibility with legacy ERSEM, we need to be able to decouple the variable
      ! from which food availability is derived from the variable that absorbs the loss due to
      ! gross food uptake. The following variables absorb the loss due to food uptake - by default
      ! they are coupled to the same variable from which available food is derived.
      type (type_model_id),allocatable,dimension(:) :: id_food_loss_source

      integer  :: nfood    
      real(rk) :: qnc,qpc
      real(rk) :: q10
      real(rk) :: hO2,rlO2
      real(rk) :: xcl,xcs,xch
      real(rk) :: su,lu,hu
      real(rk),allocatable :: pue(:),pufood(:)
      logical,allocatable ::foodispel(:),food_ll(:)
      real(rk) :: pudil
      real(rk) :: sd,sdmO2,sdc,xdc
      real(rk) :: sr,pur
      real(rk) :: ptur,pirr, dwat,dQ6
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self,configunit)
      class (type_ersem_benthic_fauna),intent(inout),target :: self
      integer,                         intent(in)           :: configunit

      integer           :: ifood
      character(len=16) :: index
      logical           :: foodispom
      real(rk)          :: pue,pueQ

      ! Initialize base model type (this will also set the internal time unit to per day)
      call self%initialize_ersem_benthic_base()

      ! Register parameters
      call self%get_parameter(self%qnc, 'qnc',  'mmol N/mg C','nitrogen to carbon ratio')
      call self%get_parameter(self%qpc, 'qpc',  'mmol P/mg C','phosphorus to carbon ratio')
      call self%get_parameter(self%q10,  'q10',  '-',          'Q_10 temperature coefficient')
      call self%get_parameter(self%rlO2, 'rlO2', 'mmol O2/m^3','minimum pelagic oxygen concentration')
      call self%get_parameter(self%hO2,  'hO2',  'mmol O2/m^3','Michaelis-Menten constant for oxygen limitation')
      call self%get_parameter(self%xcl,  'xcl',  'mg C/m^2',   'abundance above which crowding reduces food uptake')
      call self%get_parameter(self%xcs,  'xcs',  'mg C/m^2',   'Michaelis-Menten constant for the impact of crowding')
      call self%get_parameter(self%xch,  'xch',  'mg C/m^2',   'abundance determining asymptotic threshold of crowding limitation (-> xch/(Yc+xch) for Yc-> inf)')
      call self%get_parameter(self%su,   'su',   '1/d',        'specific maximum uptake at reference temperature')
      call self%get_parameter(self%lu,   'lu',   'mg C/m^2',   'Michaelis-Menten constant for food preference as function of food concentration')
      call self%get_parameter(self%hu,   'hu',   'mg C/m^2',   'Michaelis-Menten constant for gross carbon uptake')
      call self%get_parameter(pue,       'pue',  '-',          'fraction of carbon in consumed live food that goes to faeces')
      call self%get_parameter(pueQ,       'pueQ',  '-',          'fraction of carbon in consumed detritus that goes to faeces')
      call self%get_parameter(self%pudil, 'pudil', '-',          'relative nutrient content of faeces')
      call self%get_parameter(self%sd,   'sd',   '1/d',        'specific mortality at reference temperature')
      call self%get_parameter(self%sdmO2,'sdmO2','1/d',        'specific maximum additional mortality due to oxygen stress')
      call self%get_parameter(self%sdc,  'sdc',  '1/d',        'specific maximum additional mortality induced by cold temperature')
      call self%get_parameter(self%xdc,  'xdc',  '1/degree C', 'e-folding temperature factor of cold mortality response')
      call self%get_parameter(self%sr,   'sr',   '1/d',        'specific rest respiration at reference temperature')
      call self%get_parameter(self%pur,  'pur',  '-',          'fraction of assimilated food that is respired')

      ! Add carbon pool as our only state variable.
      call self%add_constituent('c',3000._rk,qn=self%qnc,qp=self%qpc)

      ! Environmental dependencies
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      call self%register_dependency(self%id_Dm,'Dm','m','depth of limiting layer for uptake')

      ! Dependencies on state variables of external modules
      call self%register_state_dependency(self%id_O2o,'O2o','mmol O_2/m^3','pelagic oxygen')
      call self%register_state_dependency(self%id_G3c,'G3c','mmol C/m^2','carbon dioxide')
      call self%register_state_dependency(self%id_G2o,'G2o','mmol O_2/m^2','oxygen')
      call self%register_state_dependency(self%id_K1p, 'K1p', 'mmol P/m^2','benthic phosphate in aerobic layer')
      call self%register_state_dependency(self%id_K4n, 'K4n', 'mmol N/m^2','benthic ammonium in aerobic layer')
      call self%register_state_dependency(self%id_K1p2,'K1p2','mmol P/m^2','benthic phosphate in anaerobic layer')
      call self%register_state_dependency(self%id_K4n2,'K4n2','mmol N/m^2','benthic ammonium in anaerobic layer')
      if (.not.legacy_ersem_compatibility) then
         call self%register_state_dependency(self%id_benTA,'benTA','mEq/m^2','benthic alkalinity in aerobic layer')
         call self%register_state_dependency(self%id_benTA2,'benTA2','mEq/m^2','benthic alkalinity in anaerobic layer')
      end if

      ! Dependencies on state variables of external modules: POM sinks that will take faeces and dead matter.
      call self%register_state_dependency(self%id_Q6c,'Q6c','mg C/m^2',   'particulate organic carbon')
      call self%register_state_dependency(self%id_Q6n,'Q6n','mmol N/m^2', 'particulate organic nitrogen')
      call self%register_state_dependency(self%id_Q6p,'Q6p','mmol P/m^2', 'particulate organic phosphorus')
      call self%register_state_dependency(self%id_Q6s,'Q6s','mmol Si/m^2','particulate organic silicate')

      ! Make it possible to hook up all POM sinks at once by coupling to a whole model "Q".
      call self%request_coupling_to_model(self%id_Q6c,'Q','c')
      call self%request_coupling_to_model(self%id_Q6n,'Q','n')
      call self%request_coupling_to_model(self%id_Q6p,'Q','p')
      call self%request_coupling_to_model(self%id_Q6s,'Q','s')

      ! Determine number of food sources
      call self%get_parameter(self%nfood, 'nfood', '', 'number of food sources',default=0)   

      ! Allocate arrays with food source specific information.
      allocate(self%id_food(self%nfood))
      allocate(self%id_foodpelc(self%nfood))
      allocate(self%id_foodpeln(self%nfood))
      allocate(self%id_foodpelp(self%nfood))
      allocate(self%id_foodpels(self%nfood))
      allocate(self%id_foodc(self%nfood))
      allocate(self%id_foodn(self%nfood))
      allocate(self%id_foodp(self%nfood))
      allocate(self%id_foods(self%nfood))
      allocate(self%id_foodc_an(self%nfood))
      allocate(self%foodispel(self%nfood))
      allocate(self%food_ll(self%nfood))
      allocate(self%id_food_loss_source(self%nfood))

      ! Allocate components of food sources
      do ifood=1,self%nfood
         write (index,'(i0)') ifood
         call self%get_parameter(self%foodispel(ifood),'food'//trim(index)//'ispel','','food source '//trim(index)//' is pelagic',default=.false.)
         call self%get_parameter(self%food_ll(ifood),'food'//trim(index)//'_ll','','availability of food source '//trim(index)//' is limited by aerobic layer height',default=.false.)
         call self%register_model_dependency(self%id_food(ifood),'food'//trim(index))
         if (self%foodispel(ifood)) then
            call self%register_dependency(self%id_foodpelc(ifood), 'food'//trim(index)//'c','mmol C/m^3','food source '//trim(index)//' carbon') 
            call self%register_dependency(self%id_foodpeln(ifood), 'food'//trim(index)//'n','mmol N/m^3','food source '//trim(index)//' nitrogen')
            call self%register_dependency(self%id_foodpelp(ifood), 'food'//trim(index)//'p','mmol P/m^3','food source '//trim(index)//' phosphorus')
            call self%register_dependency(self%id_foodpels(ifood), 'food'//trim(index)//'s','mmol Si/m^3','food source '//trim(index)//' silicate')
            call self%request_coupling_to_model(self%id_foodpelc(ifood),self%id_food(ifood),standard_variables%total_carbon)
            call self%request_coupling_to_model(self%id_foodpeln(ifood),self%id_food(ifood),standard_variables%total_nitrogen)
            call self%request_coupling_to_model(self%id_foodpelp(ifood),self%id_food(ifood),standard_variables%total_phosphorus)
            call self%request_coupling_to_model(self%id_foodpels(ifood),self%id_food(ifood),standard_variables%total_silicate)
         else
            call self%register_dependency(self%id_foodc(ifood),'food'//trim(index)//'c','mmol C/m^2','Food '//trim(index)//' carbon') 
            call self%register_dependency(self%id_foodn(ifood),'food'//trim(index)//'n','mmol N/m^2','Food '//trim(index)//' nitrogen')
            call self%register_dependency(self%id_foodp(ifood),'food'//trim(index)//'p','mmol P/m^2','Food '//trim(index)//' phosphorus')
            call self%register_dependency(self%id_foods(ifood),'food'//trim(index)//'s','mmol Si/m^2','Food '//trim(index)//' silicate')
            call self%request_coupling_to_model(self%id_foodc(ifood),self%id_food(ifood),standard_variables%total_carbon)
            call self%request_coupling_to_model(self%id_foodn(ifood),self%id_food(ifood),standard_variables%total_nitrogen)
            call self%request_coupling_to_model(self%id_foodp(ifood),self%id_food(ifood),standard_variables%total_phosphorus)
            call self%request_coupling_to_model(self%id_foods(ifood),self%id_food(ifood),standard_variables%total_silicate)
            call self%register_dependency(self%id_foodc_an(ifood),'food'//trim(index)//'c_an','mg C/m^2','food source '//trim(index)//' carbon in anaerobic layer')
            call self%request_coupling(self%id_foodc_an(ifood),'zero_hz')
         end if

         ! Legacy ERSEM computes available detritus based on the predator's depth range, but applies the detritus loss
         ! to a detritus pool with a different depth distribution. To be able to reproduce this (inconsistent!) behaviour
         ! we allow separate coupling to the pool that should absorb the detritus loss (this comes in addition to the
         ! pool representing available detritus).
         call self%register_model_dependency(self%id_food_loss_source(ifood),'food'//trim(index)//'_loss_source')
         call self%couplings%set_string('food'//trim(index)//'_loss_source','food'//trim(index))
      end do

      ! Allocate and obtain food source preferences
      allocate(self%pufood(self%nfood))
      do ifood=1,self%nfood
         write (index,'(i0)') ifood
         call self%get_parameter(self%pufood(ifood),'pufood'//trim(index),'-','preference for food source '//trim(index))
      end do

      ! Allocate excretion rates
      allocate(self%pue(self%nfood))
      do ifood=1,self%nfood
         write (index,'(i0)') ifood
         call self%get_parameter(foodispom,'food'//trim(index)//'ispom','','food source '//trim(index)//' is detritus',default=.false.)
         if (foodispom) then
            ! Use assimilation efficiency for particulate organic matter.
            self%pue(ifood) = pueQ

            ! Legacy ERSEM applies the loss of detritus due to feeding to the same pool that absorbs faeces and dead matter,
            ! even though the availability of detritus for consumption is computed over a different depth range.
            ! Implement the legacy behaviour here.
            if (legacy_ersem_compatibility.and..not.self%foodispel(ifood)) &
               call self%couplings%set_string('food'//trim(index)//'_loss_source','Q')
         else
            ! Use assimilation efficiency for living matter.
            self%pue(ifood) = pue
         end if
      end do

      if (any(self%foodispel)) &
         call self%get_parameter(self%dwat, 'dwat', 'm', 'water layer accessible for pelagic food uptake',default=1._rk)
      call self%get_parameter(self%dQ6,  'dQ6',  'm', 'depth of available sediment layer',default=0._rk)

      ! Get contribution for bioturbation and bioirrigation
      call self%get_parameter(self%ptur, 'ptur','-','relative contribution to bioturbation',default=0._rk)
      call self%get_parameter(self%pirr, 'pirr','-','relative contribution to bioirrigation',default=0._rk)
      call self%register_diagnostic_variable(self%id_biotur,'biotur','mg C/m^2/d','bioturbation activity',output=output_none,domain=domain_bottom)
      call self%register_diagnostic_variable(self%id_bioirr,'bioirr','mg C/m^2/d','bioirrigation activity',output=output_none,domain=domain_bottom)
      call self%add_to_aggregate_variable(total_bioturbation_activity, self%id_biotur)
      call self%add_to_aggregate_variable(total_bioirrigation_activity, self%id_bioirr) 

      call self%register_diagnostic_variable(self%id_fYG3c,'fYG3c','mg C/m^2/d','respiration',output=output_time_step_averaged,domain=domain_bottom)

   end subroutine initialize

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)

      class (type_ersem_benthic_fauna),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: Yc,Yn,Yp,O2o,foodP,Dm
      real(rk) :: eT,eO,eC,ETW,Y,x
      real(rk) :: rate
      real(rk) :: SYc,SYn,SYp
      real(rk),dimension(self%nfood) :: foodcP,foodnP,foodpP,foodsP,foodcP_an
      real(rk),dimension(self%nfood) :: feed, sflux, prefcorr
      real(rk) :: foodsum
      real(rk),dimension(self%nfood) :: grossfluxc,grossfluxn,grossfluxp,grossfluxs
      real(rk),dimension(self%nfood) :: netfluxc,netfluxn,netfluxp
      real(rk) :: mort_act,mort_ox,mort_cold,mortflux
      integer  :: ifood,istate
      real(rk) :: fBTYc,nfBTYc,fYG3c,p_an
      real(rk) :: excess_c,excess_n,excess_p

      _HORIZONTAL_LOOP_BEGIN_

      _GET_HORIZONTAL_(self%id_c,Yc)

      _GET_HORIZONTAL_(self%id_Dm,Dm)
      _GET_(self%id_O2o,O2o)

      _GET_(self%id_ETW,ETW)

      Yn = Yc*self%qnc
      Yp = Yc*self%qpc

      ! Calculate temperature limitation factor
      eT = self%q10**((ETW-10._rk)/10._rk) - self%q10**((ETW-32._rk)/3._rk)

      ! Calculate oxygen limitation factor
      eO = (O2o-self%rlO2)**3/((O2o-self%rlO2)**3+(self%hO2-self%rlO2)**3)

      ! Calculate limitaiton factor describing a decrease in feeding rate due to oevrcrowding.
      ! To disable any effect of overcrowding on feeding, set parameter xcl to a very large value.
      ! This can for instance be done to disable the overcrowding effect for meiobenthos/Y4, as in SSB-ERSEM.
      Y = Yc - self%xcl
      if (Y>0._rk) then
         x = Y * Y/(Y+self%xcs)
         eC = 1._rk - x/(x+self%xch)
      else
         eC = 1._rk
      end if

      ! Calculate uptake rate
      rate = self%su * Yc * eT * eO * eC

      ! Get food concentrations: benthic and pelagic!
      do ifood=1,self%nfood
         if (self%foodispel(ifood)) then
            _GET_(self%id_foodpelc(ifood),foodcP(ifood))
            _GET_(self%id_foodpeln(ifood),foodnP(ifood))
            _GET_(self%id_foodpelp(ifood),foodpP(ifood))
            _GET_(self%id_foodpels(ifood),foodsP(ifood))

            foodcP_an(ifood) = 0.0_rk
         else
            _GET_HORIZONTAL_(self%id_foodc(ifood),foodcP(ifood))
            _GET_HORIZONTAL_(self%id_foodn(ifood),foodnP(ifood))
            _GET_HORIZONTAL_(self%id_foodp(ifood),foodpP(ifood))
            _GET_HORIZONTAL_(self%id_foods(ifood),foodsP(ifood))

            _GET_HORIZONTAL_(self%id_foodc_an(ifood),foodcP_an(ifood))
         end if
      end do

      ! Prey carbon was returned in mmol (due to units of standard_variables%total_carbon); convert to mg
      foodcP = foodcP*CMass

      ! In case of pelagic food source, its availability can be limited by a near-bottom fraction available to organism.
      ! In case of benthic food source, it is possible to limit availability of food source by depth of aerobic layer,
      ! e.g. limiting aerobic bacteria as a food source for suspension-feeders by ratio of habitat depth to aerobic layer depth.
      ! Food sources limited in such a way must be specified using logical food{n}_ll in fabm.yaml file.
      do ifood=1,self%nfood
         if (self%foodispel(ifood)) then
            prefcorr(ifood) = self%pufood(ifood) * self%dwat
         else
            if (self%food_ll(ifood)) then
               prefcorr(ifood) = self%pufood(ifood) * min(1._rk,self%dQ6/Dm)
            else
               prefcorr(ifood) = self%pufood(ifood)
            end if
         end if
      end do

      ! Compute effective preferences for individual food sources: "feed".
      ! Weighting factors for original preferences increase hyperbolically from 0 at low food density to 1 at high food density.
      feed = prefcorr * (prefcorr * foodcP/(prefcorr * foodcP + self%lu))

      ! Compute specific uptake rates of the different food sources: "sflux" (1/d).
      ! These are based on a Michaelis-Menten/Type II functional response with dynamic preferences "feed".
      foodsum = sum(feed * foodcP)
      if (foodsum>0._rk) then
         sflux = rate * feed / (foodsum + self%hu)
      else
         sflux = 0._rk
      end if

      ! Gross absolute uptake fluxes (matter/m2/d) per food source.
      grossfluxc = sflux * foodcP
      grossfluxn = sflux * foodnP
      grossfluxp = sflux * foodpP
      grossfluxs = sflux * foodsP

      ! Compute net absolute uptake fluxes (matter/m2/d) per food source from
      ! gross fluxes and assimilation inefficiency.
      netfluxc = grossfluxc * (1._rk -             self%pue)
      netfluxn = grossfluxn * (1._rk - self%pudil*self%pue)
      netfluxp = grossfluxp * (1._rk - self%pudil*self%pue)

      ! Based on relative uptake rate of each food source, decrease all state variables of that food source.
      do ifood=1,self%nfood
         do istate=1,size(self%id_food(ifood)%state)
            _GET_(self%id_food(ifood)%state(istate),foodP)
            _SET_BOTTOM_EXCHANGE_(self%id_food_loss_source(ifood)%state(istate),-sflux(ifood)*foodP)   
         end do
         do istate=1,size(self%id_food(ifood)%bottom_state)
            _GET_HORIZONTAL_(self%id_food(ifood)%bottom_state(istate),foodP)
            _SET_BOTTOM_ODE_(self%id_food_loss_source(ifood)%bottom_state(istate),-sflux(ifood)*foodP)
         end do
      end do

      ! Total gross and net food uptake, summed over all food sources (matter/m2/d).
      fBTYc = sum(grossfluxc)
      nfBTYc = sum(netfluxc)

      ! Store net carbon, nitrogen and phosphorus fluxes associated with food uptake for later use,
      ! and send difference between gross and net fluxes to faeces.
      SYc = nfBTYc
      SYn = sum(netfluxn)
      SYp = sum(netfluxp)
      _SET_BOTTOM_ODE_(self%id_Q6c,fBTYc - nfBTYc)
      _SET_BOTTOM_ODE_(self%id_Q6n,sum(grossfluxn) - sum(netfluxn))
      _SET_BOTTOM_ODE_(self%id_Q6p,sum(grossfluxp) - sum(netfluxp))
      _SET_BOTTOM_ODE_(self%id_Q6s,sum(grossfluxs))

      ! Compute contribution to bioturbation and bioirrigation from gross food uptake.
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_biotur, self%ptur * fBTYc)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_bioirr, self%pirr * fBTYc) 

      ! Respiration fluxes = activity respiration + basal respiration
      fYG3c = self%sr * Yc * eT + self%pur * nfBTYc

      ! Store carbon flux resulting from respiration for later use (note: respiration does not affect nitrogen, phosphorus).
      ! Also account for its production of benthic CO2 and consumption of benthic oxygen.
      SYc = SYc - fYG3c
      _SET_BOTTOM_ODE_(self%id_G3c, fYG3c/CMass)
      _SET_BOTTOM_ODE_(self%id_G2o,-fYG3c/CMass)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fYG3c,fYG3c)

      ! Specific mortality fluxes (1/d): background mortality + mortality due to oxygen limitation + mortality due to cold.
      mort_act  = self%sd * eT
      mort_ox   = self%sdmO2 * (1._rk - eO) * eT
      mort_cold = self%sdc * exp(ETW * self%xdc)

      ! Total specific mortality rate (1/d)
      mortflux = mort_act + mort_ox + mort_cold

      ! Apply mortality to biomass, and send dead matter to particulate organic carbon pool.
      _SET_BOTTOM_ODE_(self%id_c, -mortflux*Yc)
      _SET_BOTTOM_ODE_(self%id_Q6c,mortflux*Yc)
      _SET_BOTTOM_ODE_(self%id_Q6n,mortflux*Yn)
      _SET_BOTTOM_ODE_(self%id_Q6p,mortflux*Yp)

      ! Compute excess carbon flux, given that the maximum realizable carbon flux needs to be balanced
      ! by corresponding nitrogen and phosphorus fluxes to maintain constant stoichiometry.
      excess_c = max(max(SYc - SYp/self%qpc, SYc - SYn/self%qnc), 0.0_rk)
      SYc = SYc - excess_c
      _SET_BOTTOM_ODE_(self%id_c,SYc)
      _SET_BOTTOM_ODE_(self%id_Q6c,excess_c)

      ! Compute excess nitrogen and phosphorus fluxes, based on final carbon flux.
      ! Excess nutrient will be exudated to preserve constant stoichiometry of biomass.
      excess_n = max(SYn - SYc*self%qnc,0.0_rk)
      excess_p = max(SYp - SYc*self%qpc,0.0_rk)

      ! In some cases food or part of it originates from the anaerobic layer.
      ! We distribute excess ammonium and phosphate between aerobic and anaerobic layers proportionally
      ! to the amount of food taken from them. In accordance to the legacy code, excess nutrients by default
      ! enter the aerobic layer. This behaviour can be changed by specifying anaerobic food sources as
      ! food{n}c_an. The sum of food uptake from the anaerobic domain, of those relative to total food
      ! uptake defines how much of the exudated nutrients go into the anaerobic layer.
      p_an = sum(foodcP_an/foodcP*grossfluxc)/max(fBTYc,1.e-8_rk)

      ! Send excess nutrients to phosphate, ammonium pools.
      _SET_BOTTOM_ODE_(self%id_K4n,(1._rk-p_an)*excess_n)
      _SET_BOTTOM_ODE_(self%id_K4n2,p_an       *excess_n)
      _SET_BOTTOM_ODE_(self%id_K1p,(1._rk-p_an)*excess_p)
      _SET_BOTTOM_ODE_(self%id_K1p2,p_an       *excess_p)

      if (.not.legacy_ersem_compatibility) then
         ! Alkalinity contributions: +1 for NH4, -1 for PO4
         _SET_BOTTOM_ODE_(self%id_benTA,(1._rk-p_an)*(excess_n-excess_p)) 
         _SET_BOTTOM_ODE_(self%id_benTA2,p_an       *(excess_n-excess_p))
      end if

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module
