#include "fabm_driver.h"

module ersem_benthic_fauna

   use fabm_types

   use fabm_particle
   use ersem_shared
   use ersem_benthic_base

   implicit none

   private

   type type_food
      type (type_model_id)                          :: id_source
      type (type_horizontal_dependency_id)          :: id_c,id_n,id_p,id_s
      type (type_horizontal_dependency_id)          :: id_c_an
      type (type_dependency_id)                     :: id_c_pel,id_n_pel,id_p_pel,id_s_pel
      type (type_horizontal_diagnostic_variable_id) :: id_fc,id_fn,id_fp,id_fs

      ! To achieve compatibility with legacy ERSEM, we need to be able to decouple the variable
      ! from which food availability is derived from the variable that absorbs the loss due to
      ! gross food uptake. The following variables absorb the loss due to food uptake - by default
      ! they are coupled to the same variable from which available food is derived.
      type (type_model_id) :: id_loss_source

      real(rk) :: pu
      real(rk) :: pue
      logical  :: ispel
      logical  :: ll
   end type

   type,extends(type_ersem_benthic_base),public :: type_ersem_benthic_fauna
      type (type_state_variable_id)   :: id_O2o
      type (type_bottom_state_variable_id) :: id_Q6c,id_Q6n,id_Q6p,id_Q6s,id_benTA,id_benTA2
      type (type_bottom_state_variable_id) :: id_G3c,id_G2o,id_K4n,id_K1p,id_K4n2,id_K1p2
      type (type_horizontal_dependency_id) :: id_Dm
      type (type_horizontal_diagnostic_variable_id) :: id_bioirr,id_biotur,id_fYG3c, id_fYKIn,id_fYK1p,id_fYQPc,id_fYQPn,id_fYQPp
      type (type_food), allocatable :: food(:)

      integer  :: nfood
      real(rk) :: qnc,qpc
      real(rk) :: hO2,rlO2
      real(rk) :: xcl,xcs,xch
      real(rk) :: su,lu,hu
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
      real(rk)          :: c0

      ! Initialize base model type (this will also set the internal time unit to per day)
      call self%initialize_ersem_benthic_base()

      ! Register parameters
      call self%get_parameter(self%qnc,  'qnc',  'mmol N/mg C','nitrogen to carbon ratio')
      call self%get_parameter(self%qpc,  'qpc',  'mmol P/mg C','phosphorus to carbon ratio')
      call self%get_parameter(c0,        'c0',   'mg C/m^2',   'background concentration',default=0.0_rk)
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
      call self%get_parameter(pueQ,      'pueQ', '-',          'fraction of carbon in consumed detritus that goes to faeces')
      call self%get_parameter(self%pudil,'pudil','-',          'relative nutrient content of faeces')
      call self%get_parameter(self%sd,   'sd',   '1/d',        'specific mortality at reference temperature')
      call self%get_parameter(self%sdmO2,'sdmO2','1/d',        'specific maximum additional mortality due to oxygen stress')
      call self%get_parameter(self%sdc,  'sdc',  '1/d',        'specific maximum additional mortality induced by cold temperature')
      call self%get_parameter(self%xdc,  'xdc',  '1/degree_C', 'e-folding temperature factor of cold mortality response')
      call self%get_parameter(self%sr,   'sr',   '1/d',        'specific rest respiration at reference temperature')
      call self%get_parameter(self%pur,  'pur',  '-',          'fraction of assimilated food that is respired')

      ! Add carbon pool as our only state variable.
      call self%add_constituent('c',3000._rk,c0,qn=self%qnc,qp=self%qpc)

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
      allocate(self%food(self%nfood))
      do ifood=1,self%nfood
         write (index,'(i0)') ifood
         call self%get_parameter(self%food(ifood)%ispel,'food'//trim(index)//'ispel','','food source '//trim(index)//' is pelagic',default=.false.)
         call self%get_parameter(self%food(ifood)%ll,'food'//trim(index)//'_ll','','availability of food source '//trim(index)//' is limited by aerobic layer height',default=.false.)
         call self%register_model_dependency(self%food(ifood)%id_source,'food'//trim(index))
         if (self%food(ifood)%ispel) then
            call self%register_dependency(self%food(ifood)%id_c_pel, 'food'//trim(index)//'c','mmol C/m^3', 'carbon in food source '//trim(index))
            call self%register_dependency(self%food(ifood)%id_n_pel, 'food'//trim(index)//'n','mmol N/m^3', 'nitrogen in food source '//trim(index))
            call self%register_dependency(self%food(ifood)%id_p_pel, 'food'//trim(index)//'p','mmol P/m^3', 'phosphorus in food source '//trim(index))
            call self%register_dependency(self%food(ifood)%id_s_pel, 'food'//trim(index)//'s','mmol Si/m^3','silicate in food source '//trim(index))
            call self%request_coupling_to_model(self%food(ifood)%id_c_pel,self%food(ifood)%id_source,standard_variables%total_carbon)
            call self%request_coupling_to_model(self%food(ifood)%id_n_pel,self%food(ifood)%id_source,standard_variables%total_nitrogen)
            call self%request_coupling_to_model(self%food(ifood)%id_p_pel,self%food(ifood)%id_source,standard_variables%total_phosphorus)
            call self%request_coupling_to_model(self%food(ifood)%id_s_pel,self%food(ifood)%id_source,standard_variables%total_silicate)
         else
            call self%register_dependency(self%food(ifood)%id_c,'food'//trim(index)//'c','mmol C/m^2', 'carbon in food source '//trim(index))
            call self%register_dependency(self%food(ifood)%id_n,'food'//trim(index)//'n','mmol N/m^2', 'nitrogen in food source '//trim(index))
            call self%register_dependency(self%food(ifood)%id_p,'food'//trim(index)//'p','mmol P/m^2', 'phosphorus in food source '//trim(index))
            call self%register_dependency(self%food(ifood)%id_s,'food'//trim(index)//'s','mmol Si/m^2','silicate in food source '//trim(index))
            call self%request_coupling_to_model(self%food(ifood)%id_c,self%food(ifood)%id_source,standard_variables%total_carbon)
            call self%request_coupling_to_model(self%food(ifood)%id_n,self%food(ifood)%id_source,standard_variables%total_nitrogen)
            call self%request_coupling_to_model(self%food(ifood)%id_p,self%food(ifood)%id_source,standard_variables%total_phosphorus)
            call self%request_coupling_to_model(self%food(ifood)%id_s,self%food(ifood)%id_source,standard_variables%total_silicate)
            call self%register_dependency(self%food(ifood)%id_c_an,'food'//trim(index)//'c_an','mg C/m^2','food source '//trim(index)//' carbon in anaerobic layer')
            call self%request_coupling(self%food(ifood)%id_c_an,'zero_hz')
         end if
         call self%register_diagnostic_variable(self%food(ifood)%id_fc,'fprey'//trim(index)//'c','mg C/m^2/d',   'uptake of carbon in food source '//trim(index),    source=source_do_bottom)
         call self%register_diagnostic_variable(self%food(ifood)%id_fn,'fprey'//trim(index)//'n','mmol N/m^2/d', 'uptake of nitrogen in food source '//trim(index),  source=source_do_bottom)
         call self%register_diagnostic_variable(self%food(ifood)%id_fp,'fprey'//trim(index)//'p','mmol P/m^2/d', 'uptake of phosphorus in food source '//trim(index),source=source_do_bottom)
         call self%register_diagnostic_variable(self%food(ifood)%id_fs,'fprey'//trim(index)//'s','mmol Si/m^2/d','uptake of silicate in food source '//trim(index),  source=source_do_bottom)

         ! Legacy ERSEM computes available detritus based on the predator's depth range, but applies the detritus loss
         ! to a detritus pool with a different depth distribution. To be able to reproduce this (inconsistent!) behaviour
         ! we allow separate coupling to the pool that should absorb the detritus loss (this comes in addition to the
         ! pool representing available detritus).
         call self%register_model_dependency(self%food(ifood)%id_loss_source,'food'//trim(index)//'_loss_source')
         call self%couplings%set_string('food'//trim(index)//'_loss_source','food'//trim(index))
      end do

      ! Obtain food source preferences (separate loop to ensure all preference are shown together)
      do ifood=1,self%nfood
         write (index,'(i0)') ifood
         call self%get_parameter(self%food(ifood)%pu,'pufood'//trim(index),'-','preference for food source '//trim(index))
      end do

      ! Set food-source-specific assimilation inefficiency.
      ! (fraction of ingested frood going to faeces)
      do ifood=1,self%nfood
         write (index,'(i0)') ifood
         call self%get_parameter(foodispom,'food'//trim(index)//'ispom','','food source '//trim(index)//' is detritus',default=.false.)
         if (foodispom) then
            ! Use assimilation efficiency for particulate organic matter.
            self%food(ifood)%pue = pueQ

            ! Legacy ERSEM applies the loss of detritus due to feeding to the same pool that absorbs faeces and dead matter,
            ! even though the availability of detritus for consumption is computed over a different depth range.
            ! Implement the legacy behaviour here.
            if (legacy_ersem_compatibility.and..not.self%food(ifood)%ispel) &
               call self%couplings%set_string('food'//trim(index)//'_loss_source','Q')
         else
            ! Use assimilation efficiency for living matter.
            self%food(ifood)%pue = pue
         end if
      end do

      if (any(self%food%ispel)) &
         call self%get_parameter(self%dwat, 'dwat', 'm', 'water layer accessible for pelagic food uptake',default=1._rk)
      call self%get_parameter(self%dQ6,  'dQ6',  'm', 'depth of available sediment layer',default=0._rk)

      ! Get contribution for bioturbation and bioirrigation
      call self%get_parameter(self%ptur, 'ptur','-','relative contribution to bioturbation',default=0._rk)
      call self%get_parameter(self%pirr, 'pirr','-','relative contribution to bioirrigation',default=0._rk)
      call self%register_diagnostic_variable(self%id_biotur,'biotur','mg C/m^2/d','bioturbation activity',output=output_none,domain=domain_bottom,source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_bioirr,'bioirr','mg C/m^2/d','bioirrigation activity',output=output_none,domain=domain_bottom,source=source_do_bottom)
      call self%add_to_aggregate_variable(total_bioturbation_activity, self%id_biotur)
      call self%add_to_aggregate_variable(total_bioirrigation_activity, self%id_bioirr)

      call self%register_diagnostic_variable(self%id_fYG3c,'fYG3c','mg C/m^2/d',  'respiration',                                 domain=domain_bottom,source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_fYKIn,'fYKIn','mmol N/m^2/d','release of dissolved inorganic nitrogen',     domain=domain_bottom,source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_fYK1p,'fYK1p','mmol P/m^2/d','release of dissolved inorganic phosphorus',   domain=domain_bottom,source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_fYQPc,'fYQPc','mg C/m^2/d',  'production of particulate organic carbon',    domain=domain_bottom,source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_fYQPn,'fYQPn','mmol N/m^2/d','production of particulate organic nitrogen',  domain=domain_bottom,source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_fYQPp,'fYQPp','mmol P/m^2/d','production of particulate organic phosphorus',domain=domain_bottom,source=source_do_bottom)

   end subroutine initialize

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)

      class (type_ersem_benthic_fauna),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: c,cP,O2o,foodP,Dm
      real(rk) :: eT,eO,eC,ETW,crowding_c,x
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

      _GET_HORIZONTAL_WITH_BACKGROUND_(self%id_c,c)
      _GET_HORIZONTAL_(self%id_c,cP)

      _GET_HORIZONTAL_(self%id_Dm,Dm)
      _GET_(self%id_O2o,O2o)

      _GET_(self%id_ETW,ETW)

      ! Temperature limitation factor
      eT = max(0.0_rk,self%q10**((ETW-10._rk)/10._rk) - self%q10**((ETW-32._rk)/3._rk))

      ! Oxygen limitation factor
      if (O2o>self%rlO2) then
         eO = (O2o-self%rlO2)**3/((O2o-self%rlO2)**3+(self%hO2-self%rlO2)**3)
      else
         eO = 0.0_rk
      end if

      ! Limitation factor describing a decrease in feeding rate due to overcrowding.
      ! To disable any effect of overcrowding on feeding, set parameter xcl to a very large value.
      ! This can for instance be done to disable the overcrowding effect for meiobenthos/Y4, as in SSB-ERSEM.
      crowding_c = c - self%xcl
      if (crowding_c>0._rk) then
         x = crowding_c * crowding_c/(crowding_c+self%xcs)
         eC = 1._rk - x/(x+self%xch)
      else
         eC = 1._rk
      end if

      ! Calculate maximum ingestion rate (mg C/m2/d).
      ! This incorporates temperature, oxygen, crowding limitation factors.
      rate = self%su * c * eT * eO * eC

      ! Get food concentrations: benthic and pelagic!
      do ifood=1,self%nfood
         if (self%food(ifood)%ispel) then
            _GET_(self%food(ifood)%id_c_pel,foodcP(ifood))
            _GET_(self%food(ifood)%id_n_pel,foodnP(ifood))
            _GET_(self%food(ifood)%id_p_pel,foodpP(ifood))
            _GET_(self%food(ifood)%id_s_pel,foodsP(ifood))

            foodcP_an(ifood) = 0.0_rk
         else
            _GET_HORIZONTAL_(self%food(ifood)%id_c,foodcP(ifood))
            _GET_HORIZONTAL_(self%food(ifood)%id_n,foodnP(ifood))
            _GET_HORIZONTAL_(self%food(ifood)%id_p,foodpP(ifood))
            _GET_HORIZONTAL_(self%food(ifood)%id_s,foodsP(ifood))

            _GET_HORIZONTAL_(self%food(ifood)%id_c_an,foodcP_an(ifood))
         end if
      end do

      ! Carbon content of food has units mmol (due to units of standard_variables%total_carbon); convert to mg
      foodcP = foodcP*CMass

      ! In case of pelagic food source, its availability can be limited by a near-bottom fraction available to organism.
      ! In case of benthic food source, it is possible to limit availability of food source by depth of aerobic layer,
      ! e.g. limiting aerobic bacteria as a food source for suspension-feeders by ratio of habitat depth to aerobic layer depth.
      ! Food sources limited in such a way must be specified using logical food{n}_ll in fabm.yaml file.
      do ifood=1,self%nfood
         if (self%food(ifood)%ispel) then
            prefcorr(ifood) = self%food(ifood)%pu * self%dwat
         else
            if (self%food(ifood)%ll) then
               prefcorr(ifood) = self%food(ifood)%pu * min(1._rk,self%dQ6/Dm)
            else
               prefcorr(ifood) = self%food(ifood)%pu
            end if
         end if
      end do

      ! Compute effective preferences for individual food sources: "feed".
      ! Weighting factors for original preferences increase hyperbolically from 0 at low food density to 1 at high food density.
      feed = prefcorr * (prefcorr * foodcP/(prefcorr * foodcP + self%lu))

      ! Compute specific ingestion of the different food sources: "sflux" (1/d).
      ! These are based on a Michaelis-Menten/Type II functional response with dynamic preferences "feed".
      foodsum = sum(feed * foodcP)
      if (foodsum>0._rk) then
         sflux = rate * feed / (foodsum + self%hu)
      else
         sflux = 0._rk
      end if

      do ifood=1,self%nfood
         _SET_HORIZONTAL_DIAGNOSTIC_(self%food(ifood)%id_fc,sflux(ifood)*foodcP(ifood))
         _SET_HORIZONTAL_DIAGNOSTIC_(self%food(ifood)%id_fn,sflux(ifood)*foodnP(ifood))
         _SET_HORIZONTAL_DIAGNOSTIC_(self%food(ifood)%id_fp,sflux(ifood)*foodpP(ifood))
         _SET_HORIZONTAL_DIAGNOSTIC_(self%food(ifood)%id_fs,sflux(ifood)*foodsP(ifood))
      end do

      ! Ingestion (matter/m2/d) per food source.
      grossfluxc = sflux * foodcP
      grossfluxn = sflux * foodnP
      grossfluxp = sflux * foodpP
      grossfluxs = sflux * foodsP

      ! Compute assimilation (matter/m2/d) per food source from ingestion and assimilation inefficiency.
      netfluxc = grossfluxc * (1._rk -            self%food%pue)
      netfluxn = grossfluxn * (1._rk - self%pudil*self%food%pue)
      netfluxp = grossfluxp * (1._rk - self%pudil*self%food%pue)

      ! Based on specific ingestion of each food source, decrease all state variables of that food source.
      do ifood=1,self%nfood
         do istate=1,size(self%food(ifood)%id_source%state)
            _GET_(self%food(ifood)%id_source%state(istate),foodP)
            _SET_BOTTOM_EXCHANGE_(self%food(ifood)%id_loss_source%state(istate),-sflux(ifood)*foodP)
         end do
         do istate=1,size(self%food(ifood)%id_source%bottom_state)
            _GET_HORIZONTAL_(self%food(ifood)%id_source%bottom_state(istate),foodP)
            _SET_BOTTOM_ODE_(self%food(ifood)%id_loss_source%bottom_state(istate),-sflux(ifood)*foodP)
         end do
      end do

      ! Total gross uptake (ingestion) and net uptake (assimilation), summed over all food sources (mg C/m2/d).
      fBTYc = sum(grossfluxc)
      nfBTYc = sum(netfluxc)

      ! Store net carbon, nitrogen and phosphorus fluxes associated with food uptake for later use,
      ! and send difference between ingestion and assimilation to faeces.
      SYc = nfBTYc
      SYn = sum(netfluxn)
      SYp = sum(netfluxp)
      _SET_BOTTOM_ODE_(self%id_Q6c,fBTYc - nfBTYc)
      _SET_BOTTOM_ODE_(self%id_Q6n,sum(grossfluxn) - sum(netfluxn))
      _SET_BOTTOM_ODE_(self%id_Q6p,sum(grossfluxp) - sum(netfluxp))
      _SET_BOTTOM_ODE_(self%id_Q6s,sum(grossfluxs))

      ! Compute contribution to bioturbation and bioirrigation from total carbon ingestion.
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_biotur, self%ptur * fBTYc)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_bioirr, self%pirr * fBTYc)

      ! Respiration fluxes = basal respiration (proportional to biomass)+ activity respiration (proportional to carbon assimilation)
      fYG3c = self%sr * cP * eT + self%pur * nfBTYc

      ! Store carbon flux resulting from respiration for later use (note: respiration does not affect nitrogen, phosphorus).
      ! Also account for its production of benthic CO2 and consumption of benthic oxygen.
      SYc = SYc - fYG3c
      _SET_BOTTOM_ODE_(self%id_G3c, fYG3c/CMass)
      _SET_BOTTOM_ODE_(self%id_G2o,-fYG3c/CMass)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fYG3c,fYG3c)

      ! Specific mortality (1/d): background mortality + mortality due to oxygen limitation + mortality due to cold.
      mort_act  = self%sd * eT
      mort_ox   = self%sdmO2 * (1._rk - eO) * eT
      mort_cold = self%sdc * exp(ETW * self%xdc)

      ! Total specific mortality (1/d)
      mortflux = mort_act + mort_ox + mort_cold

      ! Apply mortality to biomass, and send dead matter to particulate organic carbon pool.
      _SET_BOTTOM_ODE_(self%id_c, -mortflux*cP)
      _SET_BOTTOM_ODE_(self%id_Q6c,mortflux*cP)
      _SET_BOTTOM_ODE_(self%id_Q6n,mortflux*cP*self%qnc)
      _SET_BOTTOM_ODE_(self%id_Q6p,mortflux*cP*self%qpc)

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

      p_an = sum(foodcP_an/(foodcP+1.e-15_rk)*grossfluxc)/max(fBTYc,1.e-8_rk)

      ! Send excess nutrients to phosphate, ammonium pools.
      _SET_BOTTOM_ODE_(self%id_K4n,(1._rk-p_an)*excess_n)
      _SET_BOTTOM_ODE_(self%id_K4n2,p_an       *excess_n)
      _SET_BOTTOM_ODE_(self%id_K1p,(1._rk-p_an)*excess_p)
      _SET_BOTTOM_ODE_(self%id_K1p2,p_an       *excess_p)

      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fYKIn,excess_n)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fYK1p,excess_p)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fYQPc,mortflux*cP+excess_c)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fYQPn,mortflux*cP*self%qnc)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fYQPp,mortflux*cP*self%qpc)

      if (.not.legacy_ersem_compatibility) then
         ! Alkalinity contributions: +1 for NH4, -1 for PO4
         _SET_BOTTOM_ODE_(self%id_benTA,(1._rk-p_an)*(excess_n-excess_p))
         _SET_BOTTOM_ODE_(self%id_benTA2,p_an       *(excess_n-excess_p))
      end if

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module
