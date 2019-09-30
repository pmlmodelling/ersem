#include "fabm_driver.h"

module ersem_benthic_bacteria

   use fabm_types

   use fabm_particle
   use ersem_shared
   use ersem_benthic_base

   implicit none

   private

   type type_food
      type (type_bottom_state_variable_id)          :: id_c, id_n, id_p
      type (type_horizontal_diagnostic_variable_id) :: id_fc,id_fn,id_fp
      real(rk) :: su
      real(rk) :: suf
      real(rk) :: puinc
      real(rk) :: pue
   end type

   type,extends(type_ersem_benthic_base),public :: type_ersem_benthic_bacteria
      type (type_bottom_state_variable_id) :: id_K1p,id_K4n
      type (type_bottom_state_variable_id) :: id_Q1c,id_Q1n,id_Q1p
      type (type_bottom_state_variable_id) :: id_Q6c,id_Q6n,id_Q6p
      type (type_bottom_state_variable_id) :: id_G2o,id_G3c,id_benTA
      type (type_food), allocatable :: food(:)
      type (type_horizontal_diagnostic_variable_id) :: id_fHG3c
      type (type_horizontal_diagnostic_variable_id) :: id_fHKIn,id_fHK1p
      type (type_horizontal_diagnostic_variable_id) :: id_fHQ1c,id_fHQPc,id_fHQ1n,id_fHQPn,id_fHQ1p,id_fHQPp
      type (type_horizontal_dependency_id) :: id_Dm

      integer  :: nfood
      real(rk) :: qnc,qpc
      real(rk) :: dd
      real(rk) :: pur,sr
      real(rk) :: pdQ1
      real(rk) :: sd
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self,configunit)
      class (type_ersem_benthic_bacteria),intent(inout),target :: self
      integer,                            intent(in)           :: configunit

      real(rk)         :: c0
      integer          :: ifood
      character(len=6) :: strindex

      ! Initialize base model type (this will also set the internal time unit to per day)
      call self%initialize_ersem_benthic_base()

      ! Get first set of parameters
      call self%get_parameter(self%qnc,  'qnc',  'mmol N/mg C','nitrogen to carbon ratio')
      call self%get_parameter(self%qpc,  'qpc',  'mmol P/mg C','phosphorus to carbon ratio')
      call self%get_parameter(c0,        'c0',   'mg C/m^2',   'background concentration',default=0.0_rk)
      call self%get_parameter(self%q10,  'q10',  '-',          'Q_10 temperature coefficient')
      call self%get_parameter(self%dd,   'dd',   'm',          'Michaelis-Menten constant for oxygen limitation through layer thickness')

      call self%add_constituent('c',0.0_rk,c0,qn=self%qnc,qp=self%qpc)

      ! Environmental dependencies
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      call self%register_dependency(self%id_Dm,'Dm','m','depth interval available to bacteria')

      call self%get_parameter(self%nfood,'nfood','','number of food sources')
      allocate(self%food(self%nfood))
      do ifood=1,self%nfood
         write (strindex,'(i0)') ifood

         ! Get substrate-specific parameters.
         call self%get_parameter(self%food(ifood)%su,   'su'//trim(strindex),   'm^2/mg C/d', 'base affinity for food source '//trim(strindex))
         call self%get_parameter(self%food(ifood)%suf,  'suf'//trim(strindex),  'm^2/mg C/d', 'additional affinity for food source '//trim(strindex)//' subject to nutrient limitation', default=0.0_rk)
         call self%get_parameter(self%food(ifood)%puinc,'puinc'//trim(strindex),'-',          'affinity for nutrients relative to carbon in food source '//trim(strindex), default=1.0_rk)
         call self%get_parameter(self%food(ifood)%pue,  'pue'//trim(strindex),  '-',          'fraction of consumed food source '//trim(strindex)//' that is excreted', default=0.0_rk)

         ! Register state variables for substrate constituents.
         call self%register_state_dependency(self%food(ifood)%id_c, 'food'//trim(strindex)//'c', 'mmol C/m^2', 'carbon in food source '//trim(strindex))
         call self%register_state_dependency(self%food(ifood)%id_n, 'food'//trim(strindex)//'n', 'mmol N/m^2', 'nitrogen in food source '//trim(strindex))
         call self%register_state_dependency(self%food(ifood)%id_p, 'food'//trim(strindex)//'p', 'mmol P/m^2', 'phosphorus in food source '//trim(strindex))

         ! By default, couple substrate constituents to a "food?" module.
         call self%request_coupling_to_model(self%food(ifood)%id_c, 'food'//trim(strindex), standard_variables%total_carbon)
         call self%request_coupling_to_model(self%food(ifood)%id_n, 'food'//trim(strindex), standard_variables%total_nitrogen)
         call self%request_coupling_to_model(self%food(ifood)%id_p, 'food'//trim(strindex), standard_variables%total_phosphorus)

         ! Register diagnostics for uptake of the substrate constituents.
         call self%register_diagnostic_variable(self%food(ifood)%id_fc, 'fc'//trim(strindex), 'mg C/m^2/d',   'uptake of carbon in food source '//trim(strindex),    source=source_do_bottom)
         call self%register_diagnostic_variable(self%food(ifood)%id_fn, 'fn'//trim(strindex), 'mmol N/m^2/d', 'uptake of nitrogen in food source '//trim(strindex),  source=source_do_bottom)
         call self%register_diagnostic_variable(self%food(ifood)%id_fp, 'fp'//trim(strindex), 'mmol P/m^2/d', 'uptake of phosphorus in food source '//trim(strindex),source=source_do_bottom)
      end do

      ! Get remaining parameters.
      call self%get_parameter(self%pur,  'pur',  '-',   'fraction of consumed carbon that is respired')
      call self%get_parameter(self%sr,   'sr',   '1/d', 'specific rest respiration')
      call self%get_parameter(self%pdQ1, 'pdQ1', '-',   'fraction of dying matter that is dissolved')
      call self%get_parameter(self%sd,   'sd',   '1/d', 'specific maximum mortality related to oxygen limitation')

      ! Dependencies on state variables of external modules.
      call self%register_state_dependency(self%id_K4n,'K4n','mmol N/m^2','ammonium')
      call self%register_state_dependency(self%id_K1p,'K1p','mmol N/m^2','phosphate')
      call self%register_state_dependency(self%id_G2o,'G2o','mmol O_2/m^2','oxygen')
      call self%register_state_dependency(self%id_G3c,'G3c','mmol C/m^2','dissolved inorganic carbon')
      if (.not.legacy_ersem_compatibility) call self%register_state_dependency(self%id_benTA,'benTA','mEq/m^2','benthic alkalinity')
      call self%register_state_dependency(self%id_Q1c,'Q1c','mmol C/m^2','dissolved organic carbon')
      call self%register_state_dependency(self%id_Q1n,'Q1n','mmol N/m^2','dissolved organic nitrogen')
      call self%register_state_dependency(self%id_Q1p,'Q1p','mmol P/m^2','dissolved organic phosphorus')
      call self%register_state_dependency(self%id_Q6c,'Q6c','mmol C/m^2','particulate organic carbon')
      call self%register_state_dependency(self%id_Q6n,'Q6n','mmol N/m^2','particulate organic nitrogen')
      call self%register_state_dependency(self%id_Q6p,'Q6p','mmol P/m^2','particulate organic phosphorus')

      ! By default, couple all POM/DOM constituents wholesale to an entire module (Q1, Q6).
      ! This allows one to couple e.g. Q1 as a whole, rather than Q1c, Q1n, Q1p.
      call self%request_coupling_to_model(self%id_Q1c,'Q1',standard_variables%total_carbon)
      call self%request_coupling_to_model(self%id_Q1n,'Q1',standard_variables%total_nitrogen)
      call self%request_coupling_to_model(self%id_Q1p,'Q1',standard_variables%total_phosphorus)
      call self%request_coupling_to_model(self%id_Q6c,'Q6',standard_variables%total_carbon)
      call self%request_coupling_to_model(self%id_Q6n,'Q6',standard_variables%total_nitrogen)
      call self%request_coupling_to_model(self%id_Q6p,'Q6',standard_variables%total_phosphorus)

      call self%register_diagnostic_variable(self%id_fHG3c,'fHG3c','mg C/m^2/d',  'respiration',                                 domain=domain_bottom,source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_fHKIn,'fHKIn','mmol N/m^2/d','release of dissolved inorganic nitrogen',     domain=domain_bottom,source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_fHK1p,'fHK1p','mmol P/m^2/d','release of dissolved inorganic phosphorus',   domain=domain_bottom,source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_fHQ1c,'fHQ1c','mg C/m^2/d',  'production of dissolved organic carbon',      domain=domain_bottom,source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_fHQ1n,'fHQ1n','mmol N/m^2/d','production of dissolved organic nitrogen',    domain=domain_bottom,source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_fHQ1p,'fHQ1p','mmol P/m^2/d','production of dissolved organic phosphorus',  domain=domain_bottom,source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_fHQPc,'fHQPc','mg C/m^2/d',  'production of particulate organic carbon',    domain=domain_bottom,source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_fHQPn,'fHQPn','mmol N/m^2/d','production of particulate organic nitrogen',  domain=domain_bottom,source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_fHQPp,'fHQPp','mmol P/m^2/d','production of particulate organic phosphorus',domain=domain_bottom,source=source_do_bottom)

   end subroutine

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)

      class (type_ersem_benthic_bacteria),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      integer :: ifood
      real(rk) :: Hc,HcP,fHc,fHn,fHp
      real(rk) :: Dm,excess_c,excess_n,excess_p
      real(rk) :: ETW,eT,eOx,eN
      real(rk) :: K4a,K1a
      real(rk),dimension(self%nfood) :: Qc,Qn,Qp,sfQ,fQc,fQn,fQp
      real(rk) :: fQIHc
      real(rk) :: fK1Hp,fHG3c,fK4Hn,sfHQ1,sfHQI,sfHQ6

      _HORIZONTAL_LOOP_BEGIN_

         ! Retrieve bacterial carbon
         _GET_HORIZONTAL_(self%id_c,HcP)
         _GET_HORIZONTAL_WITH_BACKGROUND_(self%id_c,Hc)

         ! Retrieve environmental conditions
         _GET_(self%id_ETW,ETW)            ! temperature (of lowermost pelagic layer)
         _GET_HORIZONTAL_(self%id_Dm,Dm)   ! depth of the sediment layer that represents bacterial habitat

         ! Retrieve carbon, nitrogen, phosphorus in substrates.
         ! Carbon (standard_variables%total_carbon) is returned in mmol/m2, so convert to mg/m2.
         do ifood=1,self%nfood
            _GET_HORIZONTAL_(self%food(ifood)%id_c,Qc(ifood))
            _GET_HORIZONTAL_(self%food(ifood)%id_n,Qn(ifood))
            _GET_HORIZONTAL_(self%food(ifood)%id_p,Qp(ifood))
         end do
         Qc = Qc*CMass

         ! Retrieve dissolved nutrients (phosphate, ammonium)
         _GET_HORIZONTAL_(self%id_K1p,K1a)
         _GET_HORIZONTAL_(self%id_K4n,K4a)

         ! Limitation by temperature and height of habitat layer.
         eT  = max(0.0_rk,self%q10**((ETW-10._rk)/10._rk) - self%q10**((ETW-32._rk)/3._rk))
         eOx = Dm/(self%dd+Dm)

         ! Effective uptake rate (sfQ, 1/d) per substrate
         do ifood=1,self%nfood
            if (self%food(ifood)%suf==0.0_rk .or. Qn(ifood)==0.0_rk .or. Qp(ifood)==0.0_rk) then
               ! Avoid division by zero when both carbon and nitrogen/phosphorus in substrate are zero:
               ! no nutrients (or no separate nutrient-limited uptake) mean complete nutrient limitation.
               eN = 0.0_rk
            else
               eN = min(1._rk,max(0._rk,Qn(ifood)/(self%qnc*Qc(ifood)))) * min(1._rk,max(0._rk,Qp(ifood)/(self%qpc*Qc(ifood))))
            end if
            sfQ(ifood) = (self%food(ifood)%su + self%food(ifood)%suf * eN) * eT * eOX * Hc
         end do

         ! Gross carbon, nitrogen, phosphorus uptake per substrate (mg/m2/d or mmol/m2/d)
         fQc = sfQ * Qc
         fQn = sfQ * Qn * self%food%puinc
         fQp = sfQ * Qp * self%food%puinc

         ! Total gross carbon uptake (mg C/m2/d)
         fQIHc = sum(fQc)

         ! (JB: code below seems to want to infer nutrient requirement/flux,
         ! but it makes no sense since K?a and fK?H? are summed while they have different units)
         ! The "if" clauses below protect against the pathologic case where the nutrient requirement
         ! and the nutrient concentration (fK4Hn and K4a, or fK1Hp and K1a) are both 0.
         fK4Hn = fQIHc * self%qnc
         if (fK4Hn>0) fK4Hn = fK4Hn * K4a/(K4a+fK4Hn)
         fK1Hp = fQIHc * self%qpc
         if (fK1Hp>0) fK1Hp = fK1Hp * K1a/(K1a+fK1Hp)

         ! Respiration (reduction in bacterial carbon and dissolved oxygen, increase in
         ! dissolved inorganic carbon, ammonium, phosphate)
         fHG3c = self%pur * fQIHc + self%sr * HcP * eT
         _SET_BOTTOM_ODE_(self%id_G2o,-fHG3c/CMass)  ! oxygen or reduction equivalent
         _SET_BOTTOM_ODE_(self%id_G3c, fHG3c/CMass)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fHG3c,fHG3c)

         ! Mortality (partition over dissolved and particulate organic pools according to coefficient pdQ1)
         sfHQI = self%sd * (1._rk - eOx)
         sfHQ6 = sfHQI * (1._rk - self%pdQ1)
         sfHQ1 = sfHQI * self%pdQ1

         _SET_BOTTOM_ODE_(self%id_Q6c, sfHQ6 * HcP/CMass)
         _SET_BOTTOM_ODE_(self%id_Q6n, sfHQ6 * HcP * self%qnc)
         _SET_BOTTOM_ODE_(self%id_Q6p, sfHQ6 * HcP * self%qpc)

         _SET_BOTTOM_ODE_(self%id_Q1c, sfHQ1 * HcP/CMass)
         _SET_BOTTOM_ODE_(self%id_Q1n, sfHQ1 * HcP * self%qnc)
         _SET_BOTTOM_ODE_(self%id_Q1p, sfHQ1 * HcP * self%qpc)

         _SET_BOTTOM_ODE_(self%id_c, -sfHQI * HcP)

         ! Non-balanced fluxes of carbon, nitrogen and phosphorus in biomass.
         ! That includes all fluxes (uptake of organic matter and nutrients, respiration) except mortality,
         ! which is stoichiometrically balanced already.
         fHc = sum(fQc*(1._rk-self%food%pue)) - fHG3c
         fHn = sum(fQn*(1._rk-self%food%pue)) + fK4Hn
         fHp = sum(fQp*(1._rk-self%food%pue)) + fK1Hp

         ! Compute excess carbon flux, given that the maximum realizable carbon flux needs to be balanced by corresponding nitrogen and phosphorus fluxes to maintain constant stoichiometry.
         excess_c = max(max(fHc - fHp/self%qpc, fHc - fHn/self%qnc), 0.0_rk)
         fHc = fHc - excess_c
         _SET_BOTTOM_ODE_(self%id_c,fHc)

         ! Compute excess nitrogen and phosphorus fluxes, based on final carbon flux. Excess nutrient will be exudated to preserve constant stoichiometry of biomass.
         excess_n = max(fHn - fHc*self%qnc,0.0_rk)
         excess_p = max(fHp - fHc*self%qpc,0.0_rk)

         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fHKIn,-fK4Hn + excess_n)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fHK1p,-fK1Hp + excess_p)

         do ifood=1,self%nfood
            ! Destruct known substrate constituents (carbon, nitrogen, phosphorus)
            ! Unlike zooplankton and benthic fauna, we do not automatically destroy any other (unknown)
            ! constituents in our food. This ensures that particulate silicate in food is not destroyed by bacteria.
            _SET_BOTTOM_ODE_(self%food(ifood)%id_c,-fQc(ifood)/CMass)
            _SET_BOTTOM_ODE_(self%food(ifood)%id_n,-fQn(ifood))
            _SET_BOTTOM_ODE_(self%food(ifood)%id_p,-fQp(ifood))
            _SET_HORIZONTAL_DIAGNOSTIC_(self%food(ifood)%id_fc,fQc(ifood))
            _SET_HORIZONTAL_DIAGNOSTIC_(self%food(ifood)%id_fn,fQn(ifood))
            _SET_HORIZONTAL_DIAGNOSTIC_(self%food(ifood)%id_fp,fQp(ifood))
         end do

         _SET_BOTTOM_ODE_(self%id_Q6c,excess_c/CMass)
         _SET_BOTTOM_ODE_(self%id_K4n,-fK4Hn + excess_n)
         _SET_BOTTOM_ODE_(self%id_K1p,-fK1Hp + excess_p)

         ! Alkalinity contributions: +1 for NH4, -1 for PO4
         if (.not.legacy_ersem_compatibility) _SET_BOTTOM_ODE_(self%id_benTA,-fK4Hn + excess_n + fK1Hp -excess_p)

         _SET_BOTTOM_ODE_(self%id_Q1c,sum(fQc*self%food%pue)/CMass)
         _SET_BOTTOM_ODE_(self%id_Q1n,sum(fQn*self%food%pue))
         _SET_BOTTOM_ODE_(self%id_Q1p,sum(fQp*self%food%pue))

         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fHQ1c,sum(fQc*self%food%pue) + sfHQ1 * HcP)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fHQ1n,sum(fQn*self%food%pue) + sfHQ1 * HcP*self%qnc)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fHQ1p,sum(fQp*self%food%pue) + sfHQ1 * HcP*self%qpc)

         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fHQPc,sfHQ6 * HcP + excess_c)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fHQPn,sfHQ6 * HcP * self%qnc)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fHQPp,sfHQ6 * HcP * self%qpc)

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module
