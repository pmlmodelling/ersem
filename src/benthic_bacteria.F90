#include "fabm_driver.h"

module ersem_benthic_bacteria

   use fabm_types

   use fabm_particle
   use ersem_shared
   use ersem_benthic_base

   implicit none

   private

   type,extends(type_ersem_benthic_base),public :: type_ersem_benthic_bacteria
      type (type_bottom_state_variable_id) :: id_K1p,id_K4n
      type (type_bottom_state_variable_id) :: id_Q1c,id_Q1n,id_Q1p
      type (type_bottom_state_variable_id) :: id_Q6c,id_Q6n,id_Q6p
      type (type_bottom_state_variable_id) :: id_Q7c,id_Q7n,id_Q7p
      type (type_bottom_state_variable_id) :: id_G2o,id_G3c,id_benTA
      type (type_horizontal_diagnostic_variable_id) :: id_fHG3c

      type (type_horizontal_dependency_id) :: id_Dm
      type (type_dependency_id) :: id_ETW

      type (type_model_id) :: id_Q1,id_Q6,id_Q7
      
      real(rk) :: qnc,qpc
      real(rk) :: q10
      real(rk) :: dd
      real(rk) :: suQ7
      real(rk) :: suQ6f
      real(rk) :: suQ6s
      real(rk) :: suQ1
      real(rk) :: puinc
      real(rk) :: pue7,pue6
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
      real(rk) :: c0

      ! Initialize base model type (this will also set the internal time unit to per day)
      call self%initialize_ersem_benthic_base()

      ! Register parameters
      call self%get_parameter(self%qnc,  'qnc',  'mmol N/mg C','nitrogen to carbon ratio')
      call self%get_parameter(self%qpc,  'qpc',  'mmol P/mg C','phosphorus to carbon ratio')
      call self%get_parameter(c0,        'c0',   'mg C/m^2',   'background concentration',default=0.0_rk)
      call self%get_parameter(self%q10,  'q10',  '-',          'Q_10 temperature coefficient')
      call self%get_parameter(self%dd,   'dd',   'm',          'Michaelis-Menten constant for oxygen limitation through layer thickness')
      call self%get_parameter(self%suQ7, 'suQ7', 'm^2/mg C/d', 'affinity for refractory matter')
      call self%get_parameter(self%suQ6s,'suQ6s','m^2/mg C/d', 'base affinity for particulate organic matter')
      call self%get_parameter(self%suQ6f,'suQ6f','m^2/mg C/d', 'additional affinity for particulate organic matter subject to nutrient limitation')
      call self%get_parameter(self%suQ1, 'suQ1', 'm^2/mg C/d', 'affinity for dissolved organic matter')
      call self%get_parameter(self%puinc,'puinc','-',          'preference for nutrients relative to carbon in particulate organic matter')
      call self%get_parameter(self%pue6, 'pue6', '-',          'fraction of consumed particulate organic matter that is excreted')
      call self%get_parameter(self%pue7, 'pue7', '-',          'fraction of consumed refractory matter that is excreted')
      call self%get_parameter(self%pur,  'pur',  '-',          'fraction of consumed carbon that is respired')
      call self%get_parameter(self%sr,   'sr',   '1/d',        'specific rest respiration')
      call self%get_parameter(self%pdQ1, 'pdQ1', '-',          'fraction of dying matter that is dissolved')
      call self%get_parameter(self%sd,   'sd',   '1/d',        'specific maximum mortality related to oxygen limitation')

      call self%add_constituent('c',0.0_rk,c0,qn=self%qnc,qp=self%qpc)

      ! Environmental dependencies
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      call self%register_dependency(self%id_Dm,'Dm','m','depth interval available to bacteria')

      ! Dependencies on state variables of external modules.
      call self%register_state_dependency(self%id_K4n,'K4n','mmol N/m^2','ammonium')
      call self%register_state_dependency(self%id_K1p,'K1p','mmol N/m^2','phosphate')
      call self%register_state_dependency(self%id_G2o,'G2o','mmol O_2/m^2','oxygen')
      call self%register_state_dependency(self%id_G3c,'G3c','mmol C/m^2','dissolved inorganic carbon')
      if (.not.legacy_ersem_compatibility) call self%register_state_dependency(self%id_benTA,'benTA','mEq/m^2','benthic alkalinity')
      call self%register_state_dependency(self%id_Q1c,'Q1c','mg C/m^2',  'dissolved organic carbon')
      call self%register_state_dependency(self%id_Q1n,'Q1n','mmol N/m^2','dissolved organic nitrogen')
      call self%register_state_dependency(self%id_Q1p,'Q1p','mmol P/m^2','dissolved organic phosphorus')
      call self%register_state_dependency(self%id_Q6c,'Q6c','mg C/m^2',  'particulate organic carbon')
      call self%register_state_dependency(self%id_Q6n,'Q6n','mmol N/m^2','particulate organic nitrogen')
      call self%register_state_dependency(self%id_Q6p,'Q6p','mmol P/m^2','particulate organic phosphorus')
      call self%register_state_dependency(self%id_Q7c,'Q7c','mg C/m^2',  'refractory particulate organic carbon')
      call self%register_state_dependency(self%id_Q7n,'Q7n','mmol N/m^2','refractory particulate organic nitrogen')
      call self%register_state_dependency(self%id_Q7p,'Q7p','mmol P/m^2','refractory particulate organic phosphorus')

      ! Allow the user to hook up substrate constituents by bulk coupling to entire models.
      ! This allows one to couple e.g. Q1 as a whole, rather than Q1c, Q1n, Q1p.
      call self%register_model_dependency(self%id_Q1,'Q1')
      call self%register_model_dependency(self%id_Q6,'Q6')
      call self%register_model_dependency(self%id_Q7,'Q7')
      call self%request_coupling_to_model(self%id_Q1c,self%id_Q1,'c')
      call self%request_coupling_to_model(self%id_Q1n,self%id_Q1,'n')
      call self%request_coupling_to_model(self%id_Q1p,self%id_Q1,'p')
      call self%request_coupling_to_model(self%id_Q6c,self%id_Q6,'c')
      call self%request_coupling_to_model(self%id_Q6n,self%id_Q6,'n')
      call self%request_coupling_to_model(self%id_Q6p,self%id_Q6,'p')
      call self%request_coupling_to_model(self%id_Q7c,self%id_Q7,'c')
      call self%request_coupling_to_model(self%id_Q7n,self%id_Q7,'n')
      call self%request_coupling_to_model(self%id_Q7p,self%id_Q7,'p')

      call self%register_diagnostic_variable(self%id_fHG3c,'fHG3c','mg C/m^2/d','respiration',output=output_time_step_averaged,domain=domain_bottom,source=source_do_bottom)

   end subroutine

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)

      class (type_ersem_benthic_bacteria),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
     
      real(rk) :: Hc,HcP,fHc,fHn,fHp
      real(rk) :: Dm,excess_c,excess_n,excess_p
      real(rk) :: ETW,eT,eOx,eN,Limit
      real(rk) :: K4a,K1a
      real(rk) :: Q1cP,Q1nP,Q1pP
      real(rk) :: AQ6c,AQ6n,AQ6p
      real(rk) :: AQ7c,AQ7n,AQ7p
      real(rk) :: sfQ7H,sfQ6H,sfQ1H,fQ7Hc,fQ6Hc,fQ1Hc,fQIHc,fQ7Hn,fQ7Hp,fQ6Hn,fQ6Hp
      real(rk) :: fQ1Hn,fQ1Hp,fK1Hp,fHG3c,fK4Hn,sfHQ1,sfHQI,sfHQ6

      _HORIZONTAL_LOOP_BEGIN_

         ! Retrieve bacterial carbon
         _GET_HORIZONTAL_(self%id_c,HcP)
         _GET_HORIZONTAL_WITH_BACKGROUND_(self%id_c,Hc)

         ! Retrieve environmental conditions
         _GET_(self%id_ETW,ETW)            ! temperature (of lowermost pelagic layer)
         _GET_HORIZONTAL_(self%id_Dm,Dm)   ! depth of the sediment layer that represents bacterial habitat

         ! Retrieve organic carbon substrates
         _GET_HORIZONTAL_(self%id_Q1c,Q1cP)
         _GET_HORIZONTAL_(self%id_Q1n,Q1nP)
         _GET_HORIZONTAL_(self%id_Q1p,Q1pP)
         _GET_HORIZONTAL_(self%id_Q6c,AQ6c)
         _GET_HORIZONTAL_(self%id_Q6n,AQ6n)
         _GET_HORIZONTAL_(self%id_Q6p,AQ6p)
         _GET_HORIZONTAL_(self%id_Q7c,AQ7c)
         _GET_HORIZONTAL_(self%id_Q7n,AQ7n)
         _GET_HORIZONTAL_(self%id_Q7p,AQ7p)

         ! Retrieve dissolved nutrients (phosphate, ammonium)
         _GET_HORIZONTAL_(self%id_K1p,K1a)
         _GET_HORIZONTAL_(self%id_K4n,K4a)

         ! Limitation by temperature, oxygen, nutrient availability.
         eT  = self%q10**((ETW-10._rk)/10._rk) - self%q10**((ETW-32._rk)/3._rk)
         eOx = Dm/(self%dd+Dm)
         if (AQ6c<=0.0_rk) then
            ! Avoid division by zero: no carbon in substrate means no nutrient limitation
            eN = 1.0_rk
         else
            eN = min(1._rk,max(0._rk,AQ6n/(self%qnc*AQ6c))) * min(1._rk,max(0._rk,AQ6p/(self%qpc*AQ6c)))
         end if

         ! Combined limitation factor
         Limit = eT * eOX * Hc

         ! Determine effective uptake rates
         sfQ7H = self%suQ7  * Limit
         sfQ6H = self%suQ6s * Limit + self%suQ6f * Limit * eN
         sfQ1H = self%suQ1  * Limit

         ! Uptake of substrate carbon
         fQ7Hc = sfQ7H * AQ7c
         fQ6Hc = sfQ6H * AQ6c
         fQ1Hc = sfQ1H * Q1cP
         fQIHc = fQ7Hc + fQ6Hc + fQ1Hc

         ! Uptake of substrate nitrogen and phosphorus
         ! (note on Q6: affinity for n and p differs from affinity for c if puinc/=1)
         fQ7Hn = sfQ7H * AQ7n
         fQ7Hp = sfQ7H * AQ7p
         fQ6Hn = sfQ6H * AQ6n * self%puinc
         fQ6Hp = sfQ6H * AQ6p * self%puinc
         fQ1Hn = sfQ1H * Q1nP
         fQ1Hp = sfQ1H * Q1pP

         ! (JB: code below seems to want to infer nutrient requirement/flux,
         ! but it makes no sense since K?a and fK?H? are summed while they have different units)
         ! The "if" clauses below protect against the pathologic case where the nutrient requirement
         ! and the nutrient concentation (fK4Hn and K4a, or fK1Hp and K1a) are both 0.
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

         ! Mortality (distribute over dissolved organics Q1 and particulate organics Q6)
         sfHQI = self%sd * (1._rk - eOx)
         sfHQ6 = sfHQI * (1._rk - self%pdQ1)
         sfHQ1 = sfHQI * self%pdQ1

         _SET_BOTTOM_ODE_(self%id_Q6c, sfHQ6 * HcP)
         _SET_BOTTOM_ODE_(self%id_Q6n, sfHQ6 * HcP * self%qnc)
         _SET_BOTTOM_ODE_(self%id_Q6p, sfHQ6 * HcP * self%qpc)

         _SET_BOTTOM_ODE_(self%id_Q1c, sfHQ1 * HcP)
         _SET_BOTTOM_ODE_(self%id_Q1n, sfHQ1 * HcP * self%qnc)
         _SET_BOTTOM_ODE_(self%id_Q1p, sfHQ1 * HcP * self%qpc)

         _SET_BOTTOM_ODE_(self%id_c, -sfHQI * HcP)

         ! Non-balanced fluxes of carbon, nitrogen and phosphorus in biomass.
         ! That includes all fluxes (uptake of organic matter and nutrients, respiration) except mortality,
         ! which is stoichiometrically balanced already.
         fHc = fQ7Hc * (1._rk-self%pue7) + fQ6Hc * (1._rk-self%pue6) + fQ1Hc - fHG3c
         fHn = fQ7Hn * (1._rk-self%pue7) + fQ6Hn * (1._rk-self%pue6) + fQ1Hn + fK4Hn
         fHp = fQ7Hp * (1._rk-self%pue7) + fQ6Hp * (1._rk-self%pue6) + fQ1Hp + fK1Hp

         ! Compute excess carbon flux, given that the maximum realizable carbon flux needs to be balanced by corresponding nitrogen and phosphorus fluxes to maintain constant stoichiometry.
         excess_c = max(max(fHc - fHp/self%qpc, fHc - fHn/self%qnc), 0.0_rk)
         fHc = fHc - excess_c
         _SET_BOTTOM_ODE_(self%id_c,fHc)

         ! Compute excess nitrogen and phosphorus fluxes, based on final carbon flux. Excess nutrient will be exudated to preserve constant stoichiometry of biomass.
         excess_n = max(fHn - fHc*self%qnc,0.0_rk)
         excess_p = max(fHp - fHc*self%qpc,0.0_rk)

         _SET_BOTTOM_ODE_(self%id_Q6c,-fQ6Hc + excess_c)
         _SET_BOTTOM_ODE_(self%id_Q6n,-fQ6Hn)
         _SET_BOTTOM_ODE_(self%id_Q6p,-fQ6Hp)

         _SET_BOTTOM_ODE_(self%id_K4n,-fK4Hn + excess_n)
         _SET_BOTTOM_ODE_(self%id_K1p,-fK1Hp + excess_p)

         ! Alkalinity contributions: +1 for NH4, -1 for PO4
         if (.not.legacy_ersem_compatibility) _SET_BOTTOM_ODE_(self%id_benTA,-fK4Hn + excess_n + fK1Hp -excess_p)

         _SET_BOTTOM_ODE_(self%id_Q7c,-fQ7Hc)
         _SET_BOTTOM_ODE_(self%id_Q7n,-fQ7Hn)
         _SET_BOTTOM_ODE_(self%id_Q7p,-fQ7Hp)

         _SET_BOTTOM_ODE_(self%id_Q1c,-fQ1Hc + fQ7Hc*self%pue7 + fQ6Hc*self%pue6)
         _SET_BOTTOM_ODE_(self%id_Q1n,-fQ1Hn + fQ7Hn*self%pue7 + fQ6Hn*self%pue6)
         _SET_BOTTOM_ODE_(self%id_Q1p,-fQ1Hp + fQ7Hp*self%pue7 + fQ6Hp*self%pue6)

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module
