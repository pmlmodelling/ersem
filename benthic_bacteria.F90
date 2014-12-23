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
      type (type_bottom_state_variable_id) :: id_G2o,id_G3c

      type (type_horizontal_dependency_id) :: id_Dm
      type (type_dependency_id) :: id_ETW

      type (type_model_id) :: id_Q1,id_Q6,id_Q7
      
      real(rk) :: qnHIcX,qpHIcX
      real(rk) :: q10HX
      real(rk) :: ddHX
      real(rk) :: suQ7HX
      real(rk) :: suQ6fHX
      real(rk) :: suQ6sHX
      real(rk) :: suQ1HX
      real(rk) :: puincHX
      real(rk) :: p_ex7ox,p_ex6ox
      real(rk) :: purHX,srHX
      real(rk) :: pdHQ1X
      real(rk) :: sdHX
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self,configunit)
      class (type_ersem_benthic_bacteria),intent(inout),target :: self
      integer,                                 intent(in)           :: configunit

      ! Set time unit to d-1
      ! This implies that all rates (sink/source terms, vertical velocities) are given in d-1.
      self%dt = 86400._rk

      ! Register parameters
      call self%get_parameter(self%qnHIcX, 'qnc',  'mmol N/mg C','Maximum nitrogen to carbon ratio')
      call self%get_parameter(self%qpHIcX, 'qpc',  'mmol P/mg C','Maximum phosphorus to carbon ratio')
      call self%get_parameter(self%q10HX,  'q10',  '-',          'Q_10 temperature coefficient')
      call self%get_parameter(self%ddHX,   'dd',   '1/m',        'Michaelis-Menten constant for limitation through layer thickness')
      call self%get_parameter(self%suQ7HX, 'suQ7', '1/d',        'Specific not nutrient limited refractory matter uptake')
      call self%get_parameter(self%suQ6fHX,'suQ6f','1/d',        'Specific nutrient limited detritus uptake')
      call self%get_parameter(self%suQ6sHX,'suQ6s','1/d',        'Specific not nutrient limited detritus uptake')
      call self%get_parameter(self%suQ1HX, 'suQ1', '1/d',        'Specific DOC uptake')
      call self%get_parameter(self%puincHX,'puinc','1/d',        'Preference factor of nutrient content')
      call self%get_parameter(self%p_ex6ox,'pue6', '-',          'Excreted fraction of uptake of POM')
      call self%get_parameter(self%p_ex7ox,'pue7', '-',          'Excreted fraction of uptake of refractory matter')
      call self%get_parameter(self%purHX,  'pur',  '1/d',        'Fraction of carbon uptake respired')
      call self%get_parameter(self%srHX,   'sr',   '1/d',        'Specific rest respiration')
      call self%get_parameter(self%pdHQ1X, 'pdQ1', '-',          'DOM-fraction of mortality')
      call self%get_parameter(self%sdHX,   'sd',   '1/d',        'Specific maximum mortality')

      call self%add_constituent('c',0.0_rk,qn=self%qnHIcX,qp=self%qpHIcX)

      ! Environmental dependencies
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      call self%register_dependency(self%id_Dm,'Dm','m','depth interval available to bacteria')

      ! Dependencies on state variables of external modules.
      call self%register_state_dependency(self%id_K4n,'K4n','mmol N/m^2','ammonium')
      call self%register_state_dependency(self%id_K1p,'K1p','mmol N/m^2','phosphate')
      call self%register_state_dependency(self%id_G2o,'G2o','mmol O2/m^2','oxygen')
      call self%register_state_dependency(self%id_G3c,'G3c','mmol C/m^2','dissolved inorganic carbon')
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
   end subroutine

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)

      class (type_ersem_benthic_bacteria),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
     
      real(rk) :: Hc,HcP,Hn,Hp,fHc,fHn,fHp,fHn1,fHp1
      real(rk) :: Dm,SK4a,SK1a,SQ6c
      real(rk) :: ETW,eT,eOx,eN,Limit
      real(rk) :: K4a,K1a
      real(rk) :: Q1cP,Q1nP,Q1pP
      real(rk) :: AQ6c,AQ6n,AQ6p
      real(rk) :: AQ7c,AQ7n,AQ7p
      real(rk) :: sfQ7H,sfQ6H,sfQ1H,fQ7Hc,fQ6Hc,fQ1Hc,fQIHc,fQ7Hn,fQ7Hp,fQ6Hn,fQ6Hp
      real(rk) :: fQ1Hn,fQ1Hp,fK1Hn,fK1Hp,fHG3c,fK4Hn,fK4Hp,sfHQ1,sfHQI,sfHQ6

      _HORIZONTAL_LOOP_BEGIN_

         ! Retrieve bacterial carbon
         _GET_HORIZONTAL_(self%id_c,Hc)
         _GET_HORIZONTAL_WITH_BACKGROUND_(self%id_c,HcP)

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

         ! Infer bacterial nitrogen and phosphorus from carbon and prescribed fixed stoichiometry.
         Hn = Hc*self%qnHIcX
         Hp = Hc*self%qpHIcX

         ! Limitation by temperature, oxygen, nutrient availability.
         eT  = self%q10HX**((ETW-10._rk)/10._rk) - self%q10HX**((ETW-32._rk)/3._rk)
         eOx = Dm/(self%ddHX+Dm)
         if (AQ6c<=0.0_rk) then
            ! Avoid division by zero: no carbon in substrate means no nutrient limitation
            eN = 1.0_rk
         else
            eN = min(1._rk,max(0._rk,AQ6n/(self%qnHIcX*AQ6c))) * min(1._rk,max(0._rk,AQ6p/(self%qpHIcX*AQ6c)))
         end if

         ! Combined limitation factor
         Limit = eT * eOX * Hc

         ! Determine effective uptake rates
         sfQ7H = ( self%suQ7HX * Limit )
         sfQ6H = ( self%suQ6fHX * Limit * eN ) + ( self%suQ6sHX * Limit )
         sfQ1H = ( self%suQ1HX * Limit )

         ! Uptake of substrate carbon
         fQ7Hc = sfQ7H * AQ7c
         fQ6Hc = sfQ6H * AQ6c
         fQ1Hc = sfQ1H * Q1cP
         fQIHc = fQ7Hc + fQ6Hc + fQ1Hc

         ! Uptake of substrate nitrogen and phosphorus
         ! (note on Q6: affinity for n and p differs from affinity for c if puincHX/=1)
         fQ7Hn = sfQ7H * AQ7n
         fQ7Hp = sfQ7H * AQ7p
         fQ6Hn = sfQ6H * AQ6n * self%puincHX
         fQ6Hp = sfQ6H * AQ6p * self%puincHX
         fQ1Hn = sfQ1H * Q1nP
         fQ1Hp = sfQ1H * Q1pP

         ! (JB: code below seems to want to infer nutrient requirement/flux,
         ! but it makes no sense since K?a and fK?Hn are summed while they have different units)
         fK4Hn = fQIHc * self%qnHIcX
         fK4Hn = fK4Hn * K4a/(K4a+fK4Hn)
         fK1Hp = fQIHc * self%qpHIcX
         fK1Hp = fK1Hp * K1a/(K1a+fK1Hp)

         ! Respiration (reduction in bacterial carbon and dissolved oxygen, increase in 
         ! dissolved inorganic carbon, ammonium, phosphate)
         fHG3c = self%purHX * fQIHc + self%srHX * Hc * eT
         fHc = -fHG3c
         _SET_BOTTOM_ODE_(self%id_G2o,-fHG3c/CMass)  ! oxygen or reduction equivalent
         _SET_BOTTOM_ODE_(self%id_G3c, fHG3c/CMass)

         ! Mortality (distribute over dissolved organics Q1 and particulate organics Q6)
         sfHQI = self%sdHX * (1._rk - eOx)
         sfHQ6 = sfHQI * (1._rk - self%pdHQ1X)
         sfHQ1 = sfHQI * self%pdHQ1X

         _SET_BOTTOM_ODE_(self%id_Q6c, sfHQ6 * HcP)
         _SET_BOTTOM_ODE_(self%id_Q6n, sfHQ6 * HcP * self%qnHIcX)
         _SET_BOTTOM_ODE_(self%id_Q6p, sfHQ6 * HcP * self%qpHIcX)

         _SET_BOTTOM_ODE_(self%id_Q1c, sfHQ1 * HcP)
         _SET_BOTTOM_ODE_(self%id_Q1n, sfHQ1 * HcP * self%qnHIcX)
         _SET_BOTTOM_ODE_(self%id_Q1p, sfHQ1 * HcP * self%qpHIcX)

         _SET_BOTTOM_ODE_(self%id_c, -(sfHQ6 + sfHQ1) * HcP)

         ! Uptake and excretion
         fHc = fHc + fQ7Hc * (1._rk-self%p_ex7ox) + fQ6Hc * (1._rk-self%p_ex6ox) + fQ1Hc
         fHn = fQ7Hn * (1._rk-self%p_ex7ox) + fQ6Hn * (1._rk-self%p_ex6ox) + fQ1Hn + fK4Hn
         fHp = fQ7Hp * (1._rk-self%p_ex7ox) + fQ6Hp * (1._rk-self%p_ex6ox) + fQ1Hp + fK1Hp

         ! Preserve fixed bacterial stoichiometry by selectively excreting part of N,P,C biomass fluxes
         ! to ammonia, phosphate and particulate carbon, respectively.
         SK4a = 0.0_rk
         SK1a = 0.0_rk
         SQ6c = 0.0_rk
         CALL adjust_fixed_nutrients(fHc,fHn,fHp,self%qnHIcX,self%qpHIcX,SK4a,SK1a,SQ6c)

         _SET_BOTTOM_ODE_(self%id_c,fHc)

         _SET_BOTTOM_ODE_(self%id_Q6c,-fQ6Hc + SQ6c)
         _SET_BOTTOM_ODE_(self%id_Q6n,-fQ6Hn)
         _SET_BOTTOM_ODE_(self%id_Q6p,-fQ6Hp)

         _SET_BOTTOM_ODE_(self%id_Q7c,-fQ7Hc)
         _SET_BOTTOM_ODE_(self%id_Q7n,-fQ7Hn)
         _SET_BOTTOM_ODE_(self%id_Q7p,-fQ7Hp)

         _SET_BOTTOM_ODE_(self%id_Q1c,-fQ1Hc + fQ7Hc*self%p_ex7ox + fQ6Hc*self%p_ex6ox)
         _SET_BOTTOM_ODE_(self%id_Q1n,-fQ1Hn + fQ7Hn*self%p_ex7ox + fQ6Hn*self%p_ex6ox)
         _SET_BOTTOM_ODE_(self%id_Q1p,-fQ1Hp + fQ7Hp*self%p_ex7ox + fQ6Hp*self%p_ex6ox)

         _SET_BOTTOM_ODE_(self%id_K4n,-fK4Hn + SK4a)
         _SET_BOTTOM_ODE_(self%id_K1p,-fK1Hp + SK1a)

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Adjust_fixed_nutrients \label{sec:AdjustFixedNutrients}
!
! !DESCRIPTION:
!  TODO - description
!
!  This routine determines the amount of excess c,n or p in the
!  source terms and excretes the appropriate amount(s) to Q6 and
!  nutrients, so that the fixed nutrient ratio is re-established
!
!  IN:  Source terms of fixed quota bio-state.(cnp)......SXx
!       Fixed quota values (np)..........................qx
!
!  INT: Excess nutrient in bio-state (np)................ExcessX
!
!  OUT: Source terms for Predator (SX), Q6 and nutrients
!\\
!\\
! !INTERFACE:
      subroutine adjust_fixed_nutrients ( SXc, SXn, SXp,qn, qp, SKn, SKp, SQc )
!
! !INPUT/OUTPUT PARAMETERS:
       real(rk), intent(inout)  :: SXc, SXn, SXp, SKn, SKp, SQc
       real(rk), intent(in)     :: qn, qp
!
! !LOCAL PARAMETERS:
       real(rk) :: ExcessN, ExcessP, ExcessC
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
       ExcessC = max(max(SXc - SXp/qp,SXc - SXn/qn),0._rk)
!
       if (ExcessC>0.0_rk) then
         SXc = SXc - ExcessC
         SQc = SQc + ExcessC
       end if

       ExcessN = max(SXn - SXc*qn,0._rk)
       ExcessP = max(SXp - SXc*qp,0._rk)

       if (ExcessN>0.0_rk) then
         SXn = SXn - ExcessN
         SKn = SKn + ExcessN
       end if

       if (ExcessP>0.0_rk) then
         SXp = SXp - ExcessP
         SKp = SKp + ExcessP
       end if

       end subroutine adjust_fixed_nutrients
!
!EOC
!-----------------------------------------------------------------------
   
end module
