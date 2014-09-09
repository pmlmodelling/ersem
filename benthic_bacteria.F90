#include "fabm_driver.h"

module pml_ersem_benthic_bacteria

   use fabm_types

   use pml_ersem_shared

   implicit none

   private

   type,extends(type_base_model),public :: type_pml_ersem_benthic_bacteria
      type (type_bottom_state_variable_id) :: id_Hc
      type (type_bottom_state_variable_id) :: id_Dm,id_K4a,id_K1a
      type (type_dependency_id) :: id_ETW
      type (type_bottom_state_variable_id) :: id_AQ6c,id_AQ6n,id_AQ6p,id_AQ7c,id_AQ7n,id_AQ7p,id_Q1c,id_Q1n,id_Q1p
      type (type_bottom_state_variable_id) :: id_G2o,id_G3c

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
      procedure :: do_benthic_bacteria
   end type

contains

   subroutine initialize(self,configunit)
      class (type_pml_ersem_benthic_bacteria),intent(inout),target :: self
      integer,                                 intent(in)           :: configunit
      call self%get_parameter(self%qnHIcX,'qnHIcX','mmol N/mg C','Maximal nitrogen to carbon ratio of benthic bacteria')
      call self%get_parameter(self%qpHIcX,'qpHIcX','mmol P/mg C','Maximal phosphorus to carbon ratio of benthic bacteria')
      call self%get_parameter(self%q10HX,'q10HX','-','Regulating temperature factor Q10 for benthic bacteria')
      call self%get_parameter(self%ddHX,'ddHX','1/m','Michaelis-Menten constant for limitation through layer thickness')
      call self%get_parameter(self%suQ7HX,'suQ7HX','1/d','Specific not nutrient limited refractory matter uptake by benthic bacteria')
      call self%get_parameter(self%suQ6fHX,'suQ6fHX','1/d','Specific nutrient limited detritus uptake by benthic bacteria')
      call self%get_parameter(self%suQ6sHX,'suQ6sHX','1/d','Specific not nutrient limited detritus uptake by benthic bacteria')
      call self%get_parameter(self%suQ1HX,'suQ1HX','1/d','Specific DOC uptake by benthic aerobic bacteria')
      call self%get_parameter(self%puincHX,'puincHX','1/d','Preference factor of nutrient content by benthic bacteria')
      call self%get_parameter(self%p_ex7ox,'p_ex7ox','-','Excreted fraction of uptake of refractory matter by benthic bacteria')
      call self%get_parameter(self%p_ex6ox,'p_ex6ox','-','Excreted fraction of uptake of POM by benthic bacteria')
      call self%get_parameter(self%purHX,'purHX','1/d','Fraction of carbon uptake respired by benthic bacetria')
      call self%get_parameter(self%srHX,'srHX','1/d','Specific rest respiration of benthic bacetria')
      call self%get_parameter(self%pdHQ1X,'pdHQ1X','-','DOM-fraction of benthic bacteria mortality')
      call self%get_parameter(self%sdHX,'sdHX','1/d','Specific maximal mortality of benthic bacteria')
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      call self%register_state_dependency(self%id_Dm,'Dm','m','depth interval available to bacteria')
      call self%register_state_dependency(self%id_K4a,'K4a','m','depth interval available to bacteria')
      call self%register_state_dependency(self%id_K4a,'K4a','m','depth interval available to bacteria')
   end subroutine

   subroutine do_benthic_bacteria(self,_ARGUMENTS_DO_BOTTOM_)

! !INPUT PARAMETERS:
 class (type_pml_ersem_benthic_bacteria),intent(in) :: self
  _DECLARE_ARGUMENTS_DO_BOTTOM_
     
     real(rk) :: Hc,Hn,Hp,fHc,fHn,fHp,fHn1,fHp1
     real(rk) :: Q1cP,Q1nP,Q1pP,Dm,SK4a,SK1a,SQ6c
     real(rk) :: ETW,eT,eOx,eN,Limit
     real(rk) :: K4a,K1a,AQ6c,AQ6n,AQ6p,AQ7c,AQ7n,AQ7p,G3c,G2o
     real(rk) :: sfQ7H,sfQ6H,fQ7Hc,fQ6Hc,fQ1Hc,fQIHc,fQ7Hn,fQ7Hp,fQ6Hn,fQ6Hp
     real(rk) :: fQ1Hn,fQ1Hp,fK1Hn,fK1Hp,fHG3c,fK4Hn,fK4Hp,sfHQ1,sfHQI,sfHQ6

     _HORIZONTAL_LOOP_BEGIN_

     _GET_HORIZONTAL_(self%id_Hc,Hc)

     _GET_HORIZONTAL_(self%id_Q1c,Q1cP)
     _GET_HORIZONTAL_(self%id_Q1n,Q1nP)
     _GET_HORIZONTAL_(self%id_Q1p,Q1pP)

     _GET_HORIZONTAL_(self%id_Dm,Dm)

     _GET_(self%id_ETW,ETW)

     _GET_HORIZONTAL_(self%id_AQ6c,AQ6c)
     _GET_HORIZONTAL_(self%id_AQ6n,AQ6n)
     _GET_HORIZONTAL_(self%id_AQ6p,AQ6p)
     _GET_HORIZONTAL_(self%id_AQ7c,AQ7c)
     _GET_HORIZONTAL_(self%id_AQ7n,AQ7n)
     _GET_HORIZONTAL_(self%id_AQ7p,AQ7p)
     _GET_HORIZONTAL_(self%id_K4a,K4a)
     _GET_HORIZONTAL_(self%id_K1a,K1a)

     _GET_HORIZONTAL_(self%id_G2o,G2o)
     _GET_HORIZONTAL_(self%id_G3c,G3c)

      Hn = Hc * self%qnHIcX
      Hp = Hc * self%qpHIcX

      eT   = self%q10HX**((ETW-10._rk)/10._rk) - self%q10HX**((ETW-32._rk)/3._rk)
      eOx  = Dm/(self%ddHX+Dm)
      eN   = min(1._rk,max(0._rk,AQ6n/(self%qnHIcX*AQ6c))) * min(1._rk,max(0._rk,AQ6p/(self%qpHIcX*AQ6c)))

      Limit = eT * eOX * Hc
   
   sfQ7H = ( self%suQ7HX * Limit )
   sfQ6H = ( self%suQ6fHX * Limit * eN ) + ( self%suQ6sHX * Limit )
   fQ7Hc = sfQ7H * AQ7c
   fQ6Hc = sfQ6H * AQ6c
   fQ1Hc = self%suQ1HX  * Limit * Q1cP
   fQIHc = fQ7Hc + fQ6Hc + fQ1Hc

   fQ7Hn = sfQ7H * AQ7n
   fQ7Hp = sfQ7H * AQ7p
   fQ6Hn = sfQ6H * AQ6n * self%puincHX
   fQ6Hp = sfQ6H * AQ6p * self%puincHX

   fQ1Hn = self%suQ1HX * Limit * Q1nP
   fQ1Hp = self%suQ1HX * Limit * Q1pP

   fK4Hn = fQIHc * self%qnHIcX
   fK4Hn = fK4Hn * K4a/(K4a+fK4Hn)
   fK1Hp = fQIHc * self%qpHIcX
   fK1Hp = fK1Hp * K1a/(K1a+fK1Hp)

   !Respiration

   fHG3c = self%purHX * fQIHc + self%srHX * Hc * eT
   _SET_ODE_BEN_(self%id_Hc, - fHG3c)
   !Oxygen or reduction equivalent
   _SET_ODE_BEN_(self%id_G2o, - fHG3c / 12.011_rk)
   _SET_ODE_BEN_(self%id_G3c, fHG3c / 12.011_rk)
  
   !Mortality

   sfHQI = self%sdHX * (1._rk - eOx)
   sfHQ6 = sfHQI * (1._rk - self%pdHQ1X)
   sfHQ1 = sfHQI * self%pdHQ1X
   
   _SET_ODE_BEN_(self%id_AQ6c, sfHQ6 * Hc)
   _SET_ODE_BEN_(self%id_AQ6n, sfHQ6 * Hc * self%qnHIcX)
   _SET_ODE_BEN_(self%id_AQ6p, sfHQ6 * Hc * self%qpHIcX)

   _SET_ODE_BEN_(self%id_Q1c, sfHQ1 * Hc)
   _SET_ODE_BEN_(self%id_Q1n, sfHQ1 * Hc * self%qnHIcX)
   _SET_ODE_BEN_(self%id_Q1p, sfHQ1 * Hc * self%qpHIcX)

   _SET_ODE_BEN_(self%id_Hc,- (sfHQ6 + sfHQ1) * Hc )

   fHn1 = - (sfHQ6 + sfHQ1) * Hc * self%qnHIcX
   fHp1 = - (sfHQ6 + sfHQ1) * Hc * self%qpHIcX

   !Uptake and excretion

   fHc = fQ7Hc * (1._rk-self%p_ex7ox) + fQ6Hc * (1._rk-self%p_ex6ox) + fQ1Hc
   fHn = fHn1 + (fQ7Hn * (1._rk-self%p_ex7ox) + fQ6Hn * (1._rk-self%p_ex6ox) + fQ1Hn + fK4Hn)
   fHp = fHp1 + (fQ7Hp * (1._rk-self%p_ex7ox) + fQ6Hp * (1._rk-self%p_ex6ox) + fQ1Hp + fK1Hp)
  
   CALL Adjust_fixed_nutrients(fHc,fHn,fHp,self%qnHIcX,self%qpHIcX,SK4a,SK1a,SQ6c)

   _SET_ODE_BEN_(self%id_Hc,fHc)

   _SET_ODE_BEN_(self%id_AQ6c,-fQ6Hc + SQ6c)
   _SET_ODE_BEN_(self%id_AQ6n,-fQ6Hn)
   _SET_ODE_BEN_(self%id_AQ6p,-fQ6Hp)

   _SET_ODE_BEN_(self%id_AQ7c,-fQ7Hc)
   _SET_ODE_BEN_(self%id_AQ7n,-fQ7Hn)
   _SET_ODE_BEN_(self%id_AQ7p,-fQ7Hp)

   _SET_ODE_BEN_(self%id_Q1c,fQ7Hc*self%p_ex7ox + fQ6Hc*self%p_ex6ox - fQ1Hc)
   _SET_ODE_BEN_(self%id_Q1n,fQ7Hn*self%p_ex7ox + fQ6Hn*self%p_ex6ox - fQ1Hn)
   _SET_ODE_BEN_(self%id_Q1p,fQ7Hp*self%p_ex7ox + fQ6Hp*self%p_ex6ox - fQ1Hp)

   _SET_ODE_BEN_(self%id_K4a,-fK4Hn + SK4a)
   _SET_ODE_BEN_(self%id_K1a,-fK1Hp + SK1a)

     _HORIZONTAL_LOOP_END_

   end subroutine

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
      SUBROUTINE Adjust_fixed_nutrients ( SXc, SXn, SXp,qn, qp, SKn, SKp, SQc )
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
       IF ( ExcessC .GT. 0.0_rk ) THEN
         SXc = SXc - ExcessC
         SQc = SQc + ExcessC
       END If

       ExcessN = max(SXn - SXc*qn,0._rk)
       ExcessP = max(SXp - SXc*qp,0._rk)

       IF ( ExcessN .GT. 0.0_rk ) THEN
         SXn = SXn - ExcessN
         SKn = SKn + ExcessN
       END If

       IF ( ExcessP .GT. 0.0_rk ) THEN
         SXp = SXp - EXcessP
         SKp = SKp + ExcessP
       END If

       END SUBROUTINE Adjust_fixed_nutrients
!
!EOC
!-----------------------------------------------------------------------
   
end module
