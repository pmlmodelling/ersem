#include "fabm_driver.h"

module pml_ersem_benthic_bacteria

   use fabm_types

   use pml_ersem_shared

   implicit none

   private

   type,extends(type_base_model),public :: type_pml_ersem_benthic_bacteria
      type (type_bottom_state_variable_id) :: id_H1c,id_Q6c,id_Q6n,id_Q6p,id_Q6s,id_Q7c,id_Q7n,id_Q7p,id_Q1c,id_Q1n,id_Q1p
      type (type_dependency_id) :: id_ETW

      real(rk) :: qnHIcX,qpHIcX
   contains
      procedure :: initialize
      procedure :: do_aerobic_bacteria
   end type

contains

   subroutine initialize(self,configunit)
      class (type_pml_ersem_benthic_bacteria),intent(inout),target :: self
      integer,                                 intent(in)           :: configunit
      call self%get_parameter(self%qnHIcX,'qnHIcX','mmol N/mg C','Maximal nitrogen to carbon ratio of benthic bacteria')
      call self%get_parameter(self%qpHIcX,'qpHIcX','mmol P/mg C','Maximal phosphorus to carbon ratio of benthic bacteria')
      call self%get_parameter(self%q10H1X,'q10H1X','-','Regulating temperature factor Q10 for benthic aerobic bacteria')
      call self%get_parameter(self%ddH1X,'ddH1X','1/m','Michaelis-Menten constant for oxygen limitation through aerobic layer depth')
      call self%get_parameter(self%suQ7H1X,'suQ7H1X','1/d','Specific not nutrient limited refractory matter uptake by benthic aerobic bacteria')
      call self%get_parameter(self%suQ6fH1X,'suQ6fH1X','1/d','Specific nutrient limited detritus uptake by benthic aerobic bacteria')
      call self%get_parameter(self%suQ6sH1X,'suQ6sH1X','1/d','Specific not nutrient limited detritus uptake by benthic aerobic bacteria')
      call self%get_parameter(self%suQ1H1X,'suQ1H1X','1/d','Specific DOC uptake by benthic aerobic bacteria')
      call self%get_parameter(self%puincH1X,'puincH1X','1/d','Preference factor of nutrient content by benthic aerobic bacteria')
   end subroutine

   subroutine do_aerobic_bacteria(self,_ARGUMENTS_DO_BOTTOM_)

! !INPUT PARAMETERS:
 class (type_pml_ersem_benthic_bacteria),intent(in) :: self
  _DECLARE_ARGUMENTS_DO_BOTTOM_
     
     real(rk) :: H1c,Q6c,Q6n,Q6p,Q6s,Q7c,Q7n,Q7p, Q1cP, Q1nP, Q1pP
     real(rk) :: D1m, D2m
     real(rk) :: ETW, eT, eOx, eN, Limit
     real(rk) :: H1n, H1p
     real(rk) :: fH1Q6(4), fH1Q7(4)
     real(rk) :: AQ6c,AQ6n,AQ6p,AQ7c,AQ7n,AQ7p
     real(rk) :: sfQ7H1,sfQ6H1,fQ7H1c,fQ6H1c,fQ1H1c,fQIH1c,fQ7H1n,fQ7H1p,fQ6H1n,fQ6H1p,fQ1H1n,fQ1H1p

     _HORIZONTAL_LOOP_BEGIN_

     _GET_HORIZONTAL_(self%id_H1c,H1c)
     _GET_HORIZONTAL_(self%id_Q6c,Q6c)
     _GET_HORIZONTAL_(self%id_Q6n,Q6n)
     _GET_HORIZONTAL_(self%id_Q6p,Q6p)
     _GET_HORIZONTAL_(self%id_Q6s,Q6s)
     _GET_HORIZONTAL_(self%id_Q7c,Q7c)
     _GET_HORIZONTAL_(self%id_Q7n,Q7n)
     _GET_HORIZONTAL_(self%id_Q7p,Q7p)
     _GET_HORIZONTAL_(self%id_Q1c,Q1cP)
     _GET_HORIZONTAL_(self%id_Q1n,Q1nP)
     _GET_HORIZONTAL_(self%id_Q1n,Q1nP)
     _GET_HORIZONTAL_(self%id_D1m,D1m)
     _GET_HORIZONTAL_(self%id_D2m,D2m)

     _GET_(self%id_ETW,ETW)
     
      H1n = H1c * self%qnHIcX
      H1p = H1c * self%qpHIcX
   
      fH1Q6(1) = Q6c
      fH1Q6(2) = Q6n
      fH1Q6(3) = Q6p
      fH1Q6(4) = Q6s

      fH1Q7(1) = Q7c
      fH1Q7(2) = Q7n
      fH1Q7(3) = Q7p
      fH1Q7(4) = 0.0_rk

      CALL AvQ6( 0.0_rk, D1m, AQ6c, AQ6n, AQ6p )
      CALL AvQ7( 0.0_rk, D1m, AQ7c, AQ7n, AQ7p )

      eT   = self%q10H1X**((ETW(I)-10._rk)/10._rk) - self%q10H1X**((ETW(I)-32._rk)/3._rk)
      eOx  = D1m/(self%ddH1X+D1m)
      eN   = eramp(AQ6n, self%qnHIcX*AQ6c) * eramp(AQ6p, self%qpHIcX*AQ6c)

      Limit = eT * eOX * H1c
   
   sfQ7H1 = ( self%suQ7H1X * Limit )
   sfQ6H1 = ( self%suQ6fH1X * Limit * eN ) + ( self%suQ6sH1X * Limit )
   fQ7H1c = sfQ7H1 * AQ7c
   fQ6H1c = sfQ6H1 * AQ6c
   fQ1H1c = self%suQ1H1X  * Limit * Q1cP
   fQIH1c = fQ7h1c + fQ6H1c + fQ1H1c

   fQ7H1n = sfQ7H1 * AQ7n
   fQ7H1p = sfQ7H1 * AQ7p
   fQ6H1n = sfQ6H1 * AQ6n * self%puincH1X
   fQ6H1p = sfQ6H1 * AQ6p * self%puincH1X

   fQ1H1n = self%suQ1H1X * Limit * Q1nP
   fQ1H1p = self%suQ1H1X * Limit * Q1pP

   fK4H1n = fQIH1c * self%qnHIcX
   fK4H1n = fK4H1n * eMM(K4/calculate_adsorption(4,K,0._rk,D1m),fK4H1n)
   fK1H1p = fQIH1c * self%qpHIcX
   fK1H1p = fK1H1p * eMM(K1/calculate_adsorption(1,K,0._rk,D1m),fK1H1p)

     _HORIZONTAL_LOOP_END_

   end subroutine
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculate available detritus (Q6)
!
! !DESCRIPTION:
!  Calculates the detritus available in the layer between the depths
!  d\_top and d\_bot.
!
! !INTERFACE:
   subroutine AvQ6( d_top, d_bot, AQ6c, AQ6n, AQ6p )
!
! !USES:
   use benthic_variables, only: Q6cP,Q6nP,Q6pP,D6m,D7m,D8m
   use benthic_parameters, only: d_totx
!
! !INPUT PARAMETERS:
!  Compartment identifier
   integer,intent(in) :: k

!  Top limit of detritus utilisation
   real(rk),intent(in) :: d_top

!  Bottom limit of detritus utilisation
   real(rk),intent(in) :: d_bot
!
! !OUTPUT PARAMETERS:
!  Available detrital carbon
   real(rk),intent(out) :: AQ6c

!  Available detrital nitrogen
   real(rk),intent(out) :: AQ6n

!  Available detrital phosphorous
   real(rk),intent(out) :: AQ6p
!
! !REVISION HISTORY:
!  Original author(s): The ERSEM development team
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
   AQ6c = Q6cP * partQ( D6m, d_top, d_bot, d_totX )
   AQ6n = Q6nP * partQ( D7m, d_top, d_bot, d_totX )
   AQ6p = Q6pP * partQ( D8m, d_top, d_bot, d_totX )

   end subroutine AvQ6
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute available detritus (Q7)
!
! !DESCRIPTION:
!  Calculates the detritus available in the layer between the depths
!  d\_top and d\_bot.
!
! !INTERFACE:
   subroutine AvQ7( d_top, d_bot, AQ7c, AQ7n, AQ7p )
!
! !USES:
   use benthic_variables, only: Q7cP,Q7nP,Q7pP,d3m,d4m,d5m
   use benthic_parameters, only: d_totx
!
! !INPUT PARAMETERS:
!  Compartment identifier
   integer,intent(in) :: k

!  Top limit of detritus utilisation
   real(rk),intent(in) :: d_top

!  Bottom limit of detritus utilisation
   real(rk),intent(in) :: d_bot
!
! !OUTPUT PARAMETERS:
!  Available detrital carbon
   real(rk),intent(out) :: AQ7c

!  Available detrital nitrogen
   real(rk),intent(out) :: AQ7n

!  Available detrital phosphorous
   real(rk),intent(out) :: AQ7p

!
! !REVISION HISTORY:
!  Original author(s): The ERSEM development team
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
   AQ7c = Q7cP * partQ( D3m, d_top, d_bot, d_totX )
   AQ7n = Q7nP * partQ( D4m, d_top, d_bot, d_totX )
   AQ7p = Q7pP * partQ( D5m, d_top, d_bot, d_totX )

   return

   end subroutine AvQ7
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute fraction of detritus between depth levels
!
! !DESCRIPTION:
!  Computes the fraction of detritus between d\_top and d\_bot. NOTE 1: 
!  Does not treat silicate. NOTE 2: Factor 13.8. This mystic factor 
!  takes care that exp(-..) won't become too
!  small so that the simulation crashes.  -13.8 is the smallest
!  number so that exp(-13.8) is evaluated correctly. It depends
!  on the accuracy of the implemented FORTRAN. 
!
! !INTERFACE:
   real(fp8) function partQ( d_pen, d_top, d_bot, d_max )
!
! !INPUT PARAMETERS:
!  Penetration depth of detrital component
   real(rk), intent(in) :: d_pen

!  Top of detrital layer
   real(rk), intent(in) ::  d_top

!  Bottom of detrital layer
   real(rk), intent(in) ::  d_bot

!  Maximum depth of detrital layer
   real(rk), intent(in) ::  d_max
!
! !LOCAL VARIABLES:
   real(rk) :: norm

   real(rk) ::  d_top1

   real(rk) ::  d_bot1

   real(rk) ::  d_max1
!
! !REVISION HISTORY:
!  Original author(s): The ERSEM development team
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
   ! Code...................................................
   d_max1 = MIN( d_pen*13.8_rk, d_max )
   d_bot1 = MIN( d_bot, d_max1 )
   d_top1 = MIN( d_top, d_bot1 )

   if ( d_max1 .gt. 0._rk ) then
      norm = 1._rk -EXP( -d_max1/d_pen )
      partQ = ( EXP( -d_top1/d_pen ) -EXP( -d_bot1/d_pen )) / norm
   else
      if ( d_top .eq. 0._rk ) then
         partQ = 1._rk
      else
         partQ = 0._rk
      end if
   end if

   return

   end function partQ
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ratio constrained to $[0,1]$
!
! !DESCRIPTION:
!  Ratio of two variables constrained to the interval $[0,1]$.
!
! !INTERFACE:
   real(rk) function eramp(x,m)
!
! !USES:
!
! !INPUT PARAMETERS:
   real(rk), intent(in) ::  x,m
!
! !REVISION HISTORY:
!  Original author(s): The ERSEM development team
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
   if (x.lt.0._fp8) then
      eramp = 0._rk
   else
      if (x.lt.m) then
         eramp = x/m
      else
         eramp = 1._rk
      end if
   end if

   return

   end function eramp
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculate absorption
!
! !DESCRIPTION:
!  Calculate average adsorption from depth 'from' to depth 'to' for box i
!  This routine is used to calculate the adsorption according the distribution
!  of detritus. The distribution of the detritus determines the distribution
!  of the bacteria.\\[\baselineskip]
!
! !INTERFACE:
   real(rk) function calculate_adsorption(mode,K,from,to)
!
! !USES:
   use ersem_constants, only: zeroX
   use benthic_variables
   use pelagic_variables, only: n_comp, n_upperX
!
! !INPUT PARAMETERS:
   integer, intent(in) :: mode,K
   real(fp8), intent(in) :: from,to
!
! !LOCAL VARIABLES:
   real(fp8) :: r,s,a,b,alfa6c,alfa7c,pM_1,pM_2
!
! !REVISION HISTORY:
!  Original author(s): The ERSEM development team
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
   if (mode.eq.1) then
      pM_1=25._rk      +1._rk
      if (N_COMP-N_UPPERX.gt.1) pM_1=benthic_morfology(2,K)+1
      pM_2=2._rk +1._rk
   else
      pM_1=3._rk +1._rk
      pM_2=3._rk +1._rk
   end if
!
   alfa6c=1._rk/D6m(K)
   alfa7c=1._rk/D3m(K)
   a=0._rk
   b=0._rk
   calculate_adsorption=.5_rk*(pM_1+pM_2)
!
   if (from .lt. D2m(K)) then
      s=min(to,D2m(K))
      r=integral_exp(-alfa6c,s-from)+integral_exp(-alfa7c,s-from)
      a=a+r
      b=b+r*pM_1
   end if

   if (to .gt. D2m(K)) then
      s=max(from,D2m(K))
      r=exp(-alfa6c*(s-from))* integral_exp(-alfa6c,to-s)&
        +exp(-alfa7c*(s-from))* integral_exp(-alfa7c,to-s)
        a=a+r
        b=b+r*pM_2
   end if

   if (a.ge.ZeroX) calculate_adsorption=b/a

   return

   end function calculate_adsorption
!
!EOC
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Integral exp
!
! !DESCRIPTION:
!  \begin{equation}
!     f(x)=\frac{\textrm{exp}(\alpha x)-1}{\alpha}
!  \end{equation}
!  \\[\baselineskip]
!
! !INTERFACE:
   real(rk) function integral_exp(alfa,d)
!
! !INPUT PARAMETERS:
   real(rk),intent(in) :: alfa,d
!
! !REVISION HISTORY:
!  Original author(s): The ERSEM development team
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
   integral_exp= (exp(alfa*d)-1._rk)/alfa

   return

   end function integral_exp
!
!EOC
!-----------------------------------------------------------------------

end module

