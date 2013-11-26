#include "fabm_driver.h"

#define IRON

module pml_ersem_mesozoo

! Jorn TODO: overwintering, deal with implicit N and P in some prey [esp. when cannibalizing]

   use fabm_types
   use pml_ersem_base

   implicit none

   private

   type,extends(type_ersem_base_model),public :: type_pml_ersem_mesozoo
      ! Variables
      type (type_model_id),         allocatable,dimension(:) :: id_prey
      type (type_dependency_id),    allocatable,dimension(:) :: id_preyc,id_preyn,id_preyp,id_preys,id_preyf
      type (type_state_variable_id),allocatable,dimension(:) :: id_preyf_target
      type (type_state_variable_id)      :: id_O3c,id_O2o
      type (type_state_variable_id)      :: id_R1c,id_R1p,id_R1n
      type (type_state_variable_id)      :: id_R2c
      type (type_state_variable_id)      :: id_R8c,id_R8p,id_R8n
      type (type_state_variable_id)      :: id_R6s
      type (type_state_variable_id)      :: id_N1p,id_N4n
      type (type_dependency_id)          :: id_ETW,id_eO2mO2

      ! Parameters
      integer  :: nprey
      real(rk) :: qpZIcX,qnZIcX
      real(rk) :: q10Z4X,chrZ4oX,minfoodZ4X,chuZ4cX
      real(rk),allocatable :: suprey_Z4X(:)
      logical,allocatable :: preyispom(:)
      real(rk) :: sumZ4X
      real(rk) :: sdZ4oX, sdZ4X, srsZ4X
      real(rk) :: puZ4X
      real(rk) :: pu_eaZ4X,pu_eaRZ4X
      real(rk) :: pe_R1Z4X

      ! ERSEM global parameters
      real(rk) :: R1R2X,xR1pX,xR1nX,urB1_O2X
   contains
!     Model procedures
      procedure :: initialize
      procedure :: do
   end type

   real(rk),parameter :: CMass = 12._rk
   real(rk),parameter :: ZeroX = 1e-8_rk

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_pml_ersem_mesozoo),intent(inout),target :: self
      integer,                       intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
      integer           :: iprey
      character(len=16) :: index
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%nprey,     'nprey',default=0)
      call self%get_parameter(self%qpZIcX,    'qpZIcX')
      call self%get_parameter(self%qnZIcX,    'qnZIcX')
      call self%get_parameter(self%q10z4X,    'q10z4X')
      call self%get_parameter(self%chrZ4oX,   'chrZ4oX')
      call self%get_parameter(self%minfoodZ4X,'minfoodZ4X')
      call self%get_parameter(self%chuZ4cX,   'chuZ4cX')
      call self%get_parameter(self%sdZ4oX,    'sdZ4oX')
      call self%get_parameter(self%sdZ4X,     'sdZ4X')
      call self%get_parameter(self%srsZ4X,    'srsZ4X')
      call self%get_parameter(self%puZ4X,     'puZ4X')
      call self%get_parameter(self%pe_R1Z4X,  'pe_R1Z4X')
      call self%get_parameter(self%sumZ4X,    'sumZ4X')
      call self%get_parameter(self%pu_eaZ4X,  'pu_eaZ4X')
      call self%get_parameter(self%pu_eaRZ4X, 'pu_eaRZ4X')

      ! Register state variables
      call self%initialize_ersem_base(c_ini=1.e-4_rk,qn=self%qnZIcX,qp=self%qpZIcX)

      ! Register links to carbon contents of prey.
      allocate(self%id_prey(self%nprey))
      allocate(self%id_preyc(self%nprey))
      allocate(self%id_preyn(self%nprey))
      allocate(self%id_preyp(self%nprey))
      allocate(self%id_preyf(self%nprey))
      allocate(self%id_preys(self%nprey))
      allocate(self%id_preyf_target(self%nprey))
      allocate(self%suprey_Z4X(self%nprey))
      allocate(self%preyispom(self%nprey))
      do iprey=1,self%nprey
         write (index,'(i0)') iprey
         call self%register_model_dependency(self%id_prey(iprey),'prey'//trim(index))
         call self%register_dependency(self%id_preyc(iprey),  'prey'//trim(index)//'c',  'mg C/m^3',   'Prey '//trim(index)//' C')    
         call self%register_dependency(self%id_preyn(iprey),  'prey'//trim(index)//'n',  'mmol N/m^3', 'Prey '//trim(index)//' N')    
         call self%register_dependency(self%id_preyp(iprey),  'prey'//trim(index)//'p',  'mmol P/m^3', 'Prey '//trim(index)//' P')    
         call self%register_dependency(self%id_preyf(iprey),  'prey'//trim(index)//'f',  'mmol Fe/m^3','Prey '//trim(index)//' Fe')    
         call self%register_dependency(self%id_preys(iprey),  'prey'//trim(index)//'s',  'mmol Si/m^3','Prey '//trim(index)//' Si')
         call self%request_coupling('prey'//trim(index)//'c','c',source=self%id_prey(iprey))
         call self%request_coupling('prey'//trim(index)//'n','n',source=self%id_prey(iprey))
         call self%request_coupling('prey'//trim(index)//'p','p',source=self%id_prey(iprey))
         call self%request_coupling('prey'//trim(index)//'s','s',source=self%id_prey(iprey))
         call self%request_coupling('prey'//trim(index)//'f','f',source=self%id_prey(iprey))
         call self%register_state_dependency(self%id_preyf_target(iprey),'prey'//trim(index)//'f_sink','mmol Fe/m^3','Target pool for Fe of assimilated prey '//trim(index))    
         call self%get_parameter(self%suprey_Z4X(iprey),'suprey'//trim(index)//'_Z4X')
         call self%get_parameter(self%preyispom(iprey),'prey'//trim(index)//'ispom',default=.false.)
      end do

      ! Register links to external nutrient pools.
      call self%register_state_dependency(self%id_N1p,'N1p','mmol P/m^3', 'Phosphate')    
      call self%register_state_dependency(self%id_N4n,'N4n','mmol N/m^3', 'Ammonium')    

      ! Register links to external labile dissolved organic matter pools.
      call self%register_state_dependency(self%id_R1c,'R1c','mg C/m^3',  'DOC')
      call self%register_state_dependency(self%id_R1p,'R1p','mmol P/m^3','DOP')    
      call self%register_state_dependency(self%id_R1n,'R1n','mmol N/m^3','DON')    

      ! Register links to external semi-labile dissolved organic matter pools.
      call self%register_state_dependency(self%id_R2c,'R2c','mg C/m^3','Semi-labile DOC')    

      ! Register links to external particulate organic matter pools.
      call self%register_state_dependency(self%id_R8c,'RPc','mg C/m^3',   'POC')    
      call self%register_state_dependency(self%id_R8p,'RPp','mmol P/m^3', 'POP')    
      call self%register_state_dependency(self%id_R8n,'RPn','mmol N/m^3', 'PON')    
      call self%register_state_dependency(self%id_R6s,'RPs','mmol Si/m^3','POSi')    

      ! Register links to external total dissolved inorganic carbon, dissolved oxygen pools
      call self%register_state_dependency(self%id_O3c,'O3c','mmol C/m^3','Carbon Dioxide')    
      call self%register_state_dependency(self%id_O2o,'O2o','mmol O/m^3','Oxygen')    

      ! Register environmental dependencies (temperature, shortwave radiation)
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      call self%register_dependency(self%id_eO2mO2,standard_variables%fractional_saturation_of_oxygen)

   end subroutine
   
   subroutine do(self,_ARGUMENTS_DO_)

      class (type_pml_ersem_mesozoo),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

   ! !LOCAL VARIABLES:
      integer  :: iprey
      real(rk) :: ETW,eO2mO2
      real(rk) :: Z4c,Z4cP
      real(rk),dimension(self%nprey) :: preycP,preypP,preynP,preysP,preyfP
      real(rk) :: SZ4c,SZ4n,SZ4p
      real(rk) :: etZ4,CORROX,eO2Z4
      real(rk),dimension(self%nprey) :: spreyZ4,rupreyZ4c,fpreyZ4c
      real(rk) :: rumZ4,put_uZ4,rugZ4
      real(rk) :: sdZ4,rdZ4
      real(rk) :: ineffZ4
      real(rk) :: rrsZ4,rraZ4
      real(rk) :: fZ4O3c
      real(rk) :: retZ4,fZ4RDc,fZ4R8c
      real(rk) :: fZ4RIp,fZ4RDp,fZ4R8p
      real(rk) :: fZ4RIn,fZ4RDn,fZ4R8n
      real(rk) :: fZ4NIn,fZ4N1p,SR8c

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Get environment (temperature, oxygen saturation)
         _GET_(self%id_ETW,ETW)
         _GET_(self%id_eO2mO2,eO2mO2)

         ! Get own concentrations
         _GET_(self%id_c,Z4c)
         _GET_SAFE_(self%id_c,Z4cP)

         ! Get prey concentrations
         do iprey=1,self%nprey
            _GET_SAFE_(self%id_preyc(iprey),  preycP(iprey))
            _GET_SAFE_(self%id_preyp(iprey),  preypP(iprey))
            _GET_SAFE_(self%id_preyn(iprey),  preynP(iprey))
            _GET_SAFE_(self%id_preys(iprey),  preysP(iprey))
            _GET_SAFE_(self%id_preyf(iprey),  preyfP(iprey))
         end do

!..Temperature effect :

         etZ4 = self%q10Z4X**((ETW-10._rk)/10._rk) - self%q10Z4X**((ETW-32._rk)/3._rk)

!..Oxygen limitation :
         CORROX = 1._rk + self%chrZ4oX
         eO2Z4 = MIN(1._rk,CORROX*(eO2mO2/(self%chrZ4oX + eO2mO2)))

!..Available food :
         spreyZ4 = self%suprey_Z4X*preycP/(preycP+self%minfoodZ4X)
         rupreyZ4c = spreyZ4*preycP
         rumZ4 = sum(rupreyZ4c)

!..Uptake :
         put_uZ4 = self%sumZ4X/(rumZ4 + self%chuZ4cX)*etZ4*Z4c
         rugZ4 = put_uZ4*rumZ4

!..Fluxes into mesoplankton :
         fpreyZ4c = put_uZ4*rupreyZ4c
         spreyZ4 = put_uZ4*spreyZ4

!..Zooplankton Grazing
!      Z4herb(i) = fP1Z4c(i) + fP2Z4c(i) + fP3Z4c(i) +fP4Z4c(i)
!      Z4carn(i) = fB1Z4c(i) + fZ5Z4c(i) + fZ6Z4c(i) + fR6Z4c(i)

!..Mortality
         sdZ4 = ((1._rk - eO2Z4)*self%sdZ4oX + self%sdZ4X) 
         rdZ4 = sdZ4*Z4cP

!..Assimilation inefficiency:
         ineffZ4 = (1._rk - self%puZ4X)

!..Excretion
         retZ4 = 0.0_rk
         do iprey=1,self%nprey
            if (self%preyispom(iprey)) then
               retZ4 = retZ4 + ineffZ4 * fpreyZ4c(iprey) * self%pu_eaZ4X
            else
               retZ4 = retZ4 + ineffZ4 * fpreyZ4c(iprey) * self%pu_eaRZ4X
            end if
         end do
         fZ4RDc = (retZ4 + rdZ4)*self%pe_R1Z4X
         fZ4R8c = (retZ4 + rdZ4)*(1._rk - self%pe_R1Z4X)
#ifdef SAVEFLX
         fZXRDc(I) = fZXRDc(I)+fZ4RDc
#endif

!..Rest respiration, corrected for prevailing temperature
         rrsZ4 = self%srsZ4X*etZ4*Z4cP

!..Activity respiration
         rraZ4 = ineffZ4 * rugZ4 - retZ4

!..Total respiration
         fZ4O3c = rrsZ4 + rraZ4

#ifdef CALC
         fO3L2c(I) = fO3L2c(I) + &
            RainR(I) * gutdiss * (1._rk - puZ4X)* pu_eaZ4X * fP2Z4c
#endif
         rraZ4 = rugZ4*(1._rk - self%puZ4X)-retZ4

!..Source equations
         SZ4c = rugZ4 - fZ4RDc - fZ4R8c - fZ4O3c   ! Jorn: cannibalism accounted for in grazing/predation section

!..Flows from and to detritus
         _SET_ODE_(self%id_R1c, + fZ4RDc * self%R1R2X)
         _SET_ODE_(self%id_R2c, + fZ4RDc * (1._rk-self%R1R2X))
         _SET_ODE_(self%id_R8c, + fZ4R8c)

!..Grazing and predation
         do iprey=1,self%nprey
            !_SET_ODE_(self%id_preyc(iprey),   - fpreyZ4c(iprey))
            !_SET_ODE_(self%id_preyChl(iprey), - spreyZ4(iprey)*preyChlP(iprey))
         end do

!..Respiration
         _SET_ODE_(self%id_O3c, + fZ4O3c/CMass)
         _SET_ODE_(self%id_O2o, - fZ4O3c*self%urB1_O2X)

!..Phosphorus dynamics in mesozooplankton, derived from carbon flows....

         fZ4RIp = (fZ4RDc + fZ4R8c) * self%qpZIcX
         fZ4RDp = min(fZ4RIp, fZ4RDc * self%qpZIcX * self%xR1pX)
         fZ4R8p = fZ4RIp - fZ4RDp

!..Source equations
         SZ4p = sum(spreyZ4*preypP) - fZ4R8p - fZ4RDp

#ifdef IRON
!  Iron dynamics
! following Vichi et al., 2007 it is assumed that the iron fraction of the ingested phytoplankton
! is egested as particulate detritus (Luca)
         do iprey=1,self%nprey
            !_SET_ODE_(self%id_preyf(iprey),       - spreyZ4(iprey)*preyfP(iprey))
            _SET_ODE_(self%id_preyf_target(iprey),  spreyZ4(iprey)*preyfP(iprey))
         end do
#endif

!..Phosphorus flux from/to detritus
         _SET_ODE_(self%id_R1p, + fZ4RDp)
         _SET_ODE_(self%id_R8p, + fZ4R8p)

!..Phosphorus flux from prey
         do iprey=1,self%nprey
            !_SET_ODE_(self%id_preyp(iprey),- spreyZ4(iprey)*preypP(iprey))
         end do

!..Nitrogen dynamics in mesozooplankton, derived from carbon flows......

         fZ4RIn = (fZ4RDc + fZ4R8c) * self%qnZIcX
         fZ4RDn = min(fZ4RIn, fZ4RDc * self%qnZIcX * self%xR1nX)
         fZ4R8n = fZ4RIn - fZ4RDn

!..Source equations
         SZ4n = sum(spreyZ4*preynP) - fZ4R8n - fZ4RDn

!..Nitrogen flux from/to detritus
         _SET_ODE_(self%id_R1n, + fZ4RDn)
         _SET_ODE_(self%id_R8n, + fZ4R8n)

!..Nitrogen flux from prey
         do iprey=1,self%nprey
            !_SET_ODE_(self%id_preyn(iprey),- spreyZ4(iprey)*preynP(iprey))
         end do

!..Silica-flux from diatoms due to mesozooplankton grazing
         do iprey=1,self%nprey
            !_SET_ODE_(self%id_preys(iprey),- spreyZ4(iprey)*preysP(iprey))
         end do
         _SET_ODE_(self%id_R6s,sum(spreyZ4*preysP))

!..re-establish the fixed nutrient ratio in zooplankton.................

         fZ4NIn = 0.0_rk
         fZ4N1p = 0.0_rk
         SR8c = 0.0_rk
         CALL Adjust_fixed_nutrients ( SZ4c, SZ4n, SZ4p, self%qnZIcX, &
                                    self%qpZIcX, fZ4NIn, fZ4N1p, SR8c)

         _SET_ODE_(self%id_c,SZ4c)
         _SET_ODE_(self%id_N4n,fZ4NIn)
         _SET_ODE_(self%id_N1p,fZ4N1p)
         _SET_ODE_(self%id_R8c,SR8c)

      ! Leave spatial loops (if any)
      _LOOP_END_

   end subroutine do

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
       IF ( ExcessC .GT. ZeroX ) THEN
         SXc = SXc - ExcessC
         SQc = SQc + ExcessC
       END If

       ExcessN = max(SXn - SXc*qn,0._rk)
       ExcessP = max(SXp - SXc*qp,0._rk)

       IF ( ExcessN .GT. ZeroX ) THEN
         SXn = SXn - ExcessN
         SKn = SKn + ExcessN
       END If

       IF ( ExcessP .GT. ZeroX ) THEN
         SXp = SXp - EXcessP
         SKp = SKp + ExcessP
       END If

       RETURN

       END SUBROUTINE Adjust_fixed_nutrients
!
!EOC
!-----------------------------------------------------------------------
   
end module