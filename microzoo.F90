#include "fabm_driver.h"

#define IRON

module pml_ersem_microzoo

   use fabm_types
   use pml_ersem_base

   implicit none

   private

   type,extends(type_ersem_base_model),public :: type_pml_ersem_microzoo
      ! Variables
      type (type_model_id),         allocatable,dimension(:) :: id_prey
      type (type_dependency_id),    allocatable,dimension(:) :: id_preyc,id_preyn,id_preyp,id_preys,id_preyf
      type (type_state_variable_id),allocatable,dimension(:) :: id_preyf_target
      type (type_state_variable_id)      :: id_O3c, id_O2o
      type (type_state_variable_id)      :: id_R1c, id_R2c, id_R6c
      type (type_state_variable_id)      :: id_R1p, id_R6p
      type (type_state_variable_id)      :: id_R1n, id_R6n
      type (type_state_variable_id)      :: id_R6s
      type (type_state_variable_id)      :: id_N1p,id_N4n
      type (type_dependency_id)          :: id_ETW,id_eO2mO2

      ! Parameters
      integer  :: nprey
      real(rk) :: chuz5cX,sumz5X,puz5X,pe_r1z5X,q10z5X,srsz5X,chrz5oX,pu_eaz5X
      real(rk) :: sdz5oX,sdz5X,qnz5cX,qpz5cX,minfoodz5X,stempz5nX,stempz5pX
      real(rk),allocatable :: suprey_Z5X(:)

      ! ERSEM global parameters
      real(rk) :: R1R2X,xR1pX,xR1nX,urB1_O2X
   contains
!     Model procedures
      procedure :: initialize
      procedure :: do
   end type

   real(rk),parameter :: CMass = 12._rk

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_pml_ersem_microzoo),intent(inout),target :: self
      integer,                        intent(in)           :: configunit
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
      call self%get_parameter(self%chuz5cX,   'chuz5cX')
      call self%get_parameter(self%sumz5X,    'sumz5X')
      call self%get_parameter(self%puz5X,     'puz5X')
      call self%get_parameter(self%pe_r1z5X,  'pe_r1z5X')
      call self%get_parameter(self%q10z5X,    'q10z5X')
      call self%get_parameter(self%srsz5X,    'srsz5X')
      call self%get_parameter(self%chrz5oX,   'chrz5oX')
      call self%get_parameter(self%pu_eaz5X,  'pu_eaz5X')
      call self%get_parameter(self%sdz5oX,    'sdz5oX')
      call self%get_parameter(self%sdz5X,     'sdz5X')
      call self%get_parameter(self%qnz5cX,    'qnz5cX')
      call self%get_parameter(self%qpz5cX,    'qpz5cX')
      call self%get_parameter(self%minfoodz5X,'minfoodz5X')
      call self%get_parameter(self%stempz5nX, 'stempz5nX')
      call self%get_parameter(self%stempz5pX, 'stempz5pX')

      call self%get_parameter(self%R1R2X,   'R1R2X')
      call self%get_parameter(self%xR1pX,   'xR1pX')
      call self%get_parameter(self%xR1nX,   'xR1nX')
      call self%get_parameter(self%urB1_O2X,'urB1_O2X')

      ! Register state variables
      call self%initialize_ersem_base(c_ini=1.e-4_rk,p_ini=4.288e-8_rk,n_ini=1.26e-6_rk)

      ! Register links to carbon contents of prey.
      allocate(self%id_prey(self%nprey))
      allocate(self%id_preyc(self%nprey))
      allocate(self%id_preyn(self%nprey))
      allocate(self%id_preyp(self%nprey))
      allocate(self%id_preyf(self%nprey))
      allocate(self%id_preys(self%nprey))
      allocate(self%id_preyf_target(self%nprey))
      allocate(self%suprey_Z5X(self%nprey))
      do iprey=1,self%nprey
         write (index,'(i0)') iprey
         call self%get_parameter(self%suprey_Z5X(iprey),'suprey'//trim(index)//'_Z5X')
         call self%register_dependency(self%id_preyc(iprey),  'prey'//trim(index)//'c',  'mg C m-3',   'Prey '//trim(index)//' C')
         call self%register_dependency(self%id_preyn(iprey),  'prey'//trim(index)//'n',  'mmol N m-3', 'Prey '//trim(index)//' N')
         call self%register_dependency(self%id_preyp(iprey),  'prey'//trim(index)//'p',  'mmol P m-3', 'Prey '//trim(index)//' P')
         call self%register_dependency(self%id_preyf(iprey),  'prey'//trim(index)//'f',  'mmol Fe m-3','Prey '//trim(index)//' Fe')
         call self%register_dependency(self%id_preys(iprey),  'prey'//trim(index)//'s',  'mmol Si m-3','Prey '//trim(index)//' Si')
         call self%register_state_dependency(self%id_preyf_target(iprey),'prey'//trim(index)//'f_sink','umol Fe m-3','Target pool for Fe of assimilated prey '//trim(index))    

         call self%register_model_dependency(self%id_prey(iprey),'prey'//trim(index))
         call self%request_coupling(self%id_preyc(iprey),'c',source=self%id_prey(iprey))
         call self%request_coupling(self%id_preyn(iprey),'n',source=self%id_prey(iprey))
         call self%request_coupling(self%id_preyp(iprey),'p',source=self%id_prey(iprey))
         call self%request_coupling(self%id_preys(iprey),'s',source=self%id_prey(iprey))
         call self%request_coupling(self%id_preyf(iprey),'f',source=self%id_prey(iprey))
      end do

      ! Register links to external nutrient pools.
      call self%register_state_dependency(self%id_N1p,'N1p','mmol P m-3', 'Phosphate')    
      call self%register_state_dependency(self%id_N4n,'N4n','mmol N m-3', 'Ammonium')    

      ! Register links to external labile dissolved organic matter pools.
      call self%register_state_dependency(self%id_R1c,'R1c','mg C m-3',  'DOC')
      call self%register_state_dependency(self%id_R1p,'R1p','mmol P m-3','DOP')    
      call self%register_state_dependency(self%id_R1n,'R1n','mmol N m-3','DON')    

      ! Register links to external semi-labile dissolved organic matter pools.
      call self%register_state_dependency(self%id_R2c,'R2c','mg C m-3','Semi-labile DOC')    

      ! Register links to external particulate organic matter pools.
      call self%register_state_dependency(self%id_R6c,'RPc','mg C m-3',   'POC')    
      call self%register_state_dependency(self%id_R6p,'RPp','mmol P m-3', 'POP')    
      call self%register_state_dependency(self%id_R6n,'RPn','mmol N m-3', 'PON')    
      call self%register_state_dependency(self%id_R6s,'RPs','mmol Si m-3','POSi')    

      ! Register links to external total dissolved inorganic carbon, dissolved oxygen pools
      call self%register_state_dependency(self%id_O3c,'O3c','mmol C m-3','Carbon Dioxide')    
      call self%register_state_dependency(self%id_O2o,'O2o','mmol O m-3','Oxygen')    

      ! Register environmental dependencies (temperature, shortwave radiation)
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      call self%register_dependency(self%id_eO2mO2,standard_variables%fractional_saturation_of_oxygen)

   end subroutine
   
   subroutine do(self,_ARGUMENTS_DO_)

      class (type_pml_ersem_microzoo),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

   ! !LOCAL VARIABLES:
      integer  :: iprey,istate
      real(rk) :: ETW, eO2mO2
      real(rk) :: Z5c,Z5p,Z5n,Z5cP,Z5nP,Z5pP
      real(rk),dimension(self%nprey) :: preycP,preynP,preypP,preysP,preyfP
      real(rk),dimension(self%nprey) :: spreyZ5,rupreyZ5c,fpreyZ5c
      real(rk) :: etZ5,CORROX,eO2Z5
      real(rk) :: rumZ5, put_uZ5,rugZ5
      real(rk) :: sdZ5,rdZ5
      real(rk) :: ineffZ5
      real(rk) :: retZ5,fZ5RDc,fZ5R6c
#ifdef SAVEFLX
      real(rk) :: fZXRDc
#endif
      real(rk) :: rrsZ5,rraZ5
      real(rk) :: fZ5O3c
      real(rk) :: fZ5RIp,fZ5RDp,fZ5R6p,fZ5N1p
      real(rk) :: fZ5RIn,fZ5RDn,fZ5R6n,fZ5NIn
      real(rk) :: qpZ5c,qnZ5c

      real(rk) :: preyP

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

      _GET_(self%id_ETW,ETW)
      _GET_(self%id_eO2mO2,eO2mO2)

      _GET_(self%id_c,Z5c)
      _GET_(self%id_p,Z5p)
      _GET_(self%id_n,Z5n)
      _GET_SAFE_(self%id_c,Z5cP)
      _GET_SAFE_(self%id_n,Z5nP)
      _GET_SAFE_(self%id_p,Z5pP)

      do iprey=1,self%nprey
         _GET_SAFE_(self%id_preyc(iprey),preycP(iprey))
         _GET_SAFE_(self%id_preyn(iprey),preynP(iprey))
         _GET_SAFE_(self%id_preyp(iprey),preypP(iprey))
         _GET_SAFE_(self%id_preys(iprey),preysP(iprey))
         _GET_SAFE_(self%id_preyf(iprey),preyfP(iprey))
      end do

      qpZ5c = Z5p/Z5c
      qnZ5c = Z5n/Z5c

!..Temperature effect :
      etZ5 = self%q10Z5X**((ETW-10._rk)/10._rk) - self%q10Z5X**((ETW-32._rk)/3._rk)

!..Oxygen limitation :
      CORROX = 1._rk + self%chrZ5oX
      eO2Z5 = MIN(1._rk,CORROX*(eO2mO2/(self%chrZ5oX + eO2mO2)))

!..Available food :
      spreyZ5 = self%suprey_Z5X*preycP/(preycP+self%minfoodZ5X)
      rupreyZ5c = spreyZ5*preycP
      rumZ5 = sum(rupreyZ5c)

!..Uptake :
      put_uZ5 = self%sumZ5X/(rumZ5 + self%chuZ5cX)*etZ5*Z5c
      rugZ5 = put_uZ5*rumZ5

!..Fluxes into microplankton :
      fpreyZ5c = put_uZ5*rupreyZ5c
      spreyZ5 = put_uZ5*spreyZ5

!..Zooplakton Grazing
!      Z5herb(I) = fP1Z5c(I) + fP2Z5c(I) + fP3Z5c(I) + fP4Z5c(I)
!      Z5carn(I) = fB1Z5c(I) + fZ6Z5c(I)

!..Mortality
      sdZ5 = ((1._rk - eO2Z5)*self%sdZ5oX + self%sdZ5X) 
      rdZ5 = sdZ5*Z5cP

!..Assimilation inefficiency:
      ineffZ5 = (1._rk - self%puZ5X)

!..Excretion
      retZ5 = ineffZ5 * rugZ5 * self%pu_eaZ5X
      fZ5RDc = (retZ5 + rdZ5)*self%pe_R1Z5X
      fZ5R6c = (retZ5 + rdZ5)*(1._rk - self%pe_R1Z5X)
#ifdef SAVEFLX
      fZXRDc = fZXRDc+fZ5RDc
#endif

!..Rest respiration, corrected for prevailing temperature
      rrsZ5 = self%srsZ5X*etZ5*Z5cP

!..Activity respiration
      rraZ5 = ineffZ5 * rugZ5 -retZ5

!..Total respiration
      fZ5O3c = rrsZ5 + rraZ5

#ifdef CALC
      fO3L2c(I) = fO3L2c(I) + & 
         RainR(I) * gutdiss * (1._rk - puZ5X)* pu_eaZ5X * fP2Z5c(I)
#endif
!..Source equation
      _SET_ODE_(self%id_c,rugZ5 - fZ5R6c - fZ5RDc - fZ5O3c)

!..Flows from and to detritus
      _SET_ODE_(self%id_R1c,(fZ5RDc * self%R1R2X))
      _SET_ODE_(self%id_R2c,(fZ5RDc * (1._rk-self%R1R2X)))
      _SET_ODE_(self%id_R6c,fZ5R6c)

!..Respiration
      _SET_ODE_(self%id_O3c,+ fZ5O3c/CMass)
      _SET_ODE_(self%id_O2o,- fZ5O3c*self%urB1_O2X)

!..Nutrient dynamics in microzooplankton, derived from carbon flows

!..Phosphorus dynamics
      fZ5RIp = retZ5*qpZ5c+sdZ5*Z5pP
      fZ5RDp = fZ5RIp*min(1._rk,self%pe_R1Z5X*self%xR1pX)
      fZ5R6p = (fZ5RIp - fZ5RDp)
      fZ5N1p = MAX( 0._rk, Z5pP - self%qpZ5cX*Z5cP)*self%stempZ5pX

!..Source equations
      _SET_ODE_(self%id_p,sum(spreyZ5*preypP) - fZ5R6p - fZ5RDp - fZ5N1p)

#ifdef IRON
! iron dynamics
! following Vichi et al., 2007 it is assumed that the iron fraction of the ingested phytoplankton
! is egested as particulate detritus (Luca)
      do iprey=1,self%nprey
         _SET_ODE_(self%id_preyf_target(iprey),+spreyZ5(iprey)*preyfP(iprey))
      end do
#endif

!..P-flow to pool
      _SET_ODE_(self%id_N1p,+ fZ5N1p)

!..Phosphorus flux from/to detritus
      _SET_ODE_(self%id_R6p,+ fZ5R6p)
      _SET_ODE_(self%id_R1p,+ fZ5RDp)

!..Nitrogen dynamics

      fZ5RIn = retZ5*qnZ5c+sdZ5*Z5nP
      fZ5RDn = fZ5RIn*min(1._rk,self%pe_R1Z5X*self%xR1nX)
      fZ5R6n = (fZ5RIn - fZ5RDn )

      fZ5NIn = MAX( 0._rk, Z5nP - self%qnZ5cX*Z5cP)*self%stempZ5nX
      _SET_ODE_(self%id_N4n,+ fZ5NIn)
      _SET_ODE_(self%id_n,sum(spreyZ5*preynP) - fZ5R6n - fZ5RDn - fZ5NIn)

!..Nitrogen flux from/to detritus
      _SET_ODE_(self%id_R6n,+ fZ5R6n)
      _SET_ODE_(self%id_R1n,+ fZ5RDn)

!..Silica-flux from diatoms due to microzooplankton grazing
      _SET_ODE_(self%id_R6s,sum(spreyZ5*preysP))

      ! Apply specific predation rates to all state variables of every prey.
      do iprey=1,self%nprey
         do istate=1,size(self%id_prey(iprey)%model%state)
            _GET_SAFE_(self%id_prey(iprey)%model%state(istate),preyP)
            _SET_ODE_(self%id_prey(iprey)%model%state(istate),-spreyZ5(iprey)*preyP)
         end do
      end do

      ! Leave spatial loops (if any)
      _LOOP_END_

   end subroutine do

end module