#include "fabm_driver.h"

!#define IRON

module pml_ersem_microzooplankton

   use fabm_types
   use fabm_particle

   use pml_ersem_shared
   use pml_ersem_pelagic_base

   implicit none

   private

   type,extends(type_ersem_pelagic_base_model),public :: type_pml_ersem_microzooplankton
      ! Variables
      type (type_model_id)                                   :: id_RP
      type (type_model_id),         allocatable,dimension(:) :: id_prey
      type (type_dependency_id),    allocatable,dimension(:) :: id_preyc,id_preyn,id_preyp,id_preys,id_preyf,id_preyl
      type (type_state_variable_id),allocatable,dimension(:) :: id_preyf_target
      type (type_state_variable_id)      :: id_O3c, id_O2o, id_L2c
      type (type_state_variable_id)      :: id_R1c, id_R2c, id_R6c
      type (type_state_variable_id)      :: id_R1p, id_R6p
      type (type_state_variable_id)      :: id_R1n, id_R6n
      type (type_state_variable_id)      :: id_R6s
      type (type_state_variable_id)      :: id_N1p,id_N4n
      type (type_dependency_id)          :: id_ETW,id_eO2mO2

      type (type_diagnostic_variable_id) :: id_fZ5O3c

      ! Parameters
      integer  :: nprey
      real(rk) :: chuz5cX,sumz5X,puz5X,pe_r1z5X,q10z5X,srsz5X,chrz5oX,pu_eaz5X
      real(rk) :: sdz5oX,sdz5X,qnz5cX,qpz5cX,minfoodz5X,stempz5nX,stempz5pX,gutdiss
      real(rk),allocatable :: suprey(:)

      ! ERSEM global parameters
      real(rk) :: R1R2X,urB1_O2X
   contains
!     Model procedures
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_pml_ersem_microzooplankton),intent(inout),target :: self
      integer,                        intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
      integer           :: iprey
      character(len=16) :: index
      real(rk)          :: c0
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%chuz5cX,   'chuc',   'mg C/m^3',   'Michaelis-Menten constant for food uptake')
      call self%get_parameter(self%sumz5X,    'sum',    '1/d',        'maximum specific uptake at reference temperature')
      call self%get_parameter(self%puz5X,     'pu',     '-',          'assimilation efficiency')
      call self%get_parameter(self%pe_r1z5X,  'pe_r1',  '-',          'dissolved fraction of excreted/dying matter')
      call self%get_parameter(self%q10z5X,    'q10',    '-',          'Q_10 temperature coefficient')
      call self%get_parameter(self%srsz5X,    'srs',    '1/d',        'specific rest respiration at reference temperature')
      call self%get_parameter(self%pu_eaz5X,  'pu_ea',  '-',          'fraction of unassimilated prey that is excreted (not respired)')
      call self%get_parameter(self%sdz5X,     'sd',     '1/d',        'basal mortality')
      call self%get_parameter(self%sdz5oX,    'sdo',    '1/d',        'maximum mortality due to oxygen limitation')
      call self%get_parameter(self%chrz5oX,   'chro',   '-',          'Michaelis-Menten constant for oxygen limitation')
      call self%get_parameter(self%qnz5cX,    'qnc',    'mmol N/mg C','maximum nitrogen to carbon ratio')
      call self%get_parameter(self%qpz5cX,    'qpc',    'mmol P/mg C','maximum phosphorus to carbon ratio')
      call self%get_parameter(self%minfoodz5X,'minfood','mg C/m^3',   'Michaelis-Menten constant to perceive food')
      call self%get_parameter(self%stempz5nX, 'stempn', '1/d',        'specific ammonium excretion rate')
      call self%get_parameter(self%stempz5pX, 'stempp', '1/d',        'specific phosphate excretion rate')

      call self%get_parameter(self%R1R2X,   'R1R2','-','labile fraction of produced DOM')
      call self%get_parameter(self%xR1pX,   'xR1p','-','transfer of phosphorus to DOM, relative to POM')
      call self%get_parameter(self%xR1nX,   'xR1n','-','transfer of nitrogen to DOM, relative to POM')
      call self%get_parameter(self%urB1_O2X,'urB1_O2','mmol O_2/mg C','oxygen consumed per carbon respired')

      call self%get_parameter(c0,'c0','mg C/m^3','background carbon concentration')

      call self%get_parameter(self%gutdiss,'gutdiss','-','fraction of prey calcite that is dissolved after ingestion')

      ! Register state variables
      call self%initialize_ersem_base(sedimentation=.false.)
      call self%add_constituent('c',1.e-4_rk,   c0)
      call self%add_constituent('n',1.26e-6_rk, qnRPIcX*c0)
      call self%add_constituent('p',4.288e-8_rk,qpRPIcX*c0)

      ! Register links to carbon contents of prey.
      call self%get_parameter(self%nprey,'nprey','','number of prey types',default=0)
      allocate(self%id_prey(self%nprey))
      allocate(self%id_preyc(self%nprey))
      allocate(self%id_preyn(self%nprey))
      allocate(self%id_preyp(self%nprey))
      allocate(self%id_preyf(self%nprey))
      allocate(self%id_preys(self%nprey))
      allocate(self%id_preyl(self%nprey))
      allocate(self%id_preyf_target(self%nprey))
      allocate(self%suprey(self%nprey))
      do iprey=1,self%nprey
         write (index,'(i0)') iprey
         call self%get_parameter(self%suprey(iprey),'suprey'//trim(index),'-','relative affinity for prey type '//trim(index))
         call self%register_dependency(self%id_preyc(iprey),'prey'//trim(index)//'c','mmol C m-3', 'Prey '//trim(index)//' C')
         call self%register_dependency(self%id_preyn(iprey),'prey'//trim(index)//'n','mmol N m-3', 'Prey '//trim(index)//' N')
         call self%register_dependency(self%id_preyp(iprey),'prey'//trim(index)//'p','mmol P m-3', 'Prey '//trim(index)//' P')
         call self%register_dependency(self%id_preys(iprey),'prey'//trim(index)//'s','mmol Si m-3','Prey '//trim(index)//' Si')
         call self%register_dependency(self%id_preyl(iprey),'prey'//trim(index)//'l','mg C m-3',   'Prey '//trim(index)//' calcite')
#ifdef IRON
         call self%register_dependency(self%id_preyf(iprey),'prey'//trim(index)//'f','mmol Fe m-3','Prey '//trim(index)//' Fe')
         call self%register_state_dependency(self%id_preyf_target(iprey),'prey'//trim(index)//'f_sink','umol Fe m-3','sink for Fe of prey '//trim(index),required=.false.)
#endif

         call self%register_model_dependency(self%id_prey(iprey),'prey'//trim(index))
         call self%request_coupling_to_model(self%id_preyc(iprey),self%id_prey(iprey),standard_variables%total_carbon)
         call self%request_coupling_to_model(self%id_preyn(iprey),self%id_prey(iprey),standard_variables%total_nitrogen)
         call self%request_coupling_to_model(self%id_preyp(iprey),self%id_prey(iprey),standard_variables%total_phosphorus)
         call self%request_coupling_to_model(self%id_preys(iprey),self%id_prey(iprey),standard_variables%total_silicate)
#ifdef IRON
         call self%request_coupling_to_model(self%id_preyf(iprey),self%id_prey(iprey),standard_variables%total_iron)
#endif
         call self%request_coupling_to_model(self%id_preyl(iprey),self%id_prey(iprey), &
                                             type_bulk_standard_variable(name='total_calcite_in_biota',aggregate_variable=.true.))
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

      call self%register_state_dependency(self%id_L2c,'L2c','mg C m-3','Calcite',required=.false.)

      ! Allow coupling of all required particulate organic matter variables to a single source model.
      call self%register_model_dependency(self%id_RP,'RP')
      call self%request_coupling_to_model(self%id_R6c,self%id_RP,'c')
      call self%request_coupling_to_model(self%id_R6n,self%id_RP,'n')
      call self%request_coupling_to_model(self%id_R6p,self%id_RP,'p')
      call self%request_coupling_to_model(self%id_R6s,self%id_RP,'s')

      ! Register diagnostics.
      call self%register_diagnostic_variable(self%id_fZ5O3c,'fZIO3c','mg C/m^3/d','respiration',output=output_time_step_averaged)

      ! Contribute to aggregate fluxes.
      call self%add_to_aggregate_variable(zooplankton_respiration_rate,self%id_fZ5O3c)

   end subroutine

   subroutine do(self,_ARGUMENTS_DO_)

      class (type_pml_ersem_microzooplankton),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

   ! !LOCAL VARIABLES:
      integer  :: iprey,istate
      real(rk) :: ETW, eO2mO2
      real(rk) :: Z5c,Z5p,Z5n,Z5cP,Z5nP,Z5pP
      real(rk),dimension(self%nprey) :: preycP,preynP,preypP,preysP,preylP
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
      eO2mO2 = min(1.0_rk,eO2mO2)

      _GET_WITH_BACKGROUND_(self%id_c,Z5c)
      _GET_WITH_BACKGROUND_(self%id_p,Z5p)
      _GET_WITH_BACKGROUND_(self%id_n,Z5n)
      _GET_(self%id_c,Z5cP)
      _GET_(self%id_n,Z5nP)
      _GET_(self%id_p,Z5pP)

      do iprey=1,self%nprey
         _GET_(self%id_preyc(iprey), preycP(iprey))
         _GET_(self%id_preyn(iprey), preynP(iprey))
         _GET_(self%id_preyp(iprey), preypP(iprey))
         _GET_(self%id_preys(iprey), preysP(iprey))
         _GET_(self%id_preyl(iprey), preylP(iprey))
      end do
      preycP = preycP*CMass

      qpZ5c = Z5p/Z5c
      qnZ5c = Z5n/Z5c

!..Temperature effect :
      etZ5 = self%q10Z5X**((ETW-10._rk)/10._rk) - self%q10Z5X**((ETW-32._rk)/3._rk)

!..Oxygen limitation :
      CORROX = 1._rk + self%chrZ5oX
      eO2Z5 = MIN(1._rk,CORROX*(eO2mO2/(self%chrZ5oX + eO2mO2)))

!..Available food :
      spreyZ5 = self%suprey*preycP/(preycP+self%minfoodZ5X)
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
      _SET_DIAGNOSTIC_(self%id_fZ5O3c,fZ5O3c)

      if (_AVAILABLE_(self%id_L2c)) then
         _SET_ODE_(self%id_L2c, (1.0_rk-self%gutdiss)*ineffZ5*self%pu_eaZ5X*sum(spreyZ5*preylP))
         _SET_ODE_(self%id_O3c,-(1.0_rk-self%gutdiss)*ineffZ5*self%pu_eaZ5X*sum(spreyZ5*preylP)/CMass)
      end if

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
         _GET_(self%id_preyf(iprey),preyP)
         if (preyP/=0.0_rk) _SET_ODE_(self%id_preyf_target(iprey),+spreyZ5(iprey)*preyP)
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
         do istate=1,size(self%id_prey(iprey)%state)
            _GET_(self%id_prey(iprey)%state(istate),preyP)
            _SET_ODE_(self%id_prey(iprey)%state(istate),-spreyZ5(iprey)*preyP)
         end do
      end do

      ! Leave spatial loops (if any)
      _LOOP_END_

   end subroutine do

end module