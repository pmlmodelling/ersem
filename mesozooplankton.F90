#include "fabm_driver.h"

module ersem_mesozooplankton

   use fabm_types
   use fabm_particle
   use fabm_expressions
   use fabm_builtin_models

   use ersem_shared
   use ersem_pelagic_base

   implicit none

   private

   type,extends(type_ersem_pelagic_base),public :: type_ersem_mesozooplankton
      ! Variables
      type (type_model_id),         allocatable,dimension(:) :: id_prey
      type (type_model_id)                                   :: id_RP
      type (type_dependency_id),    allocatable,dimension(:) :: id_preyc,id_preyn,id_preyp,id_preys,id_preyf,id_preyl
      type (type_state_variable_id),allocatable,dimension(:) :: id_preyf_target
      type (type_state_variable_id)      :: id_O3c,id_O2o,id_L2c
      type (type_state_variable_id)      :: id_R1c,id_R1p,id_R1n
      type (type_state_variable_id)      :: id_R2c
      type (type_state_variable_id)      :: id_R8c,id_R8p,id_R8n,id_R8s
      type (type_state_variable_id)      :: id_N1p,id_N4n
      type (type_dependency_id)          :: id_ETW,id_eO2mO2,id_totprey
      type (type_horizontal_dependency_id) :: id_inttotprey

      type (type_diagnostic_variable_id) :: id_fZ4O3c

      ! Parameters
      integer  :: nprey
      real(rk) :: qpZIcX,qnZIcX
      real(rk) :: q10Z4X,chrZ4oX,minfoodZ4X,chuZ4cX
      real(rk),allocatable :: suprey(:),pu_eaZ4X(:)
      real(rk) :: sumZ4X
      real(rk) :: sdZ4oX, sdZ4X, srsZ4X
      real(rk) :: puZ4X
      real(rk) :: pe_R1Z4X
      real(rk) :: gutdiss

      real(rk) :: MinpreyX,Z4repwX,Z4mortX

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
      class (type_ersem_mesozooplankton),intent(inout),target :: self
      integer,                       intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
      integer           :: iprey
      character(len=16) :: index
      type (type_weighted_sum),pointer :: total_prey_calculator
      real(rk)          :: pu_eaZ4X,pu_eaRZ4X
      logical           :: preyispom
      real(rk)          :: c0
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%qpZIcX,    'qpc',    'mmol P/mg C','phosphorous to carbon ratio')
      call self%get_parameter(self%qnZIcX,    'qnc',    'mmol N/mg C','nitrogen to carbon ratio')
      call self%get_parameter(self%q10z4X,    'q10',    '-','Q_10 temperature coefficient')
      call self%get_parameter(self%chrZ4oX,   'chro',   '','Michaelis Menten constant for oxygen limitation')
      call self%get_parameter(self%minfoodZ4X,'minfood','mg C/m^3','Michaelis-Menten constant to perceive food')
      call self%get_parameter(self%chuZ4cX,   'chuc',   'mg C/m^3','Michaelis Menten constant for food uptake')
      call self%get_parameter(self%sumZ4X,    'sum',    '1/d','maximum specific uptake at reference temperature')
      call self%get_parameter(self%sdZ4oX,    'sdo',    '1/d','specific mortality due to oxygen limitation')
      call self%get_parameter(self%sdZ4X,     'sd',     '1/d','specific basal mortality')
      call self%get_parameter(self%srsZ4X,    'srs',    '1/d','specific rest respiration at reference temperature')
      call self%get_parameter(self%puZ4X,     'pu',     '-', 'relative assimilation efficiency')
      call self%get_parameter(pu_eaZ4X,       'pu_ea',  '-','excreted fraction of prey uptake')
      call self%get_parameter(pu_eaRZ4X,      'pu_eaR', '-','excreted fraction of POM uptake')
      call self%get_parameter(self%pe_R1Z4X,  'pe_R1',  '-','DOM fraction of excreted matter')
      call self%get_parameter(self%MinpreyX,  'Minprey','mg C/m^2','food treshold for overwintering state')
      call self%get_parameter(self%Z4repwX,   'repw',   '1/d','specific overwintering respiration')
      call self%get_parameter(self%Z4mortX,   'mort',   '1/d','specific overwintering mortality')
      call self%get_parameter(c0,             'c0',     'mg C/m^3','background concentration')
      call self%get_parameter(self%gutdiss,   'gutdiss','-','fraction of prey calcite that is dissolved after ingestion')

      call self%get_parameter(self%R1R2X,   'R1R2','-','labile fraction of produced DOM')
      call self%get_parameter(self%xR1pX,   'xR1p','-','transfer of phosphorus to DOM, relative to POM')
      call self%get_parameter(self%xR1nX,   'xR1n','-','transfer of nitrogen to DOM, relative to POM')
      call self%get_parameter(self%urB1_O2X,'urB1_O2','mmol O_2/mg C','oxygen consumed per carbon respired')

      ! Register state variables
      call self%initialize_ersem_base(sedimentation=.false.)
      call self%add_constituent('c',1.e-4_rk,c0,qn=self%qnZIcX,qp=self%qpZIcX)

      ! Create an expression that will compute the total prey
      ! (will be depth integrated to determine overwintering)
      allocate(total_prey_calculator)
      total_prey_calculator%output_units = 'mg C m-3'

      ! Determine number of prey types.
      call self%get_parameter(self%nprey,'nprey','','number of prey types',default=0)

      ! Get prey-specific parameters.
      allocate(self%suprey(self%nprey))
      allocate(self%pu_eaZ4X(self%nprey))
      do iprey=1,self%nprey
         write (index,'(i0)') iprey
         call self%get_parameter(self%suprey(iprey),'suprey'//trim(index),'-','relative affinity for prey type '//trim(index))
      end do
      do iprey=1,self%nprey
         write (index,'(i0)') iprey
         call self%get_parameter(preyispom,'prey'//trim(index)//'ispom','','prey type '//trim(index)//' is POM',default=.false.)
         if (preyispom) then
            self%pu_eaZ4X(iprey) = pu_eaRZ4X
         else
            self%pu_eaZ4X(iprey) = pu_eaZ4X
         end if
      end do

      ! Get prey-specific coupling links.
      allocate(self%id_prey(self%nprey))
      allocate(self%id_preyc(self%nprey))
      allocate(self%id_preyn(self%nprey))
      allocate(self%id_preyp(self%nprey))
      allocate(self%id_preyf(self%nprey))
      allocate(self%id_preys(self%nprey))
      allocate(self%id_preyl(self%nprey))
      allocate(self%id_preyf_target(self%nprey))
      do iprey=1,self%nprey
         write (index,'(i0)') iprey
         call self%register_dependency(self%id_preyc(iprey), 'prey'//trim(index)//'c','mmol C m-3', 'Prey '//trim(index)//' C')    
         call self%register_dependency(self%id_preyn(iprey), 'prey'//trim(index)//'n','mmol N m-3', 'Prey '//trim(index)//' N')    
         call self%register_dependency(self%id_preyp(iprey), 'prey'//trim(index)//'p','mmol P m-3', 'Prey '//trim(index)//' P')    
         call self%register_dependency(self%id_preys(iprey), 'prey'//trim(index)//'s','mmol Si m-3','Prey '//trim(index)//' Si')
         call self%register_dependency(self%id_preyl(iprey), 'prey'//trim(index)//'l','mg C m-3',   'Prey '//trim(index)//' calcite')

         call self%register_model_dependency(self%id_prey(iprey),'prey'//trim(index))
         call self%request_coupling_to_model(self%id_preyc(iprey),self%id_prey(iprey),standard_variables%total_carbon)
         call self%request_coupling_to_model(self%id_preyn(iprey),self%id_prey(iprey),standard_variables%total_nitrogen)
         call self%request_coupling_to_model(self%id_preyp(iprey),self%id_prey(iprey),standard_variables%total_phosphorus)
         call self%request_coupling_to_model(self%id_preys(iprey),self%id_prey(iprey),standard_variables%total_silicate)
         call self%request_coupling_to_model(self%id_preyl(iprey),self%id_prey(iprey), &
                                             type_bulk_standard_variable(name='total_calcite_in_biota',aggregate_variable=.true.))

         if (use_iron) then
            call self%register_dependency(self%id_preyf(iprey), 'prey'//trim(index)//'f','mmol Fe m-3','Prey '//trim(index)//' Fe')    
            call self%register_state_dependency(self%id_preyf_target(iprey),'prey'//trim(index)//'f_sink','umol Fe m-3','sink for Fe of prey '//trim(index),required=.false.)    
            call self%request_coupling_to_model(self%id_preyf(iprey),self%id_prey(iprey),standard_variables%total_iron)
         end if

         call total_prey_calculator%add_component('prey'//trim(index)//'c',self%suprey(iprey))
      end do

      ! Add the submodel that will compute total prey for us, and create a variable that will contain its depth integral.
      call self%add_child(total_prey_calculator,'totprey_calculator',configunit=-1)
      call self%register_dependency(self%id_totprey,'totprey')
      call self%request_coupling(self%id_totprey,'totprey_calculator/result')
      call self%register_expression_dependency(self%id_inttotprey,vertical_integral(self%id_totprey))

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
      call self%register_state_dependency(self%id_R8c,'RPc','mg C m-3',   'POC')    
      call self%register_state_dependency(self%id_R8p,'RPp','mmol P m-3', 'POP')    
      call self%register_state_dependency(self%id_R8n,'RPn','mmol N m-3', 'PON')    
      call self%register_state_dependency(self%id_R8s,'RPs','mmol Si m-3','POSi')    

      ! Allow coupling of all required particulate organic matter variables to a single source model.
      call self%register_model_dependency(self%id_RP,'RP')
      call self%request_coupling_to_model(self%id_R8c,self%id_RP,'c')
      call self%request_coupling_to_model(self%id_R8n,self%id_RP,'n')
      call self%request_coupling_to_model(self%id_R8p,self%id_RP,'p')
      call self%request_coupling_to_model(self%id_R8s,self%id_RP,'s')

      ! Register links to external total dissolved inorganic carbon, dissolved oxygen pools
      call self%register_state_dependency(self%id_O3c,'O3c','mmol C m-3','Carbon Dioxide')    
      call self%register_state_dependency(self%id_O2o,'O2o','mmol O m-3','Oxygen')    

      call self%register_state_dependency(self%id_L2c,'L2c','mg C m-3','Calcite',required=.false.)

      ! Register environmental dependencies (temperature, shortwave radiation)
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      call self%register_dependency(self%id_eO2mO2,standard_variables%fractional_saturation_of_oxygen)

      ! Register diagnostics
      call self%register_diagnostic_variable(self%id_fZ4O3c,'fZIO3c','mg C/m^3/d','respiration',output=output_time_step_averaged)

      ! Contribute to aggregate fluxes.
      call self%add_to_aggregate_variable(zooplankton_respiration_rate,self%id_fZ4O3c)

   end subroutine
   
   subroutine do(self,_ARGUMENTS_DO_)

      class (type_ersem_mesozooplankton),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

   ! !LOCAL VARIABLES:
      integer  :: iprey,istate
      real(rk) :: ETW,eO2mO2
      real(rk) :: Z4c,Z4cP
      real(rk),dimension(self%nprey) :: preycP,preypP,preynP,preysP,preylP
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
      real(rk) :: preyP
      real(rk) :: intprey

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         _GET_HORIZONTAL_(self%id_inttotprey,intprey)
         intprey = intprey*CMass

         if (intprey>self%MinpreyX) then
            ! Enough prey available - not overwintering

            ! Get environment (temperature, oxygen saturation)
            _GET_(self%id_ETW,ETW)
            _GET_(self%id_eO2mO2,eO2mO2)
            eO2mO2 = min(1.0_rk,eO2mO2)

            ! Get own concentrations
            _GET_WITH_BACKGROUND_(self%id_c,Z4c)
            _GET_(self%id_c,Z4cP)

            ! Get prey concentrations
            do iprey=1,self%nprey
               _GET_(self%id_preyc(iprey), preycP(iprey))
               _GET_(self%id_preyp(iprey), preypP(iprey))
               _GET_(self%id_preyn(iprey), preynP(iprey))
               _GET_(self%id_preys(iprey), preysP(iprey))
               _GET_(self%id_preyl(iprey), preylP(iprey))
            end do
            preycP = preycP*CMass

   !..Temperature effect :

            etZ4 = self%q10Z4X**((ETW-10._rk)/10._rk) - self%q10Z4X**((ETW-32._rk)/3._rk)

   !..Oxygen limitation :
            CORROX = 1._rk + self%chrZ4oX
            eO2Z4 = MIN(1._rk,CORROX*(eO2mO2/(self%chrZ4oX + eO2mO2)))

   !..Available food :
            spreyZ4 = self%suprey*preycP/(preycP+self%minfoodZ4X)
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
            retZ4 = ineffZ4 * sum(fpreyZ4c*self%pu_eaZ4X)
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

            if (_AVAILABLE_(self%id_L2c)) then
               _SET_ODE_(self%id_L2c, (1.0_rk-self%gutdiss)*ineffZ4*sum(self%pu_eaZ4X*spreyZ4*preylP))
               _SET_ODE_(self%id_O3c,-(1.0_rk-self%gutdiss)*ineffZ4*sum(self%pu_eaZ4X*spreyZ4*preylP)/CMass)
            end if
               
            rraZ4 = rugZ4*(1._rk - self%puZ4X)-retZ4

   !..Source equations
            SZ4c = rugZ4 - fZ4RDc - fZ4R8c - fZ4O3c   ! Jorn: cannibalism accounted for in grazing/predation section

   !..Flows from and to detritus
            _SET_ODE_(self%id_R1c, + fZ4RDc * self%R1R2X)
            _SET_ODE_(self%id_R2c, + fZ4RDc * (1._rk-self%R1R2X))
            _SET_ODE_(self%id_R8c, + fZ4R8c)

   !..Respiration
            _SET_ODE_(self%id_O3c, + fZ4O3c/CMass)
            _SET_ODE_(self%id_O2o, - fZ4O3c*self%urB1_O2X)

   !..Phosphorus dynamics in mesozooplankton, derived from carbon flows....

            fZ4RIp = (fZ4RDc + fZ4R8c) * self%qpZIcX
            fZ4RDp = min(fZ4RIp, fZ4RDc * self%qpZIcX * self%xR1pX)
            fZ4R8p = fZ4RIp - fZ4RDp

   !..Source equations
            SZ4p = sum(spreyZ4*preypP) - fZ4R8p - fZ4RDp

            if (use_iron) then
   !  Iron dynamics
   ! following Vichi et al., 2007 it is assumed that the iron fraction of the ingested phytoplankton
   ! is egested as particulate detritus (Luca)
               do iprey=1,self%nprey
                  _GET_(self%id_preyf(iprey),preyP)
                  if (preyP/=0.0_rk) _SET_ODE_(self%id_preyf_target(iprey),spreyZ4(iprey)*preyP)
               end do
            end if

   !..Phosphorus flux from/to detritus
            _SET_ODE_(self%id_R1p, + fZ4RDp)
            _SET_ODE_(self%id_R8p, + fZ4R8p)

   !..Nitrogen dynamics in mesozooplankton, derived from carbon flows......

            fZ4RIn = (fZ4RDc + fZ4R8c) * self%qnZIcX
            fZ4RDn = min(fZ4RIn, fZ4RDc * self%qnZIcX * self%xR1nX)
            fZ4R8n = fZ4RIn - fZ4RDn

   !..Source equations
            SZ4n = sum(spreyZ4*preynP) - fZ4R8n - fZ4RDn

   !..Nitrogen flux from/to detritus
            _SET_ODE_(self%id_R1n, + fZ4RDn)
            _SET_ODE_(self%id_R8n, + fZ4R8n)

   !..Silica-flux from diatoms due to mesozooplankton grazing
            _SET_ODE_(self%id_R8s,sum(spreyZ4*preysP))

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

            ! Apply specific predation rates to all state variables of every prey.
            do iprey=1,self%nprey
               do istate=1,size(self%id_prey(iprey)%state)
                  _GET_(self%id_prey(iprey)%state(istate),preyP)
                  _SET_ODE_(self%id_prey(iprey)%state(istate),-spreyZ4(iprey)*preyP)
               end do
            end do
         else
            ! Insufficient prey - overwintering
            _GET_(self%id_c,Z4cP)

!.. Respiration
            fZ4O3c = Z4cP * self%Z4repwX

!.. Mortality
            fZ4R8c = Z4cP * self%Z4mortX

            _SET_ODE_(self%id_R8c,fZ4R8c)
            _SET_ODE_(self%id_R8n,fZ4R8c*self%qnZIcX)
            _SET_ODE_(self%id_R8p,fZ4R8c*self%qpZIcX)

            _SET_ODE_(self%id_O3c,fZ4O3c/CMass)
            _SET_ODE_(self%id_N4n,fZ4O3c*self%qnZIcX)
            _SET_ODE_(self%id_N1p,fZ4O3c*self%qpZIcX)

            _SET_ODE_(self%id_c,- fZ4R8c - fZ4O3c)

         end if

         _SET_DIAGNOSTIC_(self%id_fZ4O3c,fZ4O3c)

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

       IF ( ExcessN .GT. 0.0_rk ) THEN ! Jorn: 0.0_rk was ZeroX before (why?), but that cause violation of conservation
         SXn = SXn - ExcessN
         SKn = SKn + ExcessN
       END If

       IF ( ExcessP .GT. 0.0_rk ) THEN ! Jorn: 0.0_rk was ZeroX before (why?), but that cause violation of conservation
         SXp = SXp - EXcessP
         SKp = SKp + ExcessP
       END If

       END SUBROUTINE Adjust_fixed_nutrients
!
!EOC
!-----------------------------------------------------------------------
   
end module