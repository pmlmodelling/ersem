#include "fabm_driver.h"

! -----------------------------------------------------------------------------
! This is a generic model for primary producers - P1, P2, P3, P4 in ERSEM.
! -----------------------------------------------------------------------------
! Silicate use (P1) is optional - it is activated by setting the "use_Si" flag 
! in the run-time configuration.
!
! Calcification (P2) is optional - it is activated by setting the "calcify"
! flag in the run-time configuration. Preprocessor symbol CALC is no longer
! used.
!
! The model can be switched from a constant ratio between labile and semi-
! labile dissolved organic matter production to a dynamic ratio by setting the
! "docdyn" flag in the run-time configuration. Preprocessor symbol DOCDYN
! is no longer used.
!
! Preprocessor symbols IRON, CENH are still used, but should be
! replaced by run-time flags provided performance does not suffer.
! -----------------------------------------------------------------------------

module pml_ersem_primary_producer

   use fabm_types
   use fabm_particle

   use pml_ersem_shared
   use pml_ersem_pelagic_base

   implicit none

   private

   type,extends(type_ersem_pelagic_base_model),public :: type_pml_ersem_primary_producer
      ! Identifiers for links to state variables of other models
      type (type_state_variable_id) :: id_O3c,id_O2o                        ! Oxygen and DIC
      type (type_state_variable_id) :: id_N5s,id_N1p,id_N3n,id_N4n,id_N7f   ! Nutrients
      type (type_state_variable_id) :: id_R1c,id_R1p,id_R1n,id_R2c          ! DOM
      type (type_state_variable_id) :: id_R6c,id_R6p,id_R6n,id_R6s,id_R6f   ! POM
      type (type_state_variable_id) :: id_L2c                               ! Free calcite (liths) - used by calcifiers only

      ! Environmental dependencies
      type (type_dependency_id)            :: id_parEIR,id_ETW   ! PAR and temperature
      type (type_dependency_id)            :: id_RainR           ! Rain ratio - used by calcifiers only
      type (type_horizontal_dependency_id) :: id_pco2a3          ! Atmospheric pCO2 - used only if CENH is active

      ! Identifiers for diagnostic variables
      type (type_diagnostic_variable_id) :: id_netP1   ! Net primary production
      type (type_diagnostic_variable_id) :: id_lD      ! Cell-bound calcite - used by calcifiers only

      ! Identifiers for coupled models
      type (type_model_id) :: id_RP   ! Model that provides all POM constituents

      ! Parameters (described in subroutine initialize, below)
      real(rk) :: sump1X
      real(rk) :: q10p1X,srsp1X,pu_eap1X,pu_rap1X,chp1sX,qnlp1cX,qplp1cX,xqcp1pX
      real(rk) :: xqcp1nX,xqnp1X,xqpp1X,qup1n3X,qup1n4X,qurp1pX,qsp1cX,esnip1X
      real(rk) :: resp1mX,sdop1X
      real(rk) :: alphaP1X,betaP1X,phimP1X,phiP1HX
      real(rk) :: R1R2X,uB1c_O2X,urB1_O2X
      real(rk) :: qflP1cX,qfRP1cX,qurP1fX
      real(rk) :: rP1mX
      integer :: LimnutX
      logical :: use_Si, calcify, docdyn

   contains

      ! Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: get_vertical_movement
      procedure :: get_sinking_rate

   end type type_pml_ersem_primary_producer

   ! Constants
   real(rk),parameter :: ChlCmin = 0.0067_rk

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_pml_ersem_primary_producer),intent(inout),target :: self
      integer,                                intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
      real(rk) :: c0,EPS
!EOP
!-----------------------------------------------------------------------
!BOC
      ! Obtain the values of all model parameters from FABM.
      ! Specify the long name and units of the parameters, which could be used by FABM (or its host)
      ! to present parameters to the user for configuration (e.g., through a GUI)
      call self%get_parameter(self%use_Si,   'use_Si', default=.true.)
      call self%get_parameter(self%sump1X,   'sum',  '1/d',        'Specific maximal productivity at reference temperature')
      call self%get_parameter(self%q10p1X,   'q10',  '-',          'Regulating temperature factor Q10')
      call self%get_parameter(self%srsp1X,   'srs',  '1/d',        'Specific rest respiration at reference temperature')
      call self%get_parameter(self%pu_eap1X, 'pu_ea','-',          'Excreted fraction of primary production')
      call self%get_parameter(self%pu_rap1X, 'pu_ra','-',          'Respired fraction of primary production')
      call self%get_parameter(self%qnlp1cX,  'qnlc', 'mmol N/mg C','Minimal nitrogen to carbon ratio')
      call self%get_parameter(self%qplp1cX,  'qplc', 'mmol P/mg C','Minimal phosphorus to carbon ratio')
      call self%get_parameter(self%xqcp1pX,  'xqcp', '-',          'Threshold for phosphorus limitation (relative to Redfield ratio)')
      call self%get_parameter(self%xqcp1nX,  'xqcn', '-',          'Threshold for nitrogen limitation (relative to Redfield ratio)')
      call self%get_parameter(self%xqpp1X,   'xqp',  '-',          'Maximal phosphorus to carbon ratio (relative to Redfield ratio)')
      call self%get_parameter(self%xqnp1X,   'xqn',  '-',          'Maximal nitrogen to carbon ratio (relative to Redfield ratio)')
      call self%get_parameter(self%qup1n3X,  'qun3', 'm^3/mg C/d', 'Nitrate affinity')
      call self%get_parameter(self%qup1n4X,  'qun4', 'm^3/mg C/d', 'Ammonium affinity')
      call self%get_parameter(self%qurp1pX,  'qurp', 'm^3/mg C/d', 'Phosphate affinity')
      if (self%use_Si) then
         call self%get_parameter(self%qsp1cX,'qsc', 'mmol Si/mg C','Maximal silicon to carbon ratio')
         call self%get_parameter(self%chp1sX,'chs', 'mmol/m^3',    'Michaelis-Menten constant for silicate limitation')
      end if
      call self%get_parameter(self%esnip1X,  'esni', '-',          'Level of nutrient limitation at which sinking commences')
      call self%get_parameter(self%resp1mX,  'resm', 'm/d',        'Maximal sinking velocity')
      call self%get_parameter(self%sdop1X,   'sdo',  '1/d',        '1.1 of minimal specific lysis rate')
      call self%get_parameter(self%alphaP1X, 'alpha','mg C m^2/mg Chl/W/d', 'Initial slope of PI-curve')
      call self%get_parameter(self%betaP1X,  'beta', 'mg C m^2 /mg Chl/W/d','Photo-inhibition parameter')
      call self%get_parameter(self%phimP1X,  'phim', 'mg Chl/mg C','Maximal effective chlorophyll to carbon photosynthesis ratio')
      call self%get_parameter(self%phiP1HX,  'phiH', 'mg Chl/mg C','Minimal effective chlorophyll to carbon photosynthesis ratio')
      call self%get_parameter(self%LimnutX,  'Limnut','',          'Nitrogen-phosphorus colimitation formulation')
      call self%get_parameter(self%docdyn,   'docdyn','','use dynamic ratio of labile to semi-labile DOM production', default=.false.)
      if (.not.self%docdyn) call self%get_parameter(self%R1R2X,'R1R2','-','Labile fraction of DOM production')
      call self%get_parameter(self%uB1c_O2X, 'uB1c_O2','mg C/mmol O2','Conversion of carbon into oxygen produced')
      call self%get_parameter(self%urB1_O2X, 'urB1_O2','mg C/mmol O2','Conversion of carbon into oxygen respired')
#ifdef IRON
      call self%get_parameter(self%qflP1cX,  'qflc','umol Fe/mg C','Minimal iron to carbon ratio')
      call self%get_parameter(self%qfRP1cX,  'qfRc','umol Fe/mg C','Maximal/optimal iron to carbon ratio')
      call self%get_parameter(self%qurP1fX,  'qurf','m^3/mg C/d',  'Specific affinity for iron')
#endif
      call self%get_parameter(EPS,           'EPS', 'm"2/mg C', 'specific extinction coefficient')
      call self%get_parameter(c0,            'c0',  'mg C/m^3', 'background concentration', default=0.0_rk)
      call self%get_parameter(self%calcify,  'calcify',default=.false.)
      call self%get_parameter(self%rP1mX,    'rm',  'm/d',      'background sinking velocity',default=0.0_rk)

      ! Register state variables
      call self%initialize_ersem_base(sedimentation=.true.)
      call self%add_constituent('c',1.e-4_rk,   c0)
      call self%add_constituent('n',1.26e-6_rk, c0*qnrpicX)
      call self%add_constituent('p',4.288e-8_rk,c0*qprpicX)
      call self%add_constituent('f',5.e-6_rk,   0.0_rk)  ! NB this does nothing if iron support is disabled.
      call self%add_constituent('chl',3.e-6_rk, c0*self%phimP1X)
      if (self%use_Si) call self%add_constituent('s',1.e-6_rk,c0*self%qsp1cX)

      ! Register links to external nutrient pools.
      call self%register_state_dependency(self%id_N1p,'N1p','mmol P/m^3', 'Phosphate')    
      call self%register_state_dependency(self%id_N3n,'N3n','mmol N/m^3', 'Nitrate')    
      call self%register_state_dependency(self%id_N4n,'N4n','mmol N/m^3', 'Ammonium')    
      if (self%use_Si) call self%register_state_dependency(self%id_N5s,'N5s','mmol Si/m^3','Silicate')    
#ifdef IRON
      call self%register_state_dependency(self%id_N7f,'N7f','umol Fe/m^3', 'Inorganic Iron')    
#endif

      ! Register links to external labile dissolved organic matter pools.
      call self%register_state_dependency(self%id_R1c,'R1c','mg C/m^3',  'DOC')
      call self%register_state_dependency(self%id_R1p,'R1p','mmol P/m^3','DOP')    
      call self%register_state_dependency(self%id_R1n,'R1n','mmol N/m^3','DON')    

      ! Register links to external semi-labile dissolved organic matter pool.
      call self%register_state_dependency(self%id_R2c,'R2c','mg C/m^3','Semi-labile DOC')    

      ! Register links to external particulate organic matter pools.
      ! At run-time, these can be coupled to any available pools (e.g., R4, R6, R8 in ERSEM)
      call self%register_state_dependency(self%id_R6c,'RPc','mg C/m^3',   'POC')    
      call self%register_state_dependency(self%id_R6p,'RPp','mmol P/m^3', 'POP')    
      call self%register_state_dependency(self%id_R6n,'RPn','mmol N/m^3', 'PON')    
      if (self%use_Si) call self%register_state_dependency(self%id_R6s,'RPs','mmol Si/m^3','POSi')    
#ifdef IRON
      call self%register_state_dependency(self%id_R6f,'RPf','umol Fe/m^3','POFe')    
#endif

      ! Automatically hook up all components of external particulate organic matter,
      ! by obtaining them from a single named model "RP". This takes away the need to couple each RP?
      ! variable individually.
      call self%register_model_dependency(self%id_RP,'RP')
      call self%request_coupling_to_model(self%id_R6c,self%id_RP,'c')
      call self%request_coupling_to_model(self%id_R6n,self%id_RP,'n')
      call self%request_coupling_to_model(self%id_R6p,self%id_RP,'p')
      if (_VARIABLE_REGISTERED_(self%id_R6s)) call self%request_coupling_to_model(self%id_R6s,self%id_RP,'s')
      if (_VARIABLE_REGISTERED_(self%id_R6f)) call self%request_coupling_to_model(self%id_R6f,self%id_RP,'f')

      ! Register links to external total dissolved inorganic carbon, dissolved oxygen pools
      call self%register_state_dependency(self%id_O3c,'O3c','mmol C/m^3','Carbon Dioxide')    
      call self%register_state_dependency(self%id_O2o,'O2o','mmol O/m^3','Oxygen')    

      ! Register diagnostic variables (i.e., model outputs)
      call self%register_diagnostic_variable(self%id_netP1,'netP1','mg C/m^3/d','net primary production',output=output_time_step_averaged)

      ! Register environmental dependencies (temperature, shortwave radiation)
      call self%register_dependency(self%id_parEIR,standard_variables%downwelling_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      if (self%calcify) then
         ! Link to rain ratio (set by calcification module; excludes effect of nutrient limitation, which is computed below)
         call self%register_dependency(self%id_RainR,'RainR','-','Rain ratio')

         ! Link to external detached liths (sink for calcite of dying phytoplankton)
         call self%register_state_dependency(self%id_L2c,'L2c','mg C/m^3','Free calcite')

         ! Create diagnostic variable for cell-bound calcite, and register its contribution to the known quantity "total_calcite_in_biota".
         ! This quantity can be used by other models (e.g., predators) to determine how much calcite is released when phytoplankton is broken down.
         call self%register_diagnostic_variable(self%id_lD,'l','mg C m-3','Bound calcite',missing_value=0._rk,output=output_none)
         call self%add_to_aggregate_variable(type_bulk_standard_variable(name='total_calcite_in_biota',aggregate_variable=.true.),self%id_lD)
      end if
#ifdef CENH
      ! Link to atmospheric CO2 (only if using CO2-enhanced primary production).
      call self%register_dependency(self%id_pco2a3,standard_variables%mole_fraction_of_carbon_dioxide_in_air)    
#endif

      ! Register contribution to light extinction
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &
         self%id_c,scale_factor=EPS,include_background=.true.)

   end subroutine

   subroutine do(self,_ARGUMENTS_DO_)

      class (type_pml_ersem_primary_producer),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

   ! !LOCAL VARIABLES:
      real(rk) :: ETW,parEIR
      real(rk) :: P1c, P1p, P1n, Chl1
      real(rk) :: P1cP,P1pP,P1nP,P1sP,Chl1P
      real(rk) :: N5s,N1pP,N3nP,N4nP
      real(rk) :: iNP1n,iNP1p,iNP1s,iNP1f,iNIP1
      real(rk) :: qpP1c,qnP1c

      real(rk) :: srsP1
      real(rk) :: fP1R6c,fP1RDc
      real(rk) :: fP1O3c,fO3P1c
      real(rk) :: fP1R6p,fP1RDp,fN1P1p
      real(rk) :: fP1R6n,fP1RDn
      real(rk) :: fNIP1n,fN3P1n,fN4P1n
      real(rk) :: fP1R6s,fP1N5s,fN5P1s

      real(rk) :: sdoP1,sumP1,sugP1,seoP1,seaP1,sraP1,sunP1,runP1,rugP1,sP1R6
      real(rk) :: runP1p,misP1p,rumP1p
      real(rk) :: runP1n,misP1n,rumP1n,rumP1n3,rumP1n4
      real(rk) :: etP1,pe_R6P1
      real(rk) :: rho,Chl_inc,Chl_loss
      real(rk) :: phi,ChlCpp
#ifdef IRON
      real(rk) :: N7fP,P1f,P1fP,qfP1c
      real(rk) :: runP1f,rumP1f,misP1f
      real(rk) :: fN7P1f,fP1R6f
#endif
      real(rk) :: fP1R1c,fP1R2c
#ifdef CENH
      real(rk) :: pco2a3,cenh
#endif
      real(rk) :: RainR, t

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve local state (carbon, phosphorus, nitrogen, chlorophyll concentrations).

         ! Concentrations including background (used in source terms)
         _GET_WITH_BACKGROUND_(self%id_c,P1c)
         _GET_WITH_BACKGROUND_(self%id_p,P1p)
         _GET_WITH_BACKGROUND_(self%id_n,P1n)
         _GET_WITH_BACKGROUND_(self%id_chl,chl1)

         ! Concentrations excluding background (used in sink terms)
         _GET_(self%id_c,P1cP)
         _GET_(self%id_p,P1pP)
         _GET_(self%id_n,P1nP)
         _GET_(self%id_chl,chl1P)

         ! Retrieve ambient nutrient concentrations
         _GET_(self%id_N1p,N1pP)
         _GET_(self%id_N3n,N3nP)
         _GET_(self%id_N4n,N4nP)

         ! Retrieve environmental dependencies (water temperature, photosynthetically active radation)
         _GET_(self%id_ETW,ETW)
         _GET_(self%id_parEIR,parEIR)

         qpP1c = P1p/P1c
         qnP1c = P1n/P1c
            
         ! Regulation factors...................................................

         iNP1p = MIN(1._rk,  &
                 MAX(0._rk, (qpP1c-self%qplP1cX) / (self%xqcP1pX*qpRPIcX-self%qplP1cX) ))
         iNP1n = MIN(1._rk,  &
                 MAX(0._rk, (qnP1c-self%qnlP1cX) / (self%xqcP1nX*qnRPIcX-self%qnlP1cX) ))

#ifdef IRON
         _GET_WITH_BACKGROUND_(self%id_f,P1f)
         qfP1c = P1f/P1c
         iNP1f = MIN(1._rk,  &
                 MAX(ZeroX, (qfP1c-self%qflP1cX) / (self%qfRP1cX-self%qflP1cX) ))
#else
         iNP1f = 1.0_rk
#endif

         if (self%use_Si) then
            _GET_WITH_BACKGROUND_(self%id_N5s,N5s) ! Jorn: kept identical to legacy ersem, but should likely exclude background
            iNP1s = MIN(1._rk, N5s/(N5s+self%chP1sX))
         else
            iNP1s = 1.0_rk
         end if

         if (self%LimnutX==0) then  ! Jorn: select case would be cleaner but makes vectorization impossible for ifort 14
            iNIP1 = (iNP1p * iNP1n)**0.5_rk
         elseif (self%LimnutX==1) then
            iNIP1 = MIN(iNP1p, iNP1n)
         else
            iNIP1 = 2.0_rk / (1._rk/iNP1p + 1._rk/iNP1n)
         end if

         ! Temperature response
         etP1 = self%q10P1X**((ETW-10._rk)/10._rk) - self%q10P1X**((ETW-32._rk)/3._rk)

         ! Production...........................................................

         ! calculate chl to C ratio.............................................

         !ChlCpp = max(min(phimP1X,chl1(I)/(p1c(I)-chl1(I)+zeroX)),ChlCmin)
         ChlCpp = chl1/p1c

         ! Gross photosynthetic activity :
         sumP1 = self%sumP1X*etP1*iNP1s*iNP1f
         phi = self%phiP1HX + (ChlCpp/self%phimP1X)*(self%phimP1X-self%phiP1HX)

         if (parEIR.gt.zeroX) THEN
            sumP1 = sumP1 * (1._rk-exp(-self%alphaP1X*parEIR*ChlCpp/sumP1)) * EXP(-self%betaP1X*parEIR*ChlCpp/sumP1)
            rho = (phi - ChlCmin) * (sumP1/(self%alphaP1X*parEIR*ChlCpp)) + ChlCmin
         else
            sumP1 = 0._rk
            rho = ChlCmin
         end if

#ifdef CENH
         ! Retrieve atmospheric pCO2 if using CO2-enhanced primary production.
         _GET_HORIZONTAL_(self%id_pco2a3,pco2a3)

         ! Enhancement factor (from MEECE D1.5) 379.48 = pco2a @ 2005
         cenh=1.+(pco2a3-379.48_rk)*0.0005_rk
#endif

         ! Nutrient-stress lysis rate :
         sdoP1 = (1._rk/(MIN( iNP1s, iNIP1 )+0.1_rk))*self%sdoP1X

         ! Excretion rate, as regulated by nutrient-stress
         seoP1 = sumP1*(1._rk-iNIP1)*(1._rk-self%pu_eaP1X)

         ! Activity-dependent excretion :
         seaP1 = sumP1*self%pu_eaP1X
         sugP1 = sumP1-seoP1-seaP1

         ! Apportioning of LOC- and DET- fraction of excretion/lysis fluxes:
         pe_R6P1 = MIN(self%qplP1cX/(qpP1c+ZeroX),self%qnlP1cX/(qnP1c+ZeroX))

         sP1R6 = pe_R6P1*sdoP1
         fP1R6c = sP1R6*P1cP

         if (self%docdyn) then
            ! Lysis produces labile DOM, excretion produces semi-labile DOM.
            fP1R1c = (1._rk-pe_R6P1)*sdoP1*P1cP
            fP1R2c = (seoP1 + seaP1)*P1c
            fP1RDc = fP1R1c+fP1R2c
         else
            ! Fixed ratio between production of labile and semi-labile DOM.
            fP1RDc = (1._rk-pe_R6P1)*sdoP1*P1cP + (seoP1 + seaP1)*P1c
            fP1R1c = fP1RDc*self%R1R2X
            fP1R2c = fP1RDc*(1._rk-self%R1R2X)
         end if

         ! Calcified matter:
         if (self%calcify) then
            _GET_(self%id_RainR,RainR)

            ! Calculate nutrient limitation impact on rain ratio:
            t=max(0._rk,ETW)  ! this is to avoid funny values of rain ratio when ETW ~ -2 degrees
            RainR = RainR * min((1._rk-iNP1p),iNP1n) * (t/(2._rk+t)) !* max(1.,P2c(I)/2.) removd as P2 is a broad class not just calicifiers
            RainR = max(RainR,0.005_rk)

            ! Compute virtual calcite attached to live cells. It is virtual in the sense that it has not been subtracted
            ! from the DIC budget yet. Thus, its carbon is not yet allocated. When consumed by predators, this is (partially)
            ! transformed into free liths. At that time, carbon used to form liths is subtracted from the DIC pool.
            _SET_DIAGNOSTIC_(self%id_lD,RainR*P1cP)

            ! For dying cells: convert virtual cell-attached calcite to actual calcite in free liths.
            ! Update DIC and lith balances accordingly (only now take away the DIC needed to form the calcite)
            _SET_ODE_(self%id_O3c,-fP1R6c*RainR/Cmass)
            _SET_ODE_(self%id_L2c, fP1R6c*RainR)

         end if

         ! Respiration..........................................................

         ! Rest respiration rate :
         srsP1 = etP1*self%srsP1X 

         ! Activity respiration rate :
         sraP1 = sugP1*self%pu_raP1X
#ifndef CENH
         ! Total respiration flux :
         fP1O3c = srsP1*P1cP+sraP1*P1c

         ! Gross production as flux from inorganic CO2 :
         fO3P1c = sumP1*P1c
#else
         fO3P1c = sumP1*P1c*cenh
         fP1O3c = srsP1*P1cP+sraP1*P1c*cenh
#endif
         rugP1 = fO3P1c

         ! Production and productivity
         sunP1 = sumP1-(seoP1+seaP1+sraP1)  ! net productivity
         runP1 = sunP1*P1c-srsP1*P1cP       ! net production

         ! To save net production
#ifndef CENH
         _SET_DIAGNOSTIC_(self%id_netP1,runP1)
#else
         _SET_DIAGNOSTIC_(self%id_netP1,(sumP1*cenh-seoP1-seaP1-sraP1*cenh)*P1c-srsP1*P1cP)
#endif

         ! Carbon Source equations

         rho = MIN(rho,0.1_rk)

         ! Chl changes (note that Chl is a component of PXc and not involved
         ! in mass balance)
         Chl_inc = rho*(sumP1-sraP1)*P1c
         Chl_loss = (sdoP1+srsP1)*Chl1P + (seoP1+seaP1)*Chl1

         _SET_ODE_(self%id_c,(fO3P1c-fP1O3c-fP1R6c-fP1RDc))

         _SET_ODE_(self%id_R1c,fP1R1c)
         _SET_ODE_(self%id_R2c,fP1R2c)
         _SET_ODE_(self%id_R6c,fP1R6c)
         _SET_ODE_(self%id_chl,(Chl_inc - Chl_loss))

         _SET_ODE_(self%id_O3c,(fP1O3c - fO3P1c)/CMass)
         _SET_ODE_(self%id_O2o,(fO3P1c*self%uB1c_O2X - fP1O3c*self%urB1_O2X))

         ! Phosphorus flux...........................................

         ! Lysis loss of phosphorus
         fP1R6p = sP1R6 * min(self%qplP1cX*P1cP,P1pP)
         fP1RDp = sdoP1 * P1pP - fP1R6p

         ! Net phosphorus uptake
         rumP1p = self%qurP1pX * N1pP * P1c
         misP1p = self%xqpP1X * qpRPIcX*P1cP - P1pP
         runP1p = sunP1*P1c * qpRPIcX*self%xqpP1X - srsP1*P1pP
         fN1P1p = MIN(rumP1p, runP1p+misP1p)

         ! Source equations
         _SET_ODE_(self%id_p,(fN1P1p-fP1RDp-fP1R6p))
         _SET_ODE_(self%id_N1p,-fN1P1p)
         _SET_ODE_(self%id_R6p,fP1R6p)
         _SET_ODE_(self%id_R1p,fP1RDp)

         ! Nitrogen flux.............................................

         ! Nitrogen loss by lysis
         fP1R6n = sP1R6 * min(self%qnlP1cX*P1cP,P1nP)
         fP1RDn = sdoP1 * P1nP - fP1R6n

         ! Net nitrogen uptake
         rumP1n3 = self%quP1n3X * N3nP * P1c
         rumP1n4 = self%quP1n4X * N4nP * P1c
         rumP1n = rumP1n3 + rumP1n4

         misP1n = self%xqnP1X * qnRPIcX*P1cP - P1nP
         runP1n = sunP1*P1c * qnRPIcX*self%xqnP1X - srsP1*P1nP
         fNIP1n = MIN(rumP1n, runP1n + misP1n)

         ! Partitioning over NH4 and NO3 uptake
         IF (fNIP1n .gt. 0._rk) THEN
            fN3P1n = fNIP1n * rumP1n3 / rumP1n
            fN4P1n = fNIP1n * rumP1n4 / rumP1n
         ELSE
            fN3P1n = 0._rk
            fN4P1n = fNIP1n
         ENDIF

         ! Source equations
         _SET_ODE_(self%id_n,(fN4P1n+fN3P1n-fP1RDn-fP1R6n))
         _SET_ODE_(self%id_N3n,-fN3P1n)
         _SET_ODE_(self%id_N4n,-fN4P1n)
         _SET_ODE_(self%id_R6n,fP1R6n)
         _SET_ODE_(self%id_R1n,fP1RDn)

         if (self%use_Si) then
            ! Silicate flux.............................................
            _GET_(self%id_s,P1sP)

            ! Excretion loss of silicate
            fP1R6s = sdoP1 * P1sP

            ! Loss of excess silicate (qsP1c > qsP1cX)
            fP1N5s = MAX ( 0._rk, P1sP-self%qsP1cX * P1cP)

            ! Net silicate uptake
            fN5P1s = MAX ( 0._rk, self%qsP1cX*runP1) - fP1N5s
#ifdef SAVEFLX
            fN5PXs(I) = fN5P1s
#endif

            ! Source equations
            _SET_ODE_(self%id_s,(fN5P1s - fP1R6s))
            _SET_ODE_(self%id_N5s,-fN5P1s)
            _SET_ODE_(self%id_R6s, fP1R6s)
         end if

#ifdef IRON
         ! Iron flux................................................

         ! Obtain internal iron concentration (umol m-3, excludes background), and external biologically available iron.
         _GET_(self%id_f,P1fP)
         _GET_(self%id_N7f,N7fP)

         ! Iron loss by lysis
         !  Because its high affinity with particles all the iron lost from phytoplankton by lysis is supposed to be 
         !  associated to organic particulate detritus. (luca)

         fP1R6f = sdoP1 * P1fP

         ! Net iron uptake
         rumP1f = self%qurP1fX * N7fP * P1c
         misP1f = self%qfRP1cX*P1cP - P1fP
         runP1f = sunP1*P1c * self%qfRP1cX - srsP1*P1fP
         fN7P1f = MIN(rumP1f, runP1f+misP1f)

         ! Source equations
         _SET_ODE_(self%id_f,(fN7P1f-fP1R6f))
         _SET_ODE_(self%id_N7f,-fN7P1f)
         _SET_ODE_(self%id_R6f,fP1R6f)
#endif

      ! Leave spatial loops (if any)
      _LOOP_END_

   end subroutine do

   function get_sinking_rate(self,_ARGUMENTS_LOCAL_) result(SDP1)
      class (type_pml_ersem_primary_producer),intent(in) :: self
      _DECLARE_ARGUMENTS_LOCAL_

      real(rk) :: P1c,P1p,P1n
      real(rk) :: qpP1c,qnP1c
      real(rk) :: SDP1
      real(rk) :: iNP1p,iNP1n,iNIP1

      _GET_WITH_BACKGROUND_(self%id_c,P1c)
      _GET_WITH_BACKGROUND_(self%id_p,P1p)
      _GET_WITH_BACKGROUND_(self%id_n,P1n)
      
      qpP1c = P1p/P1c
      qnP1c = P1n/P1c

      iNP1p = MIN(1._rk,  &
               MAX(0._rk, (qpP1c-self%qplP1cX) / (self%xqcP1pX*qpRPIcX-self%qplP1cX) ))
      iNP1n = MIN(1._rk,  &
               MAX(0._rk, (qnP1c-self%qnlP1cX) / (self%xqcP1nX*qnRPIcX-self%qnlP1cX) ))

      if (self%LimnutX==0) then  ! NB select case would be cleaner but makes vectorization impossible for ifort 14
         iNIP1 = (iNP1p * iNP1n)**0.5_rk
      elseif (self%LimnutX==1) then
         iNIP1 = MIN(iNP1p, iNP1n)
      else
         iNIP1 = 2.0_rk / (1._rk/iNP1p + 1._rk/iNP1n)
      end if

      ! Sedimentation and resting stages.....................................
      SDP1 = self%resP1mX * MAX(0._rk, (self%esNIP1X - iNIP1)) + self%rP1mX
   end function

   subroutine get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)
      class (type_pml_ersem_primary_producer),intent(in) :: self
      _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

      real(rk) :: SDP1

      _LOOP_BEGIN_

         ! Retrieve local state
         SDP1 = -get_sinking_rate(self,_ARGUMENTS_LOCAL_)
         _SET_VERTICAL_MOVEMENT_(self%id_c,SDP1)
         _SET_VERTICAL_MOVEMENT_(self%id_p,SDP1)
         _SET_VERTICAL_MOVEMENT_(self%id_n,SDP1)
         if (self%use_Si) _SET_VERTICAL_MOVEMENT_(self%id_s,SDP1)
         _SET_VERTICAL_MOVEMENT_(self%id_chl,SDP1)
#ifdef IRON
         _SET_VERTICAL_MOVEMENT_(self%id_f,SDP1)
#endif

      _LOOP_END_

   end subroutine get_vertical_movement

end module