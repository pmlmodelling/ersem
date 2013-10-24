#include "fabm_driver.h"

#define IRON
!#define CENH

module pml_ersem_vphyt

   use fabm_types

   implicit none

   private

   type,extends(type_base_model),public :: type_pml_ersem_vphyt
      ! Variable identifiers
      type (type_state_variable_id)      :: id_P1c,id_P1n,id_P1p,id_P1s,id_Chl1
      type (type_state_variable_id)      :: id_O3c,id_O2o
      type (type_state_variable_id)      :: id_N5s,id_N1p,id_N3n,id_N4n
      type (type_state_variable_id)      :: id_R1c,id_R1p,id_R1n,id_R2c,id_R6c,id_R6p,id_R6n,id_R6s
#ifdef IRON   
      type (type_state_variable_id)      :: id_P1f
      type (type_state_variable_id)      :: id_N7f,id_R6f
#endif
#ifdef CENH
      type (type_horizontal_dependency_id) :: id_pco2a3
#endif
      type (type_dependency_id)          :: id_EIR,id_ETW
      type (type_diagnostic_variable_id) :: id_netP1

      ! Parameters
      real(rk) :: qnrpicX,qprpicX,sump1X
      real(rk) :: q10p1X,srsp1X,pu_eap1X,pu_rap1X,chp1sX,qnlp1cX,qplp1cX,xqcp1pX
      real(rk) :: xqcp1nX,xqnp1X,xqpp1X,qup1n3X,qup1n4X,qurp1pX,qsp1cX,esnip1X
      real(rk) :: resp1mX,sdop1X
      real(rk) :: alphaP1X,betaP1X,phimP1X,phiP1HX
      real(rk) :: pEIR_eowX,R1R2X,uB1c_O2X,urB1_O2X
#ifdef IRON
      real(rk) :: qflP1cX,qfRP1cX,qurP1fX
#endif
      integer :: LimnutX
      logical :: use_Si
   contains
!     Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: get_vertical_movement
   end type type_pml_ersem_vphyt

   real(rk),parameter :: CMass = 12._rk
   real(rk),parameter :: ZeroX = 1e-8_rk
   real(rk),parameter :: secs_pr_day = 86400.0_rk
   real(rk),parameter :: ChlCmin=0.0067

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_pml_ersem_vphyt),intent(inout),target :: self
      integer,                     intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%use_Si,   'use_Si',  default=.true.)
      call self%get_parameter(self%qnrpicX,  'qnrpicX')
      call self%get_parameter(self%qprpicX,  'qprpicX') 
      call self%get_parameter(self%sump1X,   'sump1X')
      call self%get_parameter(self%q10p1X,   'q10p1X')
      call self%get_parameter(self%srsp1X,   'srsp1X')
      call self%get_parameter(self%pu_eap1X, 'pu_eap1X')
      call self%get_parameter(self%pu_rap1X, 'pu_rap1X')
      call self%get_parameter(self%qnlp1cX,  'qnlp1cX')
      call self%get_parameter(self%qplp1cX,  'qplp1cX')
      call self%get_parameter(self%xqcp1pX,  'xqcp1pX')
      call self%get_parameter(self%xqcp1nX,  'xqcp1nX')
      call self%get_parameter(self%xqnp1X,   'xqnp1X')
      call self%get_parameter(self%xqpp1X,   'xqpp1X')
      call self%get_parameter(self%qup1n3X,  'qup1n3X')
      call self%get_parameter(self%qup1n4X,  'qup1n4X')
      call self%get_parameter(self%qurp1pX,  'qurp1pX')
      if (self%use_Si) then
         call self%get_parameter(self%qsp1cX,'qsp1cX')
         call self%get_parameter(self%chp1sX,'chp1sX')
      end if
      call self%get_parameter(self%esnip1X,  'esnip1X')
      call self%get_parameter(self%resp1mX,  'resp1mX')
      call self%get_parameter(self%sdop1X,   'sdop1X')
      call self%get_parameter(self%alphaP1X, 'alphaP1X')
      call self%get_parameter(self%betaP1X,  'betaP1X')
      call self%get_parameter(self%phimP1X,  'phimP1X')
      call self%get_parameter(self%phiP1HX,  'phiP1HX')
      call self%get_parameter(self%pEIR_eowX,'pEIR_eowX')
      call self%get_parameter(self%LimnutX,  'LimnutX')
      call self%get_parameter(self%R1R2X,    'R1R2X')
      call self%get_parameter(self%uB1c_O2X, 'uB1c_O2X')
      call self%get_parameter(self%urB1_O2X, 'urB1_O2X')
#ifdef IRON
      call self%get_parameter(self%qflP1cX,  'qflP1cX')
      call self%get_parameter(self%qfRP1cX,  'qfRP1cX')
      call self%get_parameter(self%qurP1fX,  'qurP1fX')
#endif

      ! Register state variables
      call self%register_state_variable(self%id_P1c, 'c',  'mg C/m^3',   'C',             1.e-4_rk,   minimum=0.0_rk)
      call self%register_state_variable(self%id_P1n, 'n',  'mmol N/m^3', 'N',             1.26e-6_rk, minimum=0.0_rk)
      call self%register_state_variable(self%id_P1p, 'p',  'mmol P/m^3', 'P',             4.288e-8_rk,minimum=0.0_rk)
      if (self%use_Si) call self%register_state_variable(self%id_P1s, 's',  'mmol Si/m^3','Si',            1.e-6_rk,   minimum=0.0_rk)
      call self%register_state_variable(self%id_Chl1,'Chl','mg C/m^3',   'Chlorophyll-a', 3.e-6_rk,   minimum=0.0_rk)
#ifdef IRON   
      call self%register_state_variable(self%id_P1f, 'f',  'umol Fe/m^3','Fe',            5.e-6_rk,   minimum=0.0_rk)
#endif

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

      ! Register links to external semi-labile dissolved organic matter pools.
      call self%register_state_dependency(self%id_R2c,'R2c','mg C/m^3','Semi-labile DOC')    

      ! Register links to external particulate organic matter (medium-size) pools.
      call self%register_state_dependency(self%id_R6c,'RPc','mg C/m^3',   'POC')    
      call self%register_state_dependency(self%id_R6p,'RPp','mmol P/m^3', 'POP')    
      call self%register_state_dependency(self%id_R6n,'RPn','mmol N/m^3', 'PON')    
      if (self%use_Si) call self%register_state_dependency(self%id_R6s,'RPs','mmol Si/m^3','POSi')    
#ifdef IRON   
      call self%register_state_dependency(self%id_R6f,'RPf','umol Fe/m^3','POFe')    
#endif

      ! Register links to external total dissolved inorganic carbon, dissolved oxygen pools
      call self%register_state_dependency(self%id_O3c,'O3c','mmol C/m^3','Carbon Dioxide')    
      call self%register_state_dependency(self%id_O2o,'O2o','mmol O/m^3','Oxygen')    
#ifdef CENH   
      call self%register_dependency(self%id_pco2a3,standard_variables%mole_fraction_of_carbon_dioxide_in_air)    
#endif

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_netP1,'netP1','mg C/m^3/d','net primary production',time_treatment=time_treatment_averaged)

      ! Register environmental dependencies (temperature, shortwave radiation)
      call self%register_dependency(self%id_EIR,standard_variables%downwelling_shortwave_flux)
      call self%register_dependency(self%id_ETW,standard_variables%temperature)

   end subroutine
   
   subroutine do(self,_ARGUMENTS_DO_)

      class (type_pml_ersem_vphyt),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

   ! !LOCAL VARIABLES:
      real(rk) :: ETW,EIR
      real(rk) :: P1c, P1p, P1n, P1s, Chl1
      real(rk) :: P1cP,P1pP,P1nP,P1sP,Chl1P
      real(rk) :: N5s,N1pP,N3nP,N4nP
      real(rk) :: iNP1n,iNP1p,iNP1s,iNP1f,iNIP1
      real(rk) :: qpP1c,qnP1c,qsP1c
      real(rk) :: parEIR
    
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
#ifdef DOCDYN
      real(rk) :: fP1R1c,fP1R2c
#endif
#ifdef CENH
      real(rk) :: pco2a3,cenh
#endif

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! First retrieve local state.
         ! NB: _GET_SAFE_ retrieves the value of the state variable, minus an offset
         ! (if any) specified by the owning model during registration of the state variable.
         ! If the returned value would be negative, 0.0 is returned instead.

         _GET_(self%id_P1c,P1c)
         _GET_(self%id_P1p,P1p)
         _GET_(self%id_P1n,P1n)
#ifdef IRON      
         _GET_(self%id_P1f,P1f)
#endif
         _GET_(self%id_Chl1,chl1)

         if (self%use_Si) then
            _GET_(self%id_P1s,P1s)
            _GET_(self%id_N5s,N5s)
            _GET_SAFE_(self%id_P1s,P1sP)
         end if

         _GET_SAFE_(self%id_P1c,P1cP)
         _GET_SAFE_(self%id_P1p,P1pP)
         _GET_SAFE_(self%id_P1n,P1nP)
         _GET_SAFE_(self%id_Chl1,chl1P)
         _GET_SAFE_(self%id_N1p,N1pP)
         _GET_SAFE_(self%id_N3n,N3nP)
         _GET_SAFE_(self%id_N4n,N4nP)
#ifdef IRON      
         _GET_SAFE_(self%id_P1f,P1fP)
         _GET_SAFE_(self%id_N7f,N7fP)
#endif

         ! Get environmental dependencies (water temperature, shortwave radation)
         _GET_(self%id_ETW,ETW)
         _GET_(self%id_EIR,EIR)

#ifdef CENH      
         ! Get atmospheric pCO2
         _GET_HORIZONTAL_(self%id_pco2a3,pco2a3)
#endif

         qpP1c = P1p/P1c
         qnP1c = P1n/P1c
         if (self%use_Si) qsP1c = P1s/P1c
#ifdef IRON
         qfP1c = P1f/P1c
#endif

         iNP1p = MIN(1._rk,  &
                 MAX(0._rk, (qpP1c-self%qplP1cX) / (self%xqcP1pX*self%qpRPIcX-self%qplP1cX) ))
         iNP1n = MIN(1._rk,  &
                 MAX(0._rk, (qnP1c-self%qnlP1cX) / (self%xqcP1nX*self%qnRPIcX-self%qnlP1cX) ))
#ifdef IRON
         iNP1f = MIN(1._rk,  &
                 MAX(ZeroX, (qfP1c-self%qflP1cX) / (self%qfRP1cX-self%qflP1cX) ))
#endif
         if (self%use_Si) then
            iNP1s = MIN(1._rk, N5s/(N5s+self%chP1sX))
         else
            iNP1s = 1.0_rk
         end if

         select case (self%LimnutX)
            case (0)
               iNIP1 = (iNP1p * iNP1n)**0.5_rk
            case (1)
               iNIP1 = MIN(iNP1p, iNP1n)
            case (2)
               iNIP1 = 2.0_rk / (1._rk/iNP1p + 1._rk/iNP1n)
         end select
            
   !..Regulation factors...................................................

   !..Temperature response
         etP1 = self%q10P1X**((ETW-10._rk)/10._rk) - self%q10P1X**((ETW-32._rk)/3._rk)

   !..Production...........................................................

   !..calculate chl to C ratio.............................................

         !ChlCpp = max(min(phimP1X,chl1(I)/(p1c(I)-chl1(I)+zeroX)),ChlCmin)
         ChlCpp = chl1/p1c

   !..Gross photosynthetic activity :
#ifdef IRON
         sumP1 = self%sumP1X*etP1*iNP1s*iNP1f
#else
         sumP1 = self%sumP1X*etP1*iNP1s
#endif          
         phi = self%phiP1HX + (ChlCpp/self%phimP1X)*(self%phimP1X-self%phiP1HX)

         parEIR = self%pEIR_eowX*EIR
         if (parEIR.gt.zeroX) THEN
            sumP1 = sumP1 * (1._rk-exp(-self%alphaP1X*parEIR*ChlCpp/sumP1)) * EXP(-self%betaP1X*parEIR*ChlCpp/sumP1)
            rho = (phi - ChlCmin) * (sumP1/(self%alphaP1X*parEIR*ChlCpp)) + ChlCmin
         else
            sumP1 = 0._rk
            rho = ChlCmin
         end if
#ifdef CENH
   ! Enhancement factor (from MEECE D1.5) 379.48 = pco2a @ 2005
         cenh=1.+(pco2a3-379.48_rk)*0.0005_rk
#endif

   !..Nutrient-stress lysis rate :
         sdoP1 = (1._rk/(MIN( iNP1s, iNIP1 )+0.1_rk))*self%sdoP1X

   !..Excretion rate, as regulated by nutrient-stress
         seoP1 = sumP1*(1._rk-iNIP1)*(1._rk-self%pu_eaP1X)

   !..activity-dependent excretion :
         seaP1 = sumP1*self%pu_eaP1X
         sugP1 = sumP1-seoP1-seaP1

   !..Apportioning of LOC- and DET- fraction of excretion/lysis fluxes:
         pe_R6P1 = MIN(self%qplP1cX/(qpP1c+ZeroX),self%qnlP1cX/(qnP1c+ZeroX))

         sP1R6 = pe_R6P1*sdoP1   !Jorn: R4 for flagellates
         fP1R6c = sP1R6*P1cP     !Jorn: R4 for flagellates

#ifdef DOCDYN
         fP1R1c=  (1._rk-pe_R6P1)*sdoP1*P1cP
         fP1R2c=  (seoP1 + seaP1)*P1c
         fP1RDc = fP1R1c+fP1R2c
#else
         fP1RDc = (1._rk-pe_R6P1)*sdoP1*P1cP + (seoP1 + seaP1)*P1c
#endif      

#ifdef CALC
! Calcified matter:
      fO3L2c(I) = fO3L2c(I) + fP2R4c(I)*RainR(I)      !Jorn: flagellates only!
#endif
   !..Respiration..........................................................

   !..Rest respiration rate :
         srsP1 = etP1*self%srsP1X 

   !..Activity respiration rate :
         sraP1 = sugP1*self%pu_raP1X
#ifndef CENH
!..Total respiration flux :
         fP1O3c = srsP1*P1cP+sraP1*P1c
      
!..Gross production as flux from inorganic CO2 to P1 :
         fO3P1c = sumP1*P1c
#else
         fO3P1c = fO3P1c*cenh !Jorn: a bug? where is fO3P1c set?
         fP1O3c = srsP1*P1cP+sraP1*P1c*cenh
#endif
         rugP1 = fO3P1c

   !..Production and productivity
         sunP1 = sumP1-(seoP1+seaP1+sraP1)  ! net productivity
         runP1 = sunP1*P1c-srsP1*P1cP       ! net production

   !..To save net production
#ifndef CENH
         _SET_DIAGNOSTIC_(self%id_netP1,runP1)
#else
         _SET_DIAGNOSTIC_(self%id_netP1,(sumP1*cenh-seoP1-seaP1-sraP1*cenh)*P1c-srsP1*P1cP)
#endif

   !..Carbon Source equations
        
         rho = MIN(rho,0.1_rk)
         ! Chl changes (note that Chl is a component of PXc and not involved
         ! in mass balance)
         Chl_inc = rho*(sumP1-sraP1)*P1c
         Chl_loss = (sdoP1+srsP1)*Chl1P + (seoP1+seaP1)*Chl1

         _SET_ODE_(self%id_P1c,(fO3P1c-fP1O3c-fP1R6c-fP1RDc)/secs_pr_day) ! Jorn: R4 in flagellates
#ifdef DOCDYN
         _SET_ODE_(self%id_R1c,fP1R1c/secs_pr_day)
         _SET_ODE_(self%id_R2c,fP1R2c/secs_pr_day)
#else
         _SET_ODE_(self%id_R1c,(fP1RDc*self%R1R2X)/secs_pr_day)
         _SET_ODE_(self%id_R2c,fP1RDc*(1._rk-self%R1R2X)/secs_pr_day)
#endif
         _SET_ODE_(self%id_R6c,fP1R6c/secs_pr_day) ! Jorn: R4 in flagellates
         _SET_ODE_(self%id_Chl1,(Chl_inc - Chl_loss)/secs_pr_day)

         _SET_ODE_(self%id_O3c,(fP1O3c - fO3P1c)/CMass/secs_pr_day)
         _SET_ODE_(self%id_O2o,(fO3P1c*self%uB1c_O2X - fP1O3c*self%urB1_O2X)/secs_pr_day)
         
   !..Phosphorus flux through P1...........................................

   !..Lysis loss of phosphorus
         fP1R6p = sP1R6 * min(self%qplP1cX*P1cP,P1pP) ! Jorn: R4 for flagellates
         fP1RDp = sdoP1 * P1pP - fP1R6p               ! Jorn: R4 for flagellates

   !..Net phosphorus uptake
         rumP1p = self%qurP1pX * N1pP * P1c
         misP1p = self%xqpP1X * self%qpRPIcX*P1cP - P1pP
         runP1p = sunP1*P1c * self%qpRPIcX*self%xqpP1X - srsP1*P1pP
         fN1P1p = MIN(rumP1p, runP1p+misP1p)

   !..Source equations
         _SET_ODE_(self%id_P1p,(fN1P1p-fP1RDp-fP1R6p)/secs_pr_day) ! Jorn: R4 for flagellates
         _SET_ODE_(self%id_N1p,-fN1P1p/secs_pr_day)
         _SET_ODE_(self%id_R6p,fP1R6p/secs_pr_day) ! Jorn: R4 for flagellates
         _SET_ODE_(self%id_R1p,fP1RDp/secs_pr_day)

   !..Nitrogen flux through P1.............................................

   !..Nitrogen loss by lysis
         fP1R6n = sP1R6 * min(self%qnlP1cX*P1cP,P1nP)
         fP1RDn = sdoP1 * P1nP - fP1R6n ! Jorn: R4 for flagellates

   !..Net nitrogen uptake
         rumP1n3 = self%quP1n3X * N3nP * P1c
         rumP1n4 = self%quP1n4X * N4nP * P1c
         rumP1n = rumP1n3 + rumP1n4

         misP1n = self%xqnP1X * self%qnRPIcX*P1cP - P1nP
         runP1n = sunP1*P1c * self%qnRPIcX*self%xqnP1X - srsP1*P1nP
         fNIP1n = MIN(rumP1n, runP1n + misP1n)

   !..Partitioning over NH4 and NO3 uptake
         IF (fNIP1n .gt. 0._rk) THEN
            fN3P1n = fNIP1n * rumP1n3 / rumP1n
            fN4P1n = fNIP1n * rumP1n4 / rumP1n
         ELSE
            fN3P1n = 0._rk
            fN4P1n = fNIP1n
         ENDIF

   !..Source equations
         _SET_ODE_(self%id_P1n,(fN4P1n+fN3P1n-fP1RDn-fP1R6n)/secs_pr_day)
         _SET_ODE_(self%id_N3n,-fN3P1n/secs_pr_day)
         _SET_ODE_(self%id_N4n,-fN4P1n/secs_pr_day)
         _SET_ODE_(self%id_R6n,fP1R6n/secs_pr_day) ! Jorn: R4 for flagellates
         _SET_ODE_(self%id_R1n,fP1RDn/secs_pr_day)

         if (self%use_Si) then
   !..Silicate flux through P1.............................................

   !..Excretion loss of silicate
            fP1R6s = sdoP1 * P1sP

   !..Loss of excess silicate (qsP1c > qsP1cX)
            fP1N5s = MAX ( 0._rk, P1sP-self%qsP1cX * P1cP)

   !..Net silicate uptake
            fN5P1s = MAX ( 0._rk, self%qsP1cX*runP1) - fP1N5s
#ifdef SAVEFLX
            fN5PXs(I) = fN5P1s
#endif

   !..Source equations
            _SET_ODE_(self%id_P1s,(fN5P1s - fP1R6s)/secs_pr_day)
            _SET_ODE_(self%id_N5s,-fN5P1s/secs_pr_day)
            _SET_ODE_(self%id_R6s, fP1R6s/secs_pr_day)
         end if

#ifdef IRON
   !..Iron flux through P1................................................

   !..Iron loss by lysis
   !  Because its high affinity with particles all the iron lost from phytoplankton by lysis is supposed to be 
   !  associated to organic particulate detritus. (luca)

         fP1R6f = sdoP1 * P1fP   ! Jorn: R4 for flagellates

   !..net iron uptake
         rumP1f = self%qurP1fX * N7fP * P1c
         misP1f = self%qfRP1cX*P1cP - P1fP
         runP1f = sunP1*P1c * self%qfRP1cX - srsP1*P1fP
         fN7P1f = MIN(rumP1f, runP1f+misP1f)

   !..Source equations
         _SET_ODE_(self%id_P1f,(fN7P1f-fP1R6f)/secs_pr_day) ! Jorn: R4 for flagellates
         _SET_ODE_(self%id_N7f,-fN7P1f/secs_pr_day)
         _SET_ODE_(self%id_R6f,fP1R6f/secs_pr_day) ! Jorn: R4 for flagellates
#endif

   !..add dia 
   !      IF ((P1c(I) .LT. 0.2).AND. (parEIR .GT. 0.5)) THEN 
   !       SP1c(I) = SP1c(I) + (0.2 - P1c(I)) 
   !       SP1n(I) = SP1n(I) + (0.2 - P1c(I))*16.0/106.0 
   !       SP1p(I) = SP1p(I) + (0.2 - P1c(I))/106.0 
   !       SP1s(I) = SP1s(I) + (0.2 - P1c(I))*16.0/106.0 
   !       SChl1(I) = SChl1(I) + (0.2 - P1c(I))*ChlCpp 
   !      END IF 

      ! Leave spatial loops (if any)
      _LOOP_END_

   end subroutine do

   subroutine get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)
      class (type_pml_ersem_vphyt),intent(in) :: self
      _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

      real(rk) :: SDP1
      real(rk) :: P1c,P1p,P1n,qpP1c,qnP1c,iNP1p,iNP1n,iNIP1

      _LOOP_BEGIN_

         ! Retrieve local state
         _GET_(self%id_P1c,P1c)
         _GET_(self%id_P1p,P1p)
         _GET_(self%id_P1n,P1n)

         qpP1c = P1p/P1c
         qnP1c = P1n/P1c

         iNP1p = MIN(1._rk,  &
                 MAX(0._rk, (qpP1c-self%qplP1cX) / (self%xqcP1pX*self%qpRPIcX-self%qplP1cX) ))
         iNP1n = MIN(1._rk,  &
                 MAX(0._rk, (qnP1c-self%qnlP1cX) / (self%xqcP1nX*self%qnRPIcX-self%qnlP1cX) ))

         select case (self%LimnutX)
            case (0)
               iNIP1 = (iNP1p * iNP1n)**.5_rk
            case (1)
               iNIP1 = MIN(iNP1p, iNP1n)
            case (2)
               iNIP1 = 2.0_rk / (1._rk/iNP1p + 1._rk/iNP1n)
         end select

!..sedimentation and resting stages.....................................
         SDP1 = -self%resP1mX * MAX(0._rk, (self%esNIP1X - iNIP1))/secs_pr_day
         _SET_VERTICAL_MOVEMENT_(self%id_P1c,SDP1)
         _SET_VERTICAL_MOVEMENT_(self%id_P1p,SDP1)
         _SET_VERTICAL_MOVEMENT_(self%id_P1n,SDP1)
         if (self%use_Si) _SET_VERTICAL_MOVEMENT_(self%id_P1s,SDP1)
         _SET_VERTICAL_MOVEMENT_(self%id_Chl1,SDP1)
#ifdef IRON
         _SET_VERTICAL_MOVEMENT_(self%id_P1f,SDP1)
#endif

      _LOOP_END_

   end subroutine get_vertical_movement

end module