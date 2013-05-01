#include "fabm_driver.h"
#define IRON
module pml_ersem_vphyt1

   use fabm_types
   use fabm_driver

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public pml_ersem_vphyt1_create, type_pml_ersem_vphyt1
!
! !PRIVATE DATA MEMBERS:
   real(rk), parameter :: secs_pr_day = 86400.0_rk
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman

   type,extends(type_base_model) :: type_pml_ersem_vphyt1
!     Variable identifiers
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

      real(rk) :: qnrpicX,qprpicX,sump1X
      real(rk) :: q10p1X,srsp1X,pu_eap1X,pu_rap1X,chp1sX,qnlp1cX,qplp1cX,xqcp1pX
      real(rk) :: xqcp1nX,xqnp1X,xqpp1X,qup1n3X,qup1n4X,qurp1pX,qsp1cX,esnip1X
      real(rk) :: resp1mX,sdop1X
      real(rk) :: alphaP1X,betaP1X,phimP1X,phiP1HX
      real(rk) :: pEIR_eowX,ChlCmin,R1R2X,uB1c_O2X,urB1_O2X
#ifdef IRON
      real(rk) :: qflP1cX,qfRP1cX,qurP1fX
#endif
      integer :: LimnutX

      contains

!     Model procedures
      procedure :: do
      procedure :: get_vertical_movement

   end type type_pml_ersem_vphyt1

   real(rk),parameter :: CMass = 12._rk
      
contains
      
   function pml_ersem_vphyt1_create(configunit,name,parent) result(self)
!
! !DESCRIPTION:
!  Here, the npzd namelist is read and the variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
      integer,                          intent(in)    :: configunit
      character(len=*),                 intent(in)    :: name
      class (type_base_model),target,   intent(inout) :: parent
      class (type_pml_ersem_vphyt1),    pointer       :: self
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
      real(rk) :: qnrpicX,qprpicX,sump1X
      real(rk) :: q10p1X,srsp1X,pu_eap1X,pu_rap1X,chp1sX,qnlp1cX,qplp1cX,xqcp1pX
      real(rk) :: xqcp1nX,xqnp1X,xqpp1X,qup1n3X,qup1n4X,qurp1pX,qsp1cX,esnip1X
      real(rk) :: resp1mX,sdop1X
      real(rk) :: alphaP1X,betaP1X,phimP1X,phiP1HX
      real(rk) :: pEIR_eowX,ChlCmin,R1R2X,uB1c_O2X,urB1_O2X
      real(rk) :: qflP1cX,qfRP1cX,qurP1fX
      integer :: LimnutX
       
      character(len=64) :: O3c_variable,O2o_variable, &
                           N1p_variable,N3n_variable,N4n_variable,N5s_variable,N7f_variable, &
                           R1c_variable,R1p_variable,R1n_variable,R2c_variable,R6c_variable,R6p_variable,R6n_variable,R6s_variable,R6f_variable

      namelist /pml_ersem_vphyt1/ qnrpicX,qprpicX,sump1X, &
        q10p1X,srsp1X,pu_eap1X,pu_rap1X,chp1sX,qnlp1cX,qplp1cX,xqcp1pX, &
        xqcp1nX,xqnp1X,xqpp1X,qup1n3X,qup1n4X,qurp1pX,qsp1cX,esnip1X, &
        resp1mX,sdop1X,alphaP1X,betaP1X,phimP1X,phiP1HX, &
        pEIR_eowX,ChlCmin,LimnutX,R1R2X,uB1c_O2X,urB1_O2X, &
        qflP1cX,qfRP1cX,qurP1fX, &
        O3c_variable,N5s_variable,R1c_variable,R2c_variable, &
       R1p_variable,R6p_variable,R1n_variable,R6n_variable,R6s_variable,R6f_variable
!EOP
!-----------------------------------------------------------------------
!BOC
      allocate(self)
      call self%initialize(name,parent)

      ! Initialize namelist parameters
      pEIR_eowX = 0.5_rk
      ChlCmin   = 0.0067_rk
      uB1c_O2X  = 0.11_rk
      urB1_O2X  = 0.1_rk
      LimnutX   = 1
      O3c_variable = 'pml_ersem_gas_dynamics_O3c'
      O2o_variable = 'pml_ersem_gas_dynamics_O2o'
      N5s_variable = 'pml_ersem_nutrients_N5s'
      N1p_variable = 'pml_ersem_nutrients_N1p'
      N3n_variable = 'pml_ersem_nutrients_N3n'
      N4n_variable = 'pml_ersem_nutrients_N4n'
      N7f_variable = 'pml_ersem_nutrients_N7f'
      R1c_variable = 'pml_ersem_dom_R1c'
      R1p_variable = 'pml_ersem_dom_R1p'
      R1n_variable = 'pml_ersem_dom_R1n'
      R2c_variable = 'pml_ersem_dom_R2c'
      R6c_variable = 'pml_ersem_pom_R6c'
      R6p_variable = 'pml_ersem_pom_R6p'
      R6n_variable = 'pml_ersem_pom_R6n'
      R6s_variable = 'pml_ersem_pom_R6s'
      R6f_variable = 'pml_ersem_pom_R6f'

      ! Read the namelist
      read(configunit,nml=pml_ersem_vphyt1)

      ! Store parameter values in our own derived type
      ! NB: all rates must be provided in values per day,
      ! and are converted here to values per second.
       self%qnrpicX = qnrpicX
       self%qprpicX = qprpicX
       self%sump1X = sump1X
       self%q10p1X = q10p1X
       self%srsp1X = srsp1X
       self%pu_eap1X = pu_eap1X
       self%pu_rap1X = pu_rap1X
       self%chp1sX = chp1sX
       self%qnlp1cX = qnlp1cX
       self%qplp1cX = qplp1cX
       self%xqcp1pX = xqcp1pX
       self%xqcp1nX = xqcp1nX
       self%xqnp1X = xqnp1X
       self%xqpp1X = xqpp1X
       self%qup1n3X = qup1n3X
       self%qup1n4X = qup1n4X
       self%qurp1pX = qurp1pX
       self%qsp1cX = qsp1cX
       self%esnip1X = esnip1X
       self%resp1mX = resp1mX
       self%sdop1X = sdop1X
       self%alphaP1X = alphaP1X
       self%betaP1X = betaP1X
       self%phimP1X = phimP1X
       self%phiP1HX = phiP1HX
       self%pEIR_eowX = pEIR_eowX
       self%ChlCmin = ChlCmin
       self%LimnutX = LimnutX
       self%R1R2X = R1R2X
       self%uB1c_O2X = uB1c_O2X
       self%urB1_O2X = urB1_O2X
#ifdef IRON
       self%qflP1cX = qflP1cX
       self%qfRP1cX = qfRP1cX
       self%qurP1fX = qurP1fX
#endif

      ! Register state variables
      call self%register_state_variable(self%id_P1c, 'P1c', 'mg C/m^3',  'Diatoms C', 1.e-4_rk,    minimum=_ZERO_)
      call self%register_state_variable(self%id_P1n, 'P1n', 'mmol N/m^3','Diatoms N', 1.26e-6_rk,  minimum=_ZERO_)
      call self%register_state_variable(self%id_P1p, 'P1p', 'mmol P/m^3','Diatoms P', 4.288e-8_rk, minimum=_ZERO_)
      call self%register_state_variable(self%id_P1s, 'P1s', 'mmol S/m^3','Diatoms S', 1.e-6_rk,    minimum=_ZERO_)
      call self%register_state_variable(self%id_Chl1,'Chl1','mg C/m^3',  'Diatoms Chlorophyll-a', 3.e-6_rk, minimum=_ZERO_)
#ifdef IRON   
      call self%register_state_variable(self%id_P1f, 'P1f', 'umol F/m^3','Diatoms F', 5.e-6_rk, minimum=_ZERO_)
#endif

      ! Register links to external nutrient pools.
      call self%register_state_dependency(self%id_N1p,N1p_variable)    
      call self%register_state_dependency(self%id_N3n,N3n_variable)    
      call self%register_state_dependency(self%id_N4n,N4n_variable)    
      call self%register_state_dependency(self%id_N5s,N5s_variable)    
#ifdef IRON   
      call self%register_state_dependency(self%id_N7f,N7f_variable)    
#endif

      ! Register links to external labile dissolved organic matter pools.
      call self%register_state_dependency(self%id_R1c,R1c_variable)
      call self%register_state_dependency(self%id_R1p,R1p_variable)    
      call self%register_state_dependency(self%id_R1n,R1n_variable)    

      ! Register links to external semi-labile dissolved organic matter pools.
      call self%register_state_dependency(self%id_R2c,R2c_variable)    

      ! Register links to external particulate organic matter (medium-size) pools.
      call self%register_state_dependency(self%id_R6c,R6c_variable)    
      call self%register_state_dependency(self%id_R6p,R6p_variable)    
      call self%register_state_dependency(self%id_R6n,R6n_variable)    
      call self%register_state_dependency(self%id_R6s,R6s_variable)    
#ifdef IRON   
      call self%register_state_dependency(self%id_R6f,R6f_variable)    
#endif

      ! Register links to external total dissolved inorganic carbon, dissolved oxygen pools
      call self%register_state_dependency(self%id_O3c,O3c_variable)    
      call self%register_state_dependency(self%id_O2o,O2o_variable)    
#ifdef CENH   
      call self%register_dependency(self%id_pco2a3,'mole_fraction_of_carbon_dioxide_in_air')    
#endif

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_netP1,'netP1','mg C/m^3/d','net primary production',time_treatment=time_treatment_step_averaged)

      ! Register environmental dependencies (temperature, shortwave radiation)
      call self%register_dependency(self%id_EIR,varname_swr)
      call self%register_dependency(self%id_ETW,varname_temp)

      return

99    call fatal_error('pml_ersem_vphyt1_create','Error reading namelist pml_ersem_vphyt1.')

100   call fatal_error('pml_ersem_vphyt1_create','Namelist pml_ersem_vphyt1 was not found.')

   end function pml_ersem_vphyt1_create

   
   subroutine do(self,_FABM_ARGS_DO_RHS_)

      class (type_pml_ersem_vphyt1),intent(in) :: self
      _DECLARE_FABM_ARGS_DO_RHS_

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
      _FABM_LOOP_BEGIN_

         _GET_(self%id_P1c,P1c)
         _GET_(self%id_P1p,P1p)
         _GET_(self%id_P1n,P1n)
         _GET_(self%id_P1s,P1s)
#ifdef IRON      
         _GET_(self%id_P1f,P1f)
         _GET_SAFE_(self%id_P1f,P1fP)
         _GET_SAFE_(self%id_N7f,N7fP)
#endif
         _GET_(self%id_Chl1,chl1)
         _GET_(self%id_N5s,N5s)

         _GET_SAFE_(self%id_P1c,P1cP)
         _GET_SAFE_(self%id_P1p,P1pP)
         _GET_SAFE_(self%id_P1n,P1nP)
         _GET_SAFE_(self%id_P1s,P1sP)
         _GET_SAFE_(self%id_Chl1,chl1P)
         _GET_SAFE_(self%id_N1p,N1pP)
         _GET_SAFE_(self%id_N3n,N3nP)
         _GET_SAFE_(self%id_N4n,N4nP)

         _GET_(self%id_ETW,ETW)
         _GET_(self%id_EIR,EIR)

#ifdef CENH      
         _GET_HORIZONTAL_(self%id_pco2a3,pco2a3)
#endif

         qpP1c = P1p/P1c
         qnP1c = P1n/P1c
         qsP1c = P1s/P1c
#ifdef IRON
         qfP1c = P1f/P1c
#endif

         iNP1p = MIN(1._rk,  &
                 MAX(0._rk, (qpP1c-self%qplP1cX) / (self%xqcP1pX*self%qpRPIcX-self%qplP1cX) ))
         iNP1n = MIN(1._rk,  &
                 MAX(0._rk, (qnP1c-self%qnlP1cX) / (self%xqcP1nX*self%qnRPIcX-self%qnlP1cX) ))
         iNP1s = MIN(1._rk, N5s/(N5s+self%chP1sX))
#ifdef IRON
         iNP1f = MIN(1._rk,  &
                 MAX(0._rk, (qfP1c-self%qflP1cX) / (self%qfRP1cX-self%qflP1cX) ))
#endif

         select case (self%LimnutX)
            case (0)
               iNIP1 = (iNP1p * iNP1n)**.5d0
            case (1)
               iNIP1 = MIN(iNP1p, iNP1n)
            case (2)
               iNIP1 = 2.0d0 / (1._rk/iNP1p + 1._rk/iNP1n)
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
         if (parEIR.gt.0._rk) THEN
            sumP1 = sumP1 * (1._rk-exp(-self%alphaP1X*parEIR*ChlCpp/sumP1)) * EXP(-self%betaP1X*parEIR*ChlCpp/sumP1)
            rho = (phi - self%ChlCmin) * (sumP1/(self%alphaP1X*parEIR*ChlCpp)) + self%ChlCmin
         else
            sumP1 = 0._rk
            rho = self%ChlCmin
         end if
#ifdef CENH
   ! Enhancement factor (from MEECE D1.5) 379.48 = pco2a @ 2005
         cenh=1.+(pco2a3-379.48)*0.0005
#endif

   !..Nutrient-stress lysis rate :
         sdoP1 = (1._rk/(MIN( iNP1s, iNIP1 )+0.1_rk))*self%sdoP1X

   !..Excretion rate, as regulated by nutrient-stress
         seoP1 = sumP1*(1._rk-iNIP1)*(1._rk-self%pu_eaP1X)

   !..activity-dependent excretion :
         seaP1 = sumP1*self%pu_eaP1X
         sugP1 = sumP1-seoP1-seaP1

   !..Apportioning of LOC- and DET- fraction of excretion/lysis fluxes:
         pe_R6P1 = MIN(self%qplP1cX/(qpP1c+0.0_rk),self%qnlP1cX/(qnP1c+0.0_rk))

         sP1R6 = pe_R6P1*sdoP1
         fP1R6c = sP1R6*P1cP

#ifdef DOCDYN
         fP1R1c=  (1._rk-pe_R6P1)*sdoP1*P1cP
         fP1R2c=  (seoP1 + seaP1)*P1c
         fP1RDc = fP1R1c+fP1R2c
#else
         fP1RDc = (1._rk-pe_R6P1)*sdoP1*P1cP + (seoP1 + seaP1)*P1c
#endif      
   !..Respiration..........................................................

   !..Rest respiration rate :
         srsP1 = etP1*self%srsP1X 

   !..Activity respiration rate :
         sraP1 = sugP1*self%pu_raP1X
#ifdef CENH
         sraP1 = sraP1*cenh
#endif

   !..Total respiration flux :
         fP1O3c = srsP1*P1cP+sraP1*P1c
         
   !..Gross production as flux from inorganic CO2 to P1 :
         fO3P1c = sumP1*P1c
#ifdef CENH
         fO3P1c = fO3P1c*cenh
#endif
         rugP1 = fO3P1c

   !..Production and productivity
         sunP1 = sumP1-(seoP1+seaP1+sraP1)        ! net productivity
         runP1 = sunP1*P1c-srsP1*P1cP       ! net production

   !..To save net production
         _SET_DIAGNOSTIC_(self%id_netP1,runP1)

   !..Carbon Source equations
        
         rho = MIN(rho,0.1_rk)
         ! Chl changes (note that Chl is a component of PXc and not involved
         ! in mass balance)
         Chl_inc = rho*(sumP1-sraP1)*P1c
         Chl_loss = (sdoP1+srsP1)*Chl1P + (seoP1+seaP1)*Chl1

         _SET_ODE_(self%id_P1c,(fO3P1c-fP1O3c-fP1R6c-fP1RDc)/secs_pr_day)
#ifdef DOCDYN
         _SET_ODE_(self%id_R1c,fP1R1c/secs_pr_day)
         _SET_ODE_(self%id_R2c,fP1R2c/secs_pr_day)
#else
         _SET_ODE_(self%id_R1c,(fP1RDc*self%R1R2X)/secs_pr_day)
         _SET_ODE_(self%id_R2c,fP1RDc*(1._rk-self%R1R2X)/secs_pr_day)
#endif
         _SET_ODE_(self%id_R6c,fP1R6c/secs_pr_day)
         _SET_ODE_(self%id_Chl1,(Chl_inc - Chl_loss)/secs_pr_day)

         _SET_ODE_(self%id_O3c,(fP1O3c - fO3P1c)/CMass/secs_pr_day)
         _SET_ODE_(self%id_O2o,(fO3P1c*self%uB1c_O2X - fP1O3c*self%urB1_O2X)/secs_pr_day)
         
   !..Phosphorus flux through P1...........................................

   !..Lysis loss of phosphorus
         fP1R6p = sP1R6 * min(self%qplP1cX*P1cP,P1pP)
         fP1RDp = sdoP1 * P1pP - fP1R6p

   !..Net phosphorus uptake
         rumP1p = self%qurP1pX * N1pP * P1c
         misP1p = self%xqpP1X * self%qpRPIcX*P1cP - P1pP
         runP1p = sunP1*P1c * self%qpRPIcX*self%xqpP1X - srsP1*P1pP
         fN1P1p = MIN(rumP1p, runP1p+misP1p)

   !..Source equations
         _SET_ODE_(self%id_P1p,(fN1P1p-fP1RDp-fP1R6p)/secs_pr_day)
         _SET_ODE_(self%id_N1p,-fN1P1p/secs_pr_day)
         _SET_ODE_(self%id_R6p,fP1R6p/secs_pr_day)
         _SET_ODE_(self%id_R1p,fP1RDp/secs_pr_day)

   !..Nitrogen flux through P1.............................................

   !..Nitrogen loss by lysis
         fP1R6n = sP1R6 * min(self%qnlP1cX*P1cP,P1nP)
         fP1RDn = sdoP1 * P1nP - fP1R6n

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
         _SET_ODE_(self%id_R6n,fP1R6n/secs_pr_day)
         _SET_ODE_(self%id_R1n,fP1RDn/secs_pr_day)

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

#ifdef IRON
   !..Iron flux through P1................................................

   !..Iron loss by lysis
   !  Because its high affinity with particles all the iron lost from phytoplankton by lysis is supposed to be 
   !  associated to organic particulate detritus. (luca)

         fP1R6f = sdoP1 * P1fP

   !..net iron uptake
         rumP1f = self%qurP1fX * N7fP * P1c
         misP1f = self%qfRP1cX*P1cP - P1fP
         runP1f = sunP1*P1c * self%qfRP1cX - srsP1*P1fP
         fN7P1f = MIN(rumP1f, runP1f+misP1f)

   !..Source equations
         _SET_ODE_(self%id_P1f,(fN7P1f-fP1R6f)/secs_pr_day)
         _SET_ODE_(self%id_N7f,-fN7P1f/secs_pr_day)
         _SET_ODE_(self%id_R6f,fP1R6f/secs_pr_day)
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
      _FABM_LOOP_END_

   end subroutine do

   subroutine get_vertical_movement(self,_FABM_ARGS_GET_VERTICAL_MOVEMENT_)
      class (type_pml_ersem_vphyt1),intent(in) :: self
      _DECLARE_FABM_ARGS_GET_VERTICAL_MOVEMENT_
      
      real(rk) :: SDP1
      real(rk) :: P1c,P1p,P1n,qpP1c,qnP1c,iNP1p,iNP1n,iNIP1

      _FABM_LOOP_BEGIN_

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
               iNIP1 = (iNP1p * iNP1n)**.5d0
            case (1)
               iNIP1 = MIN(iNP1p, iNP1n)
            case (2)
               iNIP1 = 2.0d0 / (1._rk/iNP1p + 1._rk/iNP1n)
         end select

!..sedimentation and resting stages.....................................
         SDP1 = -self%resP1mX * MAX(0._rk, (self%esNIP1X - iNIP1))/secs_pr_day
         _SET_VERTICAL_MOVEMENT_(self%id_P1c,SDP1)
         _SET_VERTICAL_MOVEMENT_(self%id_P1p,SDP1)
         _SET_VERTICAL_MOVEMENT_(self%id_P1n,SDP1)
         _SET_VERTICAL_MOVEMENT_(self%id_P1s,SDP1)
         _SET_VERTICAL_MOVEMENT_(self%id_Chl1,SDP1)
#ifdef IRON
         _SET_VERTICAL_MOVEMENT_(self%id_P1f,SDP1)
#endif

      _FABM_LOOP_END_

   end subroutine get_vertical_movement
   
end module