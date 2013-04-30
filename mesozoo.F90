!-----------------------------------------------------------------------
! TODO PML-ERSEM, version 0.9
! Plymouth Marine Laboratory
! Prospect Place, The Hoe
! Plymouth, PL1 3DH, UK
!    contact: momm@pml.ac.uk
!-----------------------------------------------------------------------
#include"ppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mesozoo \label{sec:mesozoo}
!
! !DESCRIPTION:
!  TODO - Descripton
!\\
!\\
! !INTERFACE:
   MODULE mesozoo
!
! !USES:
   use global_declarations
   use general, ONLY: Adjust_Fixed_Nutrients
   use gas_dynamics, ONLY: ub1c_o2x,urB1_O2X
!
   IMPLICIT NONE
!
!  Default all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public mesozoop
!
! !PUBLIC DATA MEMBERS:
   real(fp8), public :: chuz4cX,sumz4X,puz4X,q10z4X,srsz4X,pu_eaz4X,pu_earz4X
   real(fp8), public :: minfoodz4X,sdz4oX,sdz4X,pe_r1z4X,chrz4oX,Z4repwX,Z4mortX
!
! !PRIVATE DATA MEMBERS:
   real(fp8) :: fZ4R8n,fZ4R8p,fZ4Rdc,fZ4Rdn,fZ4Rdp,ineffZ4
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
! !TO DO:
!  Document this module.
!
!EOP
!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mesozoop \label{sec:mesozoop}
!
! !DESCRIPTION:
!  TODO - description
!
!  simple mesozooplankton routine
!\\
!\\
! !INTERFACE:
   SUBROUTINE mesozoop(I)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer :: i
!
! !LOCAL VARIABLES:
   real(fp8) :: SZ4n, SZ4p, etZ4
   real(fp8) :: rumZ4, put_uZ4, rrsZ4, rraZ4, rugZ4, rdZ4, sdZ4
   real(fp8) :: sP1Z4,sP2Z4,sP3Z4,sP4Z4,sB1Z4,sR6Z4,sZ5Z4,sZ6Z4
   real(fp8) :: corrox, eO2Z4, fZ4RIp, fZ4RIn
   real(fp8) :: retZ4
   real(fp8) :: ruP1Z4c, ruP2Z4c, ruP3Z4c, ruP4Z4c
   real(fp8) :: ruZ6Z4c, ruB1Z4c, ruZ5Z4c, ruR6Z4c, ruZ4Z4c
   real(fp8) :: temp_n,temp_p
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
! !TO DO:
!  Document this module.
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
!..set up dummy source terms for nutrient components....................
       SZ4n = SZ4c(I) * qnZIcX
       SZ4p = SZ4c(I) * qpZIcX

!..Temperature effect :

       etZ4 = q10Z4X**((ETW(I)-10._fp8)/10._fp8) - q10Z4X**((ETW(I)-32._fp8)/3._fp8)

!..Oxygen limitation :
       CORROX = 1._fp8 + chrZ4oX
       eO2Z4 = MIN(1._fp8,CORROX*(eO2mO2/(chrZ4oX + eO2mO2)))

!..Available food :
       sB1Z4 = suB1_Z4X*B1cP(I)/(B1cP(I)+minfoodZ4X)
       sP1Z4 = suP1_Z4X*P1cP(I)/(P1cP(I)+minfoodZ4X)
       sP2Z4 = suP2_Z4X*P2cP(I)/(P2cP(I)+minfoodZ4X)
       sP3Z4 = suP3_Z4X*P3cP(I)/(P3cP(I)+minfoodZ4X)
       sP4Z4 = suP4_Z4X*P4cP(I)/(P4cP(I)+minfoodZ4X)
       sZ5Z4 = suZ5_Z4X*Z5cP(I)/(Z5cP(I)+minfoodZ4X)
       sZ6Z4 = suZ6_Z4X*Z6cP(I)/(Z6cP(I)+minfoodZ4X)
       sR6Z4 = suR6_Z4X*R6cP(I)/(R6cP(I)+minfoodZ4X)
       ruZ4Z4c = suZ4_Z4X*Z4cP(I)*Z4cP(I)/(Z4cP(I)+minfoodZ4X)
       ruB1Z4c = sB1Z4*B1cP(I)
       ruP1Z4c = sP1Z4*P1cP(I)
       ruP2Z4c = sP2Z4*P2cP(I)
       ruP3Z4c = sP3Z4*P3cP(I)
       ruP4Z4c = sP4Z4*P4cP(I)
       ruZ5Z4c = sZ5Z4*Z5cP(I)
       ruZ6Z4c = sZ6Z4*Z6cP(I)
       ruR6Z4c = sR6Z4*R6cP(I)

       rumZ4 = ruP1Z4c + ruP2Z4c + ruP3Z4c + ruP4Z4c + ruZ4Z4c &
     &       + ruB1Z4c + ruZ5Z4c + ruZ6Z4c + ruR6Z4c

!..Uptake :
       put_uZ4 = sumZ4X/(rumZ4 + chuZ4cX)*etZ4*Z4c(I)
       rugZ4 = put_uZ4*rumZ4

!..Fluxes into mesoplankton :
       fB1Z4c(I) = put_uZ4*ruB1Z4c
       fP1Z4c(I) = put_uZ4*ruP1Z4c
       fP2Z4c(I) = put_uZ4*ruP2Z4c
       fP3Z4c(I) = put_uZ4*ruP3Z4c
       fP4Z4c(I) = put_uZ4*ruP4Z4c
       fZ5Z4c(I) = put_uZ4*ruZ5Z4c
       fZ6Z4c(I) = put_uZ4*ruZ6Z4c
       fR6Z4c(I) = put_uZ4*ruR6Z4c
       fZ4Z4c(I) = put_uZ4*ruZ4Z4c
       sB1Z4 = put_uZ4*sB1Z4
       sP1Z4 = put_uZ4*sP1Z4
       sP2Z4 = put_uZ4*sP2Z4
       sP3Z4 = put_uZ4*sP3Z4
       sP4Z4 = put_uZ4*sP4Z4
       sZ5Z4 = put_uZ4*sZ5Z4
       sZ6Z4 = put_uZ4*sZ6Z4
       sR6Z4 = put_uZ4*sR6Z4

!..Zooplankton Grazing
!      Z4herb(i) = fP1Z4c(i) + fP2Z4c(i) + fP3Z4c(i) +fP4Z4c(i)
!      Z4carn(i) = fB1Z4c(i) + fZ5Z4c(i) + fZ6Z4c(i) + fR6Z4c(i)

!..Mortality
       sdZ4 = ((1._fp8 - eO2Z4)*sdZ4oX + sdZ4X) 
       rdZ4 = sdZ4*Z4cP(I)

!..Assimilation inefficiency:
      ineffZ4 = (1._fp8 - puZ4X)

!..Excretion
       retZ4 = ineffZ4 * ( (rugZ4-fR6Z4c(I)) * pu_eaZ4X &
             + fR6Z4c(I) * pu_eaRZ4X )
       fZ4RDc = (retZ4 + rdZ4)*pe_R1Z4X
       fZ4R8c(I) = (retZ4 + rdZ4)*(1._fp8 - pe_R1Z4X)
#ifdef SAVEFLX
       fZXRDc(I) = fZXRDc(I)+fZ4RDc
#endif

!..Rest respiration, corrected for prevailing temperature
       rrsZ4 = srsZ4X*etZ4*Z4cP(I)

!..Activity respiration
       rraZ4 = ineffZ4 * rugZ4 - retZ4

!..Total respiration
       fZ4O3c(I) = rrsZ4 + rraZ4

#ifdef CALC
       fO3L2c(I) = fO3L2c(I) + &
         RainR(I) * gutdiss * (1._fp8 - puZ4X)* pu_eaZ4X * fP2Z4c(I)
#endif
       rraZ4 = rugZ4*(1._fp8 - puZ4X)-retZ4

!..Source equations
       SZ4c(I) = SZ4c(I) + rugZ4 &
                         - fZ4RDc - fZ4R8c(I) - fZ4Z4c(I) - fZ4O3c(I)

!..Flows from and to detritus
       SR1c(I) = SR1c(I) + (fZ4RDc * R1R2X)
       SR2c(I) = SR2c(I) + (fZ4RDc * (1._fp8-R1R2X))
       SR8c(I) = SR8c(I) + fZ4R8c(I)
       SR6c(I) = SR6c(I) - fR6Z4c(I)

!..Grazing and predation
       SP1c(I) = SP1c(I) - fP1Z4c(I)
       SP2c(I) = SP2c(I) - fP2Z4c(I)
       SP3c(I) = SP3c(I) - fP3Z4c(I)
       SP4c(I) = SP4c(I) - fP4Z4c(I)
       SZ5c(I) = SZ5c(I) - fZ5Z4c(I)
       SZ6c(I) = SZ6c(I) - fZ6Z4c(I)
       SB1c(I) = SB1c(I) - fB1Z4c(I)
       SChl1(I) = SChl1(I) - sP1Z4*chl1P(I)
       SChl2(I) = SChl2(I) - sP2Z4*chl2P(I)
       SChl3(I) = SChl3(I) - sP3Z4*chl3P(I)
       SChl4(I) = SChl4(I) - sP4Z4*chl4P(I)

!..Respiration
       SO3c(I) = SO3c(I) + fZ4O3c(I)/CMass
       SO2o(I) = SO2o(I) - fZ4O3c(I)*urB1_O2X

!..Phosphorus dynamics in mesozooplankton, derived from carbon flows....

       fZ4RIp = (fZ4RDc + fZ4R8c(I)) * qpZIcX
       fZ4RDp = min(fZ4RIp, fZ4RDc * qpZIcX * xR1pX)
       fZ4R8p = fZ4RIp - fZ4RDp

!..Source equations
       SZ4p = SZ4p + sP1Z4*P1pP(I) &
                   + sP2Z4*P2pP(I) &
                   + sP3Z4*P3pP(I) &
                   + sP4Z4*P4pP(I) &
                   + sZ5Z4*Z5pP(I) &
                   + sZ6Z4*Z6pP(I) &
                   + sB1Z4*B1pP(I) &
                   + sR6Z4*R6pP(I) &
                   - fZ4R8p - fZ4RDp

#ifdef IRON
!  Iron dynamics
! following Vichi et al., 2007 it is assumed that the iron fraction of the ingested phytoplankton
! is egested as particulate detritus (Luca)

   SR6f(I)=SR6f(I) + sP1Z4*P1fP(I)+sP4Z4*P4fP(I)-sR6Z4*R6fP(I)
   SR4f(I)=SR4f(I) + sP2Z4*P2fP(I)+sP3Z4*P3fP(I)
   SN7f(I)=SN7f(I)+sR6Z4*R6fP(I)

   SP1f(I)=SP1f(I)- sP1Z4*P1fP(I)
   SP4f(I)=SP4f(I)- sP4Z4*P4fP(I)
   SP2f(I)=SP2f(I)- sP2Z4*P2fP(I)
   SP3f(I)=SP3f(I)- sP3Z4*P3fP(I) 
#endif

!..Phosphorus flux from/to detritus
       SR1p(I) = SR1p(I) + fZ4RDp
       SR6p(I) = SR6p(I) - sR6Z4*R6pP(I)
       SR8p(I) = SR8p(I) + fZ4R8p

!..Phosphorus flux from prey
       SP1p(I) = SP1p(I) - sP1Z4*P1pP(I)
       SP2p(I) = SP2p(I) - sP2Z4*P2pP(I)
       SP3p(I) = SP3p(I) - sP3Z4*P3pP(I)
       SP4p(I) = SP4p(I) - sP4Z4*P4pP(I)
       SZ5p(I) = SZ5p(I) - sZ5Z4*Z5pP(I)
       SZ6p(I) = SZ6p(I) - sZ6Z4*Z6pP(I)
       SB1p(I) = SB1p(I) - sB1Z4*B1pP(I)

!..Nitrogen dynamics in mesozooplankton, derived from carbon flows......

       fZ4RIn = (fZ4RDc + fZ4R8c(I)) * qnZIcX
       fZ4RDn = min(fZ4RIn, fZ4RDc * qnZIcX * xR1nX)
       fZ4R8n = fZ4RIn - fZ4RDn

!..Source equations
       SZ4n = SZ4n + sP1Z4*P1nP(I) &
                   + sP2Z4*P2nP(I) &
                   + sP3Z4*P3nP(I) &
                   + sP4Z4*P4nP(I) &
                   + sZ5Z4*Z5nP(I) &
                   + sZ6Z4*Z6nP(I) &
                   + sB1Z4*B1nP(I) &
                   + sR6Z4*R6nP(I) &
                   - fZ4R8n - fZ4RDn

!..Nitrogen flux from/to detritus
       SR1n(I) = SR1n(I) + fZ4RDn
       SR6n(I) = SR6n(I) - sR6Z4*R6nP(I)
       SR8n(I) = SR8n(I) + fZ4R8n

!..Nitrogen flux from prey
       SP1n(I) = SP1n(I) - sP1Z4*P1nP(I)
       SP2n(I) = SP2n(I) - sP2Z4*P2nP(I)
       SP3n(I) = SP3n(I) - sP3Z4*P3nP(I)
       SP4n(I) = SP4n(I) - sP4Z4*P4nP(I)
       SZ5n(I) = SZ5n(I) - sZ5Z4*Z5nP(I)
       SZ6n(I) = SZ6n(I) - sZ6Z4*Z6nP(I)
       SB1n(I) = SB1n(I) - sB1Z4*B1nP(I)

!..Silica-flux from diatoms due to mesozooplankton grazing
       SP1s(I) = SP1s(I) - sP1Z4 * P1sP(I)
       SR8s(I) = SR8s(I) + sP1Z4 * P1sP(I)

!..re-establish the fixed nutrient ratio in zooplankton.................

    temp_p = SN1p(I)
    temp_n = SN4n(I)
    CALL Adjust_fixed_nutrients ( SZ4c(I), SZ4n, SZ4p, qnZIcX, &
                                qpZIcX, SN4n(I), SN1p(I), SR8c(I) )
    fZ4N1p(I) = SN1p(I)-temp_p
    fZ4NIn(I) = SN4n(I)-temp_n
    
   RETURN

   END SUBROUTINE mesozoop
!
!EOC
!-----------------------------------------------------------------------

   END MODULE mesozoo

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
