#include "fabm_driver.h"
module ersem_carbonate

   use fabm_types
   use fabm_builtin_models
   use ersem_shared

   implicit none

   private

   type,extends(type_base_model),public :: type_ersem_carbonate
!     Variable identifiers
      type (type_state_variable_id)     :: id_O3c,id_TA,id_bioalk
      type (type_dependency_id)         :: id_ETW, id_X1X, id_dens, id_pres
      type (type_dependency_id)         :: id_Carb_in,id_pco2_in
      type (type_horizontal_dependency_id) :: id_wnd,id_PCO2A

      type (type_diagnostic_variable_id) :: id_ph,id_pco2,id_CarbA, id_BiCarb, id_Carb, id_TA_diag
      type (type_diagnostic_variable_id) :: id_Om_cal,id_Om_arg

      type (type_horizontal_diagnostic_variable_id) :: id_fairmg

      integer :: iswCO2X,iswtalk,iswASFLUX
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
   end type

   public :: co2dyn, CaCO3_Saturation

contains

   subroutine initialize(self,configunit)
!
! !INPUT PARAMETERS:
      class (type_ersem_carbonate), intent(inout), target :: self
      integer,                      intent(in)            :: configunit

      integer :: iswbioalk
!
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%iswCO2X,'iswCO2','','carbonate system diagnostics (0: off, 1: on)',default=1)
      call self%get_parameter(self%iswASFLUX,'iswASFLUX','','air-sea CO2 exchange (0: none, 1: Nightingale and Liss, 2: Wanninkhof 1992 without chemical enhancement, 3: Wanninkhof 1992 with chemical enhancement, 4: Wanninkhof 1999)',default=1)
      call self%get_parameter(self%iswtalk,'iswtalk','','alkalinity formulation (1-4: from salinity and temperature, 5: dynamic alkalinity)',default=5)
      if (self%iswtalk<1.or.self%iswtalk>5) call self%fatal_error('initialize','iswtalk takes values between 1 and 5 only')

      call self%register_state_variable(self%id_O3c,'c','mmol C/m^3','total dissolved inorganic carbon', 2200._rk,minimum=0._rk)
      call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_O3c)

      if (self%iswtalk==5) then
         ! Total alkalinity is a state variable.
         call self%register_state_variable(self%id_TA,'TA','umol/kg','total alkalinity',2300._rk,minimum=1.e-4_rk, &
            standard_variable=standard_variables%alkalinity_expressed_as_mole_equivalent)
      else
         ! Total alkalinity is a diagnostic variable, parameterized as function of salinity and temperature.
         call self%register_diagnostic_variable(self%id_TA_diag,'TA','umol/kg','total alkalinity', act_as_state_variable=.true., &
            standard_variable=standard_variables%alkalinity_expressed_as_mole_equivalent)

         call self%get_parameter(iswbioalk,'iswbioalk','','use bioalkalinity (0: off, 1: on)',default=1)
         if (iswbioalk==1) then
            ! Register state variable to track "bioalkalinity", i.e., the difference between parameterized
            ! and actual alkalinity that is created by biogeochemical processes modifying alkalinity.
            call self%register_state_variable(self%id_bioalk,'bioalk','umol/kg','bioalkalinity')

            ! Redirect all source terms and boundary fluxes imposed on total alkalinity [a diagnostic]
            ! to bioalkalinity [a state variable]
            call copy_fluxes(self,self%id_TA_diag,self%id_bioalk)
         end if
      end if

      if (self%iswCO2X==1) then
         call self%register_diagnostic_variable(self%id_ph,    'pH',    '-',      'pH',standard_variable=standard_variables%ph_reported_on_total_scale)
         call self%register_diagnostic_variable(self%id_pco2,  'pCO2',  '1e-6',    'partial pressure of CO2')
         call self%register_diagnostic_variable(self%id_CarbA, 'CarbA', 'mmol/m^3','carbonic acid concentration')
         call self%register_diagnostic_variable(self%id_BiCarb,'BiCarb','mmol/m^3','bicarbonate concentration')
         call self%register_diagnostic_variable(self%id_Carb,  'Carb',  'mmol/m^3','carbonate concentration',standard_variable=standard_variables%mole_concentration_of_carbonate_expressed_as_carbon)

         call self%register_diagnostic_variable(self%id_Om_cal,'Om_cal','-','calcite saturation')
         call self%register_diagnostic_variable(self%id_Om_arg,'Om_arg','-','aragonite saturation')
      end if
      if (self%iswASFLUX>=1) then
         call self%register_diagnostic_variable(self%id_fairmg,'fairmg','mmol C/m^2/d','air-sea flux of CO2')
      end if

      call self%register_dependency(self%id_ETW, standard_variables%temperature)
      call self%register_dependency(self%id_X1X, standard_variables%practical_salinity)
      call self%register_dependency(self%id_dens,standard_variables%density)
      call self%register_dependency(self%id_pres,standard_variables%pressure)
      call self%register_dependency(self%id_pco2_in,'pCO2','1e-6','previous pCO2')
      call self%register_dependency(self%id_Carb_in,'Carb','mmol/m^3','previous carbonate concentration')

      if (self%iswASFLUX>=1) then
         call self%register_dependency(self%id_wnd,  standard_variables%wind_speed)
         call self%register_dependency(self%id_PCO2A,standard_variables%mole_fraction_of_carbon_dioxide_in_air)
      end if

      self%dt = 3600._rk*24._rk

   end subroutine

   function approximate_alkalinity(iswtalk,T,S) result(TA)
      integer, intent(in) :: iswtalk
      real(rk),intent(in) :: T,S
      real(rk)            :: TA

      ! NB total alkalinity must be given in umol kg-1

      if (iswtalk==1) then       !!!!   Old formulation
         IF (S .GE. 34.65_rk) THEN
            TA =  66.96_rk*S - 36.803_rk   ! Bellerby & Canoba
         ELSE
            TA = 3887._rk - 46.25_rk*S     ! Borges & Frankignoulle & Canoba
         END IF
      elseif (iswtalk==2) then
         TA = 520.1_rk + 51.24_rk*S        ! from Millero et al, MarChem, 1998
      elseif (iswtalk==3) then
         IF (T.lt.20._rk) then
            TA = S/35._rk*(2291._rk-2.69_rk*(T-20._rk)-0.046_rk*(T-20._rk)**2)   ! Millero et al, MarChem, 1998
         else
            TA = 520.1_rk+51.24_rk*S                                             ! Millero et al, MarChem, 1998
         ENDIF
      elseif (iswtalk==4) then
         IF (T.lt.20._rk) then
            TA = 2305._rk+53.97_rk*(S-35._rk)+2.74_rk*(S-35._rk)**2-1.16_rk*(T-20._rk)-0.04_rk*(T-20._rk)**2   ! Lee et al., Geophys Res Lett, 1998
         else
            TA = 2305._rk+58.66_rk*(S-35._rk)+2.32_rk*(S-35._rk)**2-1.41_rk*(T-20._rk)+0.04_rk*(T-20._rk)**2   ! Lee et al., Geophys Res Lett, 1998
         ENDIF
      end if
   end function

   subroutine do(self,_ARGUMENTS_DO_)
      class (type_ersem_carbonate), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: O3c,ETW,X1X,density,pres
      real(rk) :: TA,bioalk,Ctot
      real(rk) :: pH,PCO2,H2CO3,HCO3,CO3,k0co2
      real(rk) :: Om_cal,Om_arg
      logical  :: success

      IF (self%iswCO2X .NE. 1 ) RETURN

      _LOOP_BEGIN_
         _GET_(self%id_O3C,O3C)
         _GET_(self%id_ETW,ETW)
         _GET_(self%id_X1X,X1X)
         _GET_(self%id_dens,density)
         _GET_(self%id_pres,pres)

         ! Calculate total alkalinity
         if (self%iswtalk/=5) then
            ! Alkalinity is parameterized as function of salinity and temperature.
            TA = approximate_alkalinity(self%iswtalk,ETW,X1X)
            if (_VARIABLE_REGISTERED_(self%id_bioalk)) then
               ! We separately track bioalkalinity - include this in total alkalinity.
               _GET_(self%id_bioalk,bioalk)
               TA = TA + bioalk
            end if
            _SET_DIAGNOSTIC_(self%id_TA_diag,TA) ! in umol/kg
         else
            ! Alkalinity is a state variable.
            _GET_(self%id_TA,TA)
         end if

         TA = TA/1.0e6_rk                ! from umol kg-1 to mol kg-1
         Ctot  = O3C / 1.e3_rk / density ! from mmol m-3 to mol kg-1

         CALL CO2DYN (ETW,X1X,pres*0.1_rk,ctot,TA,pH,PCO2,H2CO3,HCO3,CO3,k0co2,success)   ! NB pressure from dbar to bar

         if (success) then
            ! Carbonate system iterative scheme converged -  save associated diagnostics.
            ! Convert outputs from fraction to ppm (pCO2) and from mol kg-1 to mmol m-3 (concentrations).
            _SET_DIAGNOSTIC_(self%id_ph,pH)
            _SET_DIAGNOSTIC_(self%id_pco2,PCO2*1.e6_rk)
            _SET_DIAGNOSTIC_(self%id_CarbA, H2CO3*1.e3_rk*density)
            _SET_DIAGNOSTIC_(self%id_Bicarb,HCO3*1.e3_rk*density)
            _SET_DIAGNOSTIC_(self%id_Carb,  CO3*1.e3_rk*density)
         else
            ! Carbonate system iterative scheme did not converge.
            ! All diagnostics retain their previous value.
            ! Use previous carbonate concentration (but current environment) for carbonate saturation states.
            _GET_(self%id_Carb_in,CO3)
            CO3 = CO3/1.e3_rk/density  ! from mmol/m3 to mol/kg
         end if

         ! Call carbonate saturation state subroutine
         CALL CaCO3_Saturation (ETW, X1X, pres*1.e4_rk, CO3, Om_cal, Om_arg)  ! NB pressure from dbar to Pa

         _SET_DIAGNOSTIC_(self%id_Om_cal,Om_cal)
         _SET_DIAGNOSTIC_(self%id_Om_arg,Om_arg)
      _LOOP_END_
   end subroutine

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_ersem_carbonate), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: O3c,T,S,PRSS,density,bioalk
      real(rk) :: wnd,PCO2A
      real(rk) :: sc,fwind,UPTAKE,FAIRCO2

      real(rk) :: ctot,TA,pH,PCO2,H2CO3,HCO3,CO3,k0co2
      logical  :: success

      if (self%iswASFLUX<=0) return

      _HORIZONTAL_LOOP_BEGIN_
         _GET_(self%id_O3c,O3c)
         _GET_(self%id_ETW,T)
         _GET_(self%id_X1X,S)
         _GET_(self%id_pres,PRSS)
         _GET_(self%id_dens,density)
         _GET_HORIZONTAL_(self%id_wnd,wnd)
         _GET_HORIZONTAL_(self%id_PCO2A,PCO2A)

         if (self%iswtalk/=5) then
            ! Alkalinity is parameterized as function of salinity and temperature.
            TA = approximate_alkalinity(self%iswtalk,T,S)
            if (_VARIABLE_REGISTERED_(self%id_bioalk)) then
               ! We separately track bioalkalinity - include this in total alkalinity.
               _GET_(self%id_bioalk,bioalk)
               TA = TA + bioalk
            end if
         else
            ! Alkalinity is a state variable.
            _GET_(self%id_TA,TA)
         end if

         TA = TA/1.0e6_rk                ! from umol kg-1 to mol kg-1
         Ctot  = O3C / 1.e3_rk / density ! from mmol m-3 to mol kg-1

!  for surface box only calculate air-sea flux
!..Only call after 2 days, because the derivation of instability in the
!..
         CALL CO2dyn(T, S, PRSS*0.1_rk,ctot,TA,pH,PCO2,H2CO3,HCO3,CO3,k0co2,success)
         if (.not.success) then
            _GET_(self%id_pco2_in,PCO2)
            PCO2 = PCO2*1.e-6_rk
         end if
         if (legacy_ersem_compatibility) then
           ! Old formulation for the Schmidt number, valid only for T<30
           ! left for documentation and back compatibility
           sc=2073.1_rk-125.62_rk*T+3.6276_rk*T**2._rk-0.043219_rk*T**3
         else
           !new formulation for the Schmidt number following Wanninkof, 2014
           T=max(T,40._rk)
           sc=2116.8_rk-136.25_rk*T+4.7353_rk*T**2._rk-0.092307_rk*T**3+0.0007555_rk*T**4._rk
         endif
         if   (self%ISWASFLUX==1) then        !Nightingale and Liss parameterisation
               fwind =  (0.222_rk *wnd**2.0_rk + 0.333_rk * wnd)*(sc/600._rk)**(-0.5_rk)
         elseif (self%iswASFLUX==2) then      ! Wanninkhof 1992 without chemical enhancement
               fwind=0.31_rk*wnd**2.0_rk*(sc/660._rk)**(-0.5_rk)
         elseif (self%iswASFLUX==3)THEN       ! Wanninkhof 1992 with chemical enhancement
               fwind=(2.5_rk*(0.5246_rk+0.016256_rk*T+0.00049946_rk*T**2.0_rk)+0.3_rk*wnd**2.0_rk)*(sc/660._rk)**(-0.5_rk)
         elseif(self%iswASFLUX==4)THEN             ! Wanninkhof 1999 
              fwind=0.0283_rk*wnd**3.0_rk*(sc/660._rk)**(-0.5_rk)
         endif

         fwind=fwind*24._rk/100._rk   ! convert to m/day
         UPTAKE = fwind * k0co2 * ( PCO2A/1.e6_rk - PCO2 )

         FAIRCO2 = UPTAKE * 1.e3_rk * density

         _SET_SURFACE_EXCHANGE_(self%id_O3c,FAIRCO2)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fairmg,FAIRCO2)
      _HORIZONTAL_LOOP_END_
   end subroutine

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CO2dyn \label{sec:CO2dyn}
!
! !DESCRIPTION:
!  TODO - description
!
!     This subroutine calculates the partial pressure of CO2 (pCO2) at
!     the ambient salinity, temperature, alkalinity and total CO2 and
!     hence the CO2 exchange across the air sea interface.
!\\
!\\
! !INTERFACE:
      SUBROUTINE CO2dyn ( T, S, PRSS,ctot,TA,pH,PCO2,H2CO3,HCO3,CO3,k0co2,success)
!
! !LOCAL VARIABLES:
!     ! TODO - SORT THESE!
      real(rk),intent(in) :: T, S, PRSS   ! NB PRSS is pressure in bar
      real(rk),intent(inout) :: ctot,TA
      real(rk),intent(out) :: pH,PCO2,H2CO3,HCO3,CO3,k0co2
      logical, intent(out) :: success

      real(rk) :: k1co2,k2co2,kb
      real(rk) :: Tmax, Btot
      INTEGER  :: ICALC
      LOGICAL  :: BORON
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC

      ICALC = 1               ! use DIC and TA to calculate CO2sys
      Tmax = max(T,0._rk)

      BORON=.True.
      IF(BORON) THEN
        BTOT=0.0004128_rk*S/35._rk
      ENDIF

      CALL CO2SET(PRSS,Tmax,S,k0co2,k1co2,k2co2,kb)
      CALL CO2CLC(k0co2,k1co2,k2co2,kb,ICALC,BORON,BTOT,ctot,TA,pH,PCO2,H2CO3,HCO3,CO3,success)

      END SUBROUTINE CO2DYN
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CO2SET \label{sec:CO2SET}
!
! !DESCRIPTION:
!  TODO - CHECK THIS!
!
! Routine to calculate CO2 system constants under the conditions set by
! P,S,T     (NOTE: PRESSURE <> 1ATM IS NOT YET CODED)
! I. Calculate constants at P=1 and S=0 using
!      ln K0  =  A + B/TK + C ln TK
!                                   (where TK is in Kelvin)
! II. Calculate constants at P=1 and salinity S using
!    ln K  =  ln K0 + (a0 + a1/TK + a2 ln TK) S**1/2
!                 + (b0 + b1TK + b2TK**2) S
! The sources of the coefficients are as follows:
!  IC=                  1                    2                  3
!               (NBS pH scale)        (SWS pH scale       (SWS pH scale
!                                       with no F)
!  KP            WEISS (1974)           WEISS(1974)          WEISS(1974
!  K1C )      MEHRBACH ACC. TO       HANSSON ACC. TO    HANSSON AND MEH
!  K2C )       MILLERO (1979)         MILLERO (1979)      ACC. TO DICKS
!   KB )                                                 AND MILLERO (1
!                                                         (K1C AND K2C
!                                                           HANSSON ACC
!                                                            MILLERO (1
!                                                              (KB ONLY
!\\
!\\
! !INTERFACE:
      SUBROUTINE CO2SET(P,T,S,k0co2,k1co2,k2co2,kb)
         real(rk),intent(in)  :: P,T,S
         real(rk),intent(out) :: k0co2,k1co2,k2co2,kb
!
! !USES:
!
! !LOCAL VARIABLES
! 
      real(rk)              :: TK, delta, kappa
      real(rk)              :: dlogTK, S2, S15, sqrtS,TK100
      real(rk),parameter    :: Rgas = 83.131_rk

!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!

!  Derive simple terms used more than once
       TK=T+273.15_rk
       dlogTK = log(TK)
       TK100=TK/100._rk
       S2 = S*S
       sqrtS = sqrt(S)
       S15 = S**1.5_rk

!  Calculation of constants as used in the OCMIP process for ICONST = 3 or 6
!  see http://www.ipsl.jussieu.fr/OCMIP/
! k0co2 = CO2/fCO2 from Weiss 1974
        k0co2 = exp(93.4517/tk100 - 60.2409 + 23.3585 * log(tk100) + &
        &       s * (.023517 - 0.023656 * tk100 + 0.0047036 * tk100 ** 2._rk))

! the commented definition of K0co2 herebelow is the one giving bit-identical output to the old CO2sys module
!        VAL=-167.8108_rk + 9345.17_rk/TK + 23.3585_rk*LOG(TK)
!        VAL=VAL + (0.023517_rk-2.3656e-4_rk*TK+4.7036e-7_rk*TK*TK)*S
!        k0co2=EXP(VAL)
! end of old definition

! Correction for high pressure of K0 (from Millero 1995)
! added YA 04/10/2010
        delta=-25.6_rk+0.2324_rk*T-0.0036246_rk*T**2
        kappa=(-5.13_rk+0.0794_rk*T)/1000._rk
        k0co2=k0co2*exp((-delta+0.5_rk*kappa*P)*P/(Rgas*TK))

! k1co2 = [H][HCO3]/[H2CO3]
! Millero p.664 (1995) using Mehrbach et al. data on seawater scale 
        k1co2=10._rk**(-1._rk*(3670.7_rk/TK - 62.008_rk + 9.7944_rk*dlogTK - &
     &     0.0118_rk * S + 0.000116_rk*S2))
! Correction for high pressure (from Millero 1995)
! added YA 04/10/2010
        delta=-25.5_rk+0.1271_rk*T
        kappa=(-3.08_rk+0.0877_rk*T)/1000._rk
        k1co2=k1co2*exp((-delta+0.5_rk*kappa*P)*P/(Rgas*TK))

! k2co2 = [H][CO3]/[HCO3]
! Millero p.664 (1995) using Mehrbach et al. data on seawater scale 
        k2co2=10._rk**(-1._rk*(1394.7_rk/TK + 4.777_rk - &
     &     0.0184_rk*S + 0.000118_rk*S2))
! Correction for high pressure (from Millero 1995)
! added YA 04/10/2010
        delta=-15.82_rk-0.0219_rk*T
        kappa=(1.13_rk-0.1475_rk*T)/1000._rk
        k2co2=k2co2*exp((-delta+0.5_rk*kappa*P)*P/(Rgas*TK))


! kb = [H][BO2]/[HBO2]
! Millero p.669 (1995) using data from Dickson (1990)
        kb=exp((-8966.9_rk - 2890.53_rk*sqrtS - 77.942_rk*S + &
     &      1.728_rk*S15 - 0.0996_rk*S2)/TK + &
     &      (148.0248_rk + 137.1942_rk*sqrtS + 1.62142_rk*S) + &
     &      (-24.4344_rk - 25.085_rk*sqrtS - 0.2474_rk*S) * &
     &      dlogTK + 0.053105_rk*sqrtS*TK)
! Correction for high pressure (from Millero, 1995)
! added YA 04/10/2010
        delta=-29.48_rk+0.1622_rk*T-0.002608_rk*T**2._rk
        kappa=-2.84_rk/1000._rk
        kb=kb*exp((-delta+0.5_rk*kappa*P)*P/(Rgas*TK))

      END SUBROUTINE CO2SET
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CO2CLC \label{sec:COCCLC}
!
! !DESCRIPTION:
!  TODO - check this.
!
! ROUTINE TO CARRY OUT CO2 CALCULATIONS WITH 2 FIXED PARAMETERS ACCORDI
! THE EQUATIONS GIVEN BY PARKS(1969) AND SKIRROW (1975)
! WITH ADDITIONS FOR INCLUDING BORON IF BORON=.TRUE.
!\\
!\\
! !INTERFACE:
      SUBROUTINE CO2CLC(k0co2,k1co2,k2co2,kb,ICALC,BORON,BTOT,ctot,TA,pH,PCO2,H2CO3,HCO3,CO3,success)
!
! !USES:
!
! !LOCAL VARIABLES:
      real(rk),intent(in)    :: k0co2,k1co2,k2co2,kb
      real(rk),intent(inout) :: ctot,TA,pH,PCO2,H2CO3,HCO3,CO3
      logical,intent(out)    :: success
      INTEGER ICALC, KARL, LQ
!put counter in to check duration in convergence loop
      INTEGER                :: COUNTER,C_CHECK,C_SW
      real(rk)              :: ALKC, ALKB,BTOT
      real(rk)              :: AKR,AHPLUS
      real(rk)              :: PROD,tol1,tol2,tol3,tol4,steg,fak
      real(rk)              :: STEGBY,Y,X,W,X1,Y1,X2,Y2,FACTOR,TERM,Z
      real(rk)              :: CO2sys_memory(7)
      LOGICAL               :: BORON,DONE
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
!  DERIVING PH REQUIRES FOLLOWING LOOP TO CONVERGE.
!  THIS SUBROUTINE RELIES ON CONVERGENCE.  IF THE ENVIRONMENTAL
!  CONDITIONS DO NOT ALLOW FOR CONVERGENCE (IN 3D MODEL THIS IS
!  LIKELY TO OCCUR NEAR LOW SALINITY REGIONS) THE MODEL WILL 
!  BE STUCK IN THE LOOP.  TO AVOID THIS A CONVERGENCE CONDITION
!  IS PUT IN PLACE TO SET A FLAGG OF -99 IN THE PH VAR FOR NON CONVEGENCE. 
!  THE MODEL IS THEN ALLOWED TO CONTINUE. 'COUNTER, C_SW,C_CHECK' ARE 
!  THE LOCAL VARS USED.
! C_SW = condition of convergence 0=yes, 1= no
! COUNTER = number of iterations
! C_CHECK = maximum number of iterations
      success = .true.
      DONE = .false.
! SET COUNTER AND SWITCH TO ZERO AND OFF
      COUNTER=0
      C_SW=0
! FROM EXPERIENCE IF THE ITERATIONS IN THE FOLLOWING DO LOOP
! EXCEEDS 15 CONVERGENCE WILL NOT OCCUR.  THE OVERHEAD OF 25 ITERATIONS 
! IS OK FOR SMALL DOMAINS WITH 1/10 AND 1/15 DEG RESOLUTION.
! I RECOMMEND A LOWER VALUE OF 15 FOR HIGHER RESOLUTION OR LARGER DOMAINS.
      C_CHECK=25
! IF CONVERGENCE IS NOT ACHIEVED THE LOCAL ARRAY CONCS MUST BE STORED TO 
! ALLOW THE MODEL TO CONTINUE. THEREFORE ....
       CO2sys_memory(1) = Ctot
       CO2sys_memory(2) = TA
       CO2sys_memory(3) = pH
       CO2sys_memory(4) = PCO2
       CO2sys_memory(5) = H2CO3
       CO2sys_memory(6) = HCO3
       CO2sys_memory(7) = CO3

!      DO II=1,NKVAL
 !       AKVAL2(II)=AKVAL(II)
  !    END DO

      AKR = k1co2/k2co2
      AHPLUS=10._rk**(-PH)
      PROD=AKR*k0co2*PCO2

      IF(BORON) THEN

        IF(ICALC.EQ.1.OR.ICALC.EQ.4) THEN
!         *** TA, BTOT AND CTOT OR PCO2 FIXED ***
!         *** ITERATIVE CALCULATION NECESSARY HERE

!         SET INITIAL GUESSES AND TOLERANCE

          H2CO3=PCO2*k0co2
          CO3=TA/10._rk
          AHPLUS=1.e-8_rk
          ALKB=BTOT
          TOL1=TA/1.e5_rk
          TOL2=H2CO3/1.e5_rk
          TOL3=CTOT/1.e5_rk
          TOL4=BTOT/1.e5_rk
!         HALTAFALL iteration to determine CO3, ALKB, AHPLUS
          KARL=1
          STEG=2._rk
          FAK=1._rk
          STEGBY=0.4_rk

          DO WHILE (.not.DONE)
            DONE=.TRUE.

!
! SET COUNTER UPDATE. 
            COUNTER=COUNTER+1

! CHECK IF CONVERGENCE HAS OCCURED IN THE NUMBER OF 
! ACCEPTABLE ITTERATIONS.  
            if(counter.ge.c_check)then
!!        IF(MASTER)THEN
!!! LOG FILE TO SHOW WHEN AND WHERE NON CONVERGENCE OCCURS.
!!           PPWRITELOG 'ERRORLOG ',III,' ',(CONCS2(II),II=1,NCONC)
!!        ENDIF
! IF NON CONVERGENCE, THE MODEL REQUIRES CONCS TO CONTAIN USABLE VALUES.
! BEST OFFER BEING THE OLD CONCS VALUES WHEN CONVERGENCE HAS BEEN 
! ACHIEVED
              Ctot = CO2sys_memory(1)
              TA   = CO2sys_memory(2)
              pH   = CO2sys_memory(3)
              PCO2 = CO2sys_memory(4)
              H2CO3= CO2sys_memory(5)
              HCO3 = CO2sys_memory(6)
              CO3  = CO2sys_memory(7)
              success = .false.

!RESET SWITCH FOR NEXT CALL TO THIS SUBROUTINE
              C_SW=0
              RETURN
!
            endif

            IF(ICALC.EQ.4) THEN
!         *** PCO2 IS FIXED ***
              Y=AHPLUS*AHPLUS*CO3/(k1co2*k2co2)
              IF(ABS(Y-H2CO3).GT.TOL2) THEN
                CO3=CO3*H2CO3/Y
                DONE=.FALSE.
              ENDIF
            ELSEIF(ICALC.EQ.1) THEN
!           *** CTOT IS FIXED ***
              Y=CO3*(1._rk+AHPLUS/k2co2+AHPLUS*AHPLUS/(k1co2*k2co2))
              IF(ABS(Y-CTOT).GT.TOL3) THEN
                CO3=CO3*CTOT/Y
                DONE=.FALSE.
              ENDIF
            ENDIF
            Y=ALKB*(1._rk+AHPLUS/kb)
            IF(ABS(Y-BTOT).GT.TOL4) THEN
              ALKB=ALKB*BTOT/Y
              DONE=.FALSE.
            ENDIF

! Alkalinity is equivalent to -(total H+), so the sign of W is opposite
! to that normally used

            Y=CO3*(2._rk+AHPLUS/k2co2)+ALKB
            IF(ABS(Y-TA).GT.TOL1) THEN
              DONE=.FALSE.
              X=LOG(AHPLUS)
              W=SIGN(1._rk,Y-TA)
              IF(W.GE.0._rk) THEN
                X1=X
                Y1=Y
              ELSE
                X2=X
                Y2=Y
              ENDIF
              LQ=KARL
              IF(LQ.EQ.1) THEN
                KARL=2*NINT(W)
              ELSEIF(IABS(LQ).EQ.2.AND.(LQ*W).LT.0._rk) THEN
                FAK=0.5_rk
                KARL=3
              ENDIF
              IF(KARL.EQ.3.AND.STEG.LT.STEGBY) THEN
                W=(X2-X1)/(Y2-Y1)
                X=X1+W*(TA-Y1)
              ELSE
                STEG=STEG*FAK
                X=X+STEG*W
              ENDIF
              AHPLUS=EXP(X)
            ENDIF
! LOOP BACK UNTIL CONVERGENCE HAS BEEN ACHIEVED 
! OR MAX NUMBER OF ITERATIONS (C_CHECK) HAS BEEN REACHED.
          ENDDO

          HCO3=CO3*AHPLUS/k2co2
          IF(ICALC.EQ.4) THEN
            CTOT=H2CO3+HCO3+CO3
          ELSEIF(ICALC.EQ.1) THEN
            H2CO3=HCO3*AHPLUS/k1co2
            PCO2=H2CO3/k0co2
          ENDIF
          PH=-LOG10(AHPLUS)
          ALKC=TA-ALKB
        ELSEIF(ICALC.EQ.2) THEN
!         *** CTOT, PCO2, AND BTOT FIXED ***
          Y=SQRT(PROD*(PROD-4._rk*k0co2*PCO2+4._rk*CTOT))
          H2CO3=PCO2*k0co2
          HCO3=(Y-PROD)/2._rk
          CO3=CTOT-H2CO3-HCO3
          ALKC=HCO3+2._rk*CO3
          AHPLUS=k1co2*H2CO3/HCO3
          PH=-LOG10(AHPLUS)
          ALKB=BTOT/(1._rk+AHPLUS/kb)
          TA=ALKC+ALKB
        ELSEIF(ICALC.EQ.3) THEN
!         *** CTOT, PH AND BTOT FIXED ***
          FACTOR=CTOT/(AHPLUS*AHPLUS+k1co2*AHPLUS+k1co2*k2co2)
          CO3=FACTOR*k1co2*k2co2
          HCO3=FACTOR*k1co2*AHPLUS
          H2CO3=FACTOR*AHPLUS*AHPLUS
          PCO2=H2CO3/k0co2
          ALKC=HCO3+2._rk*CO3
          ALKB=BTOT/(1._rk+AHPLUS/kb)
          TA=ALKC+ALKB
        ELSEIF(ICALC.EQ.5) THEN
!         *** TA, PH AND BTOT FIXED ***
          ALKB=BTOT/(1._rk+AHPLUS/kb)
          ALKC=TA-ALKB
          HCO3=ALKC/(1._rk+2._rk*k2co2/AHPLUS)
          CO3=HCO3*k2co2/AHPLUS
          H2CO3=HCO3*AHPLUS/k1co2
          PCO2=H2CO3/k0co2
          CTOT=H2CO3+HCO3+CO3
        ELSEIF(ICALC.EQ.6) THEN
!         *** PCO2, PH AND BTOT FIXED ***
          ALKB=BTOT/(1._rk+AHPLUS/kb)
          H2CO3=PCO2*k0co2
          HCO3=H2CO3*k1co2/AHPLUS
          CO3=HCO3*k2co2/AHPLUS
          CTOT=H2CO3+HCO3+CO3
          ALKC=HCO3+2._rk*CO3
          TA=ALKC+ALKB
        ENDIF
      ELSE
        IF(ICALC.EQ.1) THEN
!         *** CTOT AND TA FIXED ***
          TERM=4._rk*TA+CTOT*AKR-TA*AKR
          Z=SQRT(TERM*TERM+4._rk*(AKR-4._rk)*TA*TA)
          CO3=(TA*AKR-CTOT*AKR-4._rk*TA+Z)/(2._rk*(AKR-4._rk))
          HCO3=(CTOT*AKR-Z)/(AKR-4._rk)
          H2CO3=CTOT-TA+CO3
          PCO2=H2CO3/k0co2
          PH=-LOG10(k1co2*H2CO3/HCO3)
        ELSEIF(ICALC.EQ.2) THEN
!         *** CTOT AND PCO2 FIXED ***
          Y=SQRT(PROD*(PROD-4._rk*k0co2*PCO2+4._rk*CTOT))
          H2CO3=PCO2*k0co2
          HCO3=(Y-PROD)/2._rk
          CO3=CTOT-H2CO3-HCO3
          TA=HCO3+2._rk*CO3
          PH=-LOG10(k1co2*H2CO3/HCO3)
        ELSEIF(ICALC.EQ.3) THEN
!         *** CTOT AND PH FIXED ***
          FACTOR=CTOT/(AHPLUS*AHPLUS+k1co2*AHPLUS+k1co2*k2co2)
          CO3=FACTOR*k1co2*k2co2
          HCO3=FACTOR*k1co2*AHPLUS
          H2CO3=FACTOR*AHPLUS*AHPLUS
          PCO2=H2CO3/k0co2
          TA=HCO3+2._rk*CO3
        ELSEIF(ICALC.EQ.4) THEN
!         *** TA AND PCO2 FIXED ***
          TERM=SQRT((8._rk*TA+PROD)*PROD)
          CO3=TA/2._rk+PROD/8._rk-TERM/8._rk
          HCO3=-PROD/4._rk+TERM/4._rk
          H2CO3=PCO2*k0co2
          CTOT=CO3+HCO3+H2CO3
          PH=-LOG10(k1co2*H2CO3/HCO3)
        ELSEIF(ICALC.EQ.5) THEN
!         *** TA AND PH FIXED ***
          HCO3=TA/(1._rk+2._rk*k2co2/AHPLUS)
          CO3=HCO3*k2co2/AHPLUS
          H2CO3=HCO3*AHPLUS/k1co2
          PCO2=H2CO3/k0co2
          CTOT=H2CO3+HCO3+CO3
        ELSEIF(ICALC.EQ.6) THEN
!         *** PCO2 AND PH FIXED ***
          H2CO3=PCO2*k0co2
          HCO3=H2CO3*k1co2/AHPLUS
          CO3=HCO3*k2co2/AHPLUS
          CTOT=H2CO3+HCO3+CO3
          TA=HCO3+2._rk*CO3
        ENDIF
      ENDIF

      END SUBROUTINE CO2CLC
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CaCO3_Saturation \label{sec:CaCO3Saturation}
!
! !DESCRIPTION:
!  TODO - check this.
!
! Routine to calculate the saturation state of calcite and aragonite    
! Inputs:                                                               
!               Tc      Temperature (C)                                 
!               S       Salinity                                        
!               Pr      Pressure (Pa)                                       
!               CO3     Carbonate ion concentration (mol.-l ie /1D6)    
!                                                                       
! Outputs                                                               
!               Om\_cal  Calite saturation                               
!               Om\_arg  Aragonite saturation                            
!                                                                       
! Intermediates                                                         
!               K\_cal   Stoichiometric solubility product for calcite   
!               K\_arg   Stoichiometric solubility product for aragonite 
!               Ca      Calcium 2+ concentration (mol.kg-1) 
!               P       Pressure (bars)
!                                                                       
! Source                                                                
!       Zeebe and Wolf-Gladrow 2001 following Mucci (1983)                
!       with pressure corrections from Millero (1995)                   
!       Code tested against reference values given in Z and W-G   
!\\
!\\
! !INTERFACE:
      SUBROUTINE CaCO3_Saturation (Tc, S, Pr, CO3, Om_cal, Om_arg)        
         real(rk),intent(in)  :: Tc, S, Pr, CO3
         real(rk),intent(out) :: Om_cal, Om_arg
!
! !LOCAL VARIABLES:
!       ! TODO - these are ...
        real(rk) Tk, Ca
        real(rk) logKspc, Kspc
        real(rk) logKspa, Kspa
        real(rk) tmp1, tmp2, tmp3
        real(rk) dV, dK, P, R

        real(rk),parameter :: Kelvin = 273.15_rk
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
! setup
        Tk = Tc + Kelvin 
        Ca = 0.01028_rk    ! Currently oceanic mean value at S=25, needs refining)
        R = 83.131_rk      !(cm3.bar.mol-1.K-1)
        P = Pr*1.e-5_rk    !pressure in bars

! calculate K for calcite
        tmp1 = -171.9065_rk - (0.077993_rk*Tk) + (2839.319_rk/Tk) + 71.595_rk*log10(Tk) 
        tmp2 = + (-0.77712_rk + (0.0028426_rk*Tk) + (178.34_rk/Tk))*SQRT(S) 
        tmp3 = - (0.07711_rk*S) + (0.0041249_rk*(S**1.5_rk))
        logKspc = tmp1 + tmp2 + tmp3
        Kspc = 10._rk**logKspc
      
! correction for pressure for calcite
        dV = -48.76_rk + 0.5304_rk*Tc
        dK = -11.76_rk/1.e3_rk + (0.3692_rk/1.e3_rk) * Tc
        tmp1 = -(dV/(R*Tk))*P + (0.5_rk*dK/(R*Tk))*P*P
        Kspc = Kspc*exp(tmp1)
        logKspc = log10(Kspc)
        
! calculate K for aragonite
        tmp1 = -171.945_rk - 0.077993_rk*Tk + 2903.293_rk / Tk + 71.595_rk* log10(Tk)
        tmp2 = + (-0.068393_rk + 0.0017276_rk*Tk + 88.135_rk/Tk)*SQRT(S)    
        tmp3 = - 0.10018_rk*S + 0.0059415_rk*S**1.5_rk
        logKspa = tmp1 + tmp2 + tmp3
        Kspa = 10._rk**logKspa
     
! correction for pressure for aragonite
        dV = -46._rk + 0.530_rk*Tc
        dK = -11.76_rk/1.e3_rk + (0.3692_rk/1.e3_rk) * Tc
        tmp1 = -(dV/(R*Tk))*P + (0.5_rk*dK/(R*Tk))*P*P
        Kspa = Kspa*exp(tmp1)
        logKspa = log10(Kspa)
        
! calculate saturation states
        Om_cal = (CO3 * Ca) / Kspc
        Om_arg = (CO3 * Ca) / Kspa

      END SUBROUTINE CaCO3_SATURATION
!
!EOC
   end module
