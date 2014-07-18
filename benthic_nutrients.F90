#include "fabm_driver.h"

module pml_ersem_benthic_nutrients

   use fabm_types

   use pml_ersem_shared

   implicit none

   private

   type,extends(type_base_model),public :: type_pml_ersem_benthic_nutrients
      type (type_bottom_state_variable_id) :: id_K3n,id_K4n,id_G2o,id_D1m,id_D2m
      type (type_state_variable_id)        :: id_N3n,id_N4n,id_O2o
      type (type_dependency_id) :: id_ETW,id_phx
      
      real(rk) :: q10nitX,hM4M3X,sM4M3X,xno3X
      real(rk) :: EDZ_1X,EDZ_2X,EDZ_3X,d_totX
      real(rk) :: qPWX,M1adsX,M4adsX,M11adsX,relax_oX,relax_mX,EDZ_mixX
      integer :: ISWphx
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type
   
contains   
   
   subroutine initialize(self,configunit)
      class (type_pml_ersem_benthic_nutrients),intent(inout),target :: self
      integer,                                 intent(in)           :: configunit
      
      ! Set time unit to d-1. This implies that all rates (sink/source terms) are given in d-1.
      self%dt = 86400._rk

      call self%get_parameter(self%q10nitX,'q10nitX','-','Q_10 temperature coefficient')
      call self%get_parameter(self%hM4M3X,'hM4M3X','mmol/m^3','Michaelis-Menten constant for nitrate limitation of nitrification')
      call self%get_parameter(self%ISWphx,'ISWphx','','pH influence on nitrification',default=0)
      call self%get_parameter(self%sM4M3X,'sM4M3X','1/d','maximum nitrification rate')
      call self%get_parameter(self%xno3X,'xno3X','mol O_2/mol N','oxygen consumed per nitrate produced in nitrification')
      call self%get_parameter(self%EDZ_1X,'EDZ_1X','m^2/d','diffusivity in 1st (oxygenated) layer')
      call self%get_parameter(self%EDZ_2X,'EDZ_2X','m^2/d','diffusivity in 2nd (oxidized) layer')
      call self%get_parameter(self%EDZ_3X,'EDZ_3X','m^2/d','diffusivity in 3rd (anoxic) layer')
      call self%get_parameter(self%qPWX,'qPWX','-','fraction of pore water in the sediment')
      call self%get_parameter(self%M1adsX,'M1adsX','-','fraction of phosphorus adsorbed to particulates in 1st layer')
      call self%get_parameter(self%M4adsX,'M4adsX','-','fraction of ammonium adsorbed to particulates')
      call self%get_parameter(self%M11adsX,'M11adsX','-','fraction of phosphorus adsorbed to particulates in 2nd layer')
      call self%get_parameter(self%relax_oX,'relax_oX','1/d','oxygen diffusion time scale for benthic/pelagic interface')
      call self%get_parameter(self%relax_mX,'relax_mX','1/d','nitrate diffusion time scale for benthic/pelagic interface')
      call self%get_parameter(self%EDZ_mixX,'EDZ_mixX','d/m','equilibrium diffusive speed between sediment surface water')
      call self%get_parameter(self%d_totX,'d_totX','m','depth of sediment column')
      
      call self%register_state_variable(self%id_K3n,'K3n','mmol/m^2','benthic nitrate')
      call self%register_state_variable(self%id_K4n,'K4n','mmol/m^2','benthic ammonium')
      call self%register_state_variable(self%id_G2o,'G2o','mmol O_2/m^2','benthic oxygen')
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_K3n)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_K4n)

      call self%register_state_variable(self%id_D1m,'D1m','m','depth of bottom interface of 1st layer')
      call self%register_state_variable(self%id_D2m,'D2m','m','depth of bottom interface of 2nd layer')

      call self%register_state_dependency(self%id_N3n,'N3n','mmol N/m^3','nitrate')
      call self%register_state_dependency(self%id_N4n,'N4n','mmol N/m^3','ammonium')
      call self%register_state_dependency(self%id_O2o,'O2o','mmol O_2/m^3','oxygen')
      
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      call self%register_dependency(self%id_phx,standard_variables%ph_reported_on_total_scale)
   end subroutine

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !INPUT PARAMETERS:
   class (type_pml_ersem_benthic_nutrients),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
   
      real(rk) :: K3n,K3nP,K4nP,G2oP
      real(rk) :: N4n,N3nP,N4nP,O2oP
      real(rk) :: ETW,phx
      real(rk) :: D1m,D2m
      real(rk) :: irrenh
      real(rk) :: MU_m,eT,eN,Fph,jM4M3n
      real(rk) :: diff1,diff2,diff3,cmix
      real(rk) :: d1,d2,d3
      real(rk) :: poro,M1ads
      real(rk) :: vG2,vM3,vM4,vM1,vM11,vM5,vG3
      real(rk) :: jmun
      real(rk) :: jmu(7),jmi(7),jmn(7)
      real(rk) :: profO2(15),profNO3(15),profN(15),profP(15)
      real(rk) :: profS(15),profCO2(15),profD(15)
      real(rk) :: D1m0,D2m0
      real(rk) :: dum
      real(rk) :: H1_eq,H2_eq
      real(rk) :: c_bot_eq, c_int1_eq, c_int2_eq, c_int_eq
   
      _HORIZONTAL_LOOP_BEGIN_

      _GET_HORIZONTAL_(self%id_K3n,K3n) !TODO background
      _GET_HORIZONTAL_(self%id_K3n,K3nP)
      _GET_HORIZONTAL_(self%id_K4n,K4nP)
      _GET_HORIZONTAL_(self%id_G2o,G2oP)

      _GET_HORIZONTAL_(self%id_D1m,D1m)
      _GET_HORIZONTAL_(self%id_D2m,D2m)

      _GET_(self%id_N3n,N3nP)
      _GET_(self%id_N4n,N4nP)
      _GET_WITH_BACKGROUND_(self%id_N4n,N4n)
      _GET_(self%id_O2o,O2oP)

      _GET_(self%id_ETW,ETW)
      _GET_(self%id_phx,phx)
      
      jmu = 0.0_rk
      jmi = 0.0_rk
      jmn = 0.0_rk
      
      irrenh = 1._rk !TODO: retrieve from benthic fauna
      
      D1m0 = 0.0_rk
      D2m0 = 0.0_rk

!* Nitrification:
! Reaction: (NH4+,OH-) + 2*O2 -> (NO3-,H+) + 2*H2O
! treated as first order reaction for NH4 in the oxygenated layer
! Average nitrate density (mmol/m3): MU_m
        
      MU_m = max(0.0_rk,K3n)/(D1m+(D2m-D1m)/3.0_rk)

      eT = EXP(LOG(self%q10nitX)*(ETW-10._rk)/10._rk)
      eN = self%hM4M3X/(self%hM4M3X+MU_m)
!Ph influence on nitrification - empirical equation
      if(self%ISWphx.eq.1)then
        Fph = MIN(2._rk,MAX(0._rk,0.6111_rk*phx-3.8889_rk))
        jM4M3n = Fph * self%sM4M3X * K4nP * (D1m/self%d_totX) * eT * eN
      else
        jM4M3n = self%sM4M3X*K4nP*(D1m/self%d_totX)*eT*eN
      endif

      _SET_ODE_BEN_(self%id_K4n, -jM4M3n)
      _SET_ODE_BEN_(self%id_K3n,  jM4M3n)
      _SET_ODE_BEN_(self%id_G2o, -self%xno3X*jM4M3n)

      jMU(1) = jMU(1) - jM4M3n
      jMU(6) = jMU(6) + jM4M3n
      jMU(4) = jMU(4) - self%xno3X*jM4M3n

!!* Denitrification:
!! The part pammonX (maximal) of anaerobic respiration reduces  NO3.
!! From this NO3 reduction the part pdenitX goes to N2-gas.
!
!      MU_m = K3nP(K)/(D1m(K)+(D2m(K)-D1m(K))/3._rk)
!      eN = eMM(MU_m/3._rk,hM3G4X)
!! "borrowed" oxygen consumption in anaerobic layer:
!
!! to avoid numerical problems, use the recalculated value from the record
!        jMIo2=jMI(4)
!
!! "expression in nitrate reduction (100%):
!      jMIno3 = -jMIo2/(xno3X - xn2X*pdenitX)
!      jM3M4n = jMIno3*eN*pammonX*(1._rk-pdenitX)
!!      jjM3M4n(i) = jMIno3*eN*pammonX*(1.0d0-pdenitX)
!      jM3G4n = jMIno3*eN*pammonX*     pdenitX
!
!      SK14n = SK14n +       jM3M4n
!      SK3n(K)  = SK3n(K)  -       jM3M4n -      jM3G4n
!      SG4n(K)  = SG4n(K)  +       jM3G4n
!      SK6e  = SK6e  + xno3X*jM3M4n + (xno3X-xn2X)*jM3G4n
!
!      jMI(1)   = jMI(1)   +       jM3M4n
!      jMI(6) = jMI(6) -       jM3M4n -      jM3G4n
!      jMI(7)  = jMI(7)  +       jM3G4n
!      jMI(4)  = jMI(4)  + xno3X*jM3M4n + (xno3X-xn2X)*jM3G4n

!!--------------------------------------------
!! Silicateregeneration - preliminary
!!--------------------------------------------
!
!!aerobic layer: 0 .. D1m(K)
!      rsil = sQ6M5X*partq(D9m(K),0._rk,D1m(K),d_totX)
!
!      SQ6s(K) = SQ6s(K) - rsil*Q6sP(K)
!      SK5s(K) = SK5s(K) + rsil*Q6sP(K)
!
!      jMU(3) = jMU(3) + rsil*Q6sP(K)
!
!      SD9m(K) = SD9m(K) + (D1m(K)/2._rk-D9m(K))*rsil
!#ifdef ERSEMDEBUG
!      if(ersem_debugger.gt.0 .and. k .eq. kp_dbg) then
!       PPWRITEDBGALLPRCS 'sil regen:',SD9m(K),D1m(K),D9m(K),rsil
!      endif
!#endif
!
!!lower layer: D1m(K) .. d_totX
!      rsil = sQ6M5X*partq(D9m(K),D1m(K),d_totX,d_totX)
!
!      SQ6s(K) = SQ6s(K) - rsil*Q6sP(K)
!      SK5s(K) = SK5s(K) + rsil*Q6sP(K)
!
!      jMI(3) = jMI(3) + rsil*Q6sP(K)
!
!      SD9m(K) = SD9m(K) + D1m(K)*rsil
!------------------------------------------
! Nutrient profiles and surface gradients
!------------------------------------------

! Preparing the profile:
      diff1  = self%EDZ_1X *irrenh
      diff2  = self%EDZ_2X *irrenh
      diff3  = self%EDZ_3X *irrenh
      cmix   = 0._rk
      d1  = D1m
      d2  = D2m
      d3  = self%d_totX

! Volume factors:
      poro=self%qPWX
      M1ads=self%M1adsX

      !IF (N_COMP-N_UPPERX.GT.0) THEN
      !   poro=benthic_morfology(1,K) ! porosity factor
      !   M1ads=benthic_morfology(2,K) ! adsorption factor
      !ENDIF

      vG2  = poro
      vM3  = poro
      vM4  = poro*self%M4adsX
      vM1  = poro*M1ads
      vM11 = poro*self%M11adsX
      vM5  = poro
      vG3  = poro

! profile parameters:

      jmun = jmu(1)
      IF (jmun<0._rk) THEN
         jmu(1) = jmun*N4n/(N4n+0.5_rk)
         jmi(1) = jmi(1) + (jmun - jmu(1))
      ENDIF
      jmu(5) = 0.0_rk !TODO: SG3c(K)
      jmi(5) = 0._rk

      CALL Prof_Parameter(profO2,  jMU(4), jMI(4), 0._rk, vG2,vG2,vG2)
      CALL Prof_Parameter(profNO3, jMU(6),jMI(6),0._rk, vM3,vM3,vM3)
      CALL Prof_Parameter(profN,   jMU(1),  jMI(1),  0._rk, vM4,vM4,vM4)
      CALL Prof_Parameter(profP,   jMU(2),  jMI(2),  0._rk, vM1,vM1,vM11)
      CALL Prof_Parameter(profS,   jMU(3),  jMI(3),  0._rk, vM5,vM5,vM5)
      CALL Prof_Parameter(profCO2, jMU(5),jMI(5),0._rk, vG3,vG3,vG3)
      CALL Prof_Parameter(profD,   d1,d2-d1,d3-d2, 1._rk,1._rk,1._rk)

      CALL EquProfile(0._rk,1._rk,profD,dum,cmix,d1,d2,d3,diff1,diff2,diff3)

      cmix   = self%EDZ_mixX

! oxygen

      ! -----------------------------------------------------------------------------------
      ! Oxygen
      ! -----------------------------------------------------------------------------------
      ! Layer 1: compute steady-state layer height and layer integral
      call compute_parabola_end(diff1,modconc(O2oP,jMU(4),cmix),jMU(4),self%d_totX,H1_eq,c_int1_eq)

      !CALL EndProfile(O2oP,G2oP,profO2,jMN(4),cmix,d1,d2,d3,diff1,diff2,diff3)

      ! Relax benthic oxygen towards equilibrium value
      c_int_eq = c_int1_eq*poro
      jMN(4) = jMN(4) + (c_int_eq-G2oP)/self%relax_oX

      ! Relax depth of first/oxygenated layer towards equilibrium value (H1_eq)
      _SET_ODE_BEN_(self%id_D1m,(max(D1m0,H1_eq)-D1m)/self%relax_oX)

      ! -----------------------------------------------------------------------------------
      ! Nitrate
      ! -----------------------------------------------------------------------------------
      ! Layer 1: compute steady-state concentration at bottom interface and layer integral
      call compute_parabola(d1,diff1,modconc(N3nP,jMU(6)+jMI(6),cmix),jMU(6),jMI(6),c_bot_eq,c_int1_eq)
      ! Layer 2: compute steady-state layer height and layer integral
      call compute_parabola_end(diff2,c_bot_eq,jMI(4),self%d_totX-d1,H2_eq,c_int2_eq)

! nitrate

      !CALL EndProfile(N3nP,K3nP,profNO3,jMN(6),cmix,d1,d2,d3,diff1,diff2,diff3)

      ! Relax benthic nitrate towards equilibrium value
      c_int_eq = (c_int1_eq+c_int2_eq)*poro
      jMN(6) = jMN(6) + (c_int_eq-K3nP)/self%relax_mX

      ! Relax depth of bottom interface of second/oxidised layer towards equilibrium value (H1_eq+H2_eq)
      _SET_ODE_BEN_(self%id_D2m,(max(D2m0,H1_eq+H2_eq)-D2m)/self%relax_mX)

! ammonium

      CALL EquProfile(N4nP,K4nP,profN,jMN(1),cmix,d1,d2,d3,diff1,diff2,diff3)
      ! Non-equilibrium correction:
      CALL NonEquFlux(profN,profD,jMN(1))
      jMN(1) = profN(1) + profN(2) + profN(3)

!! phospate
!
!      CALL EquProfile(N1pP(I),K1pP(K),profP,jMN(2),k)
!      ! Non-equilibrium correction:
!      CALL NonEquFlux(profP,profD,jMN(2))
!      jMN(2) = profP(1) + profP(2) + profP(3)
!
!#ifdef ERSEMDEBUG
!      if(ersem_debugger.gt.0 .and. k .eq. kp_dbg) then
!        PPWRITEDBGALLPRCS 'jmn(`2)',jmn(2)
!        PPWRITEDBGALLPRCS profP(1)
!        PPWRITEDBGALLPRCS profP(2)
!        PPWRITEDBGALLPRCS profP(3)
!        PPWRITEDBGALLPRCS N1pP(I)
!        PPWRITEDBGALLPRCS K1pP(K)
!      endif
!#endif
!
!! silicate
!      
!      CALL EquProfile(N5sP(I),K5sP(K),profS,jMN(3),k)
!      ! Non-equilibrium correction:
!      CALL NonEquFlux(profS,profD,jMN(3))
!      jMN(3) = profS(1) + profS(2) + profS(3)
!
!! carbon dioxide
!     
!      CALL EquProfile(O3cP(I),G3cP(K),profCO2,jMN(5),k)
!      ! Non-equilibrium correction:
!      CALL NonEquFlux(profCO2,profD,jMN(5))
!      jMN(5) = profCO2(1) + profCO2(2) + profCO2(3)
!
!! nitrate gas:
!      
!      jmn(7) = 0._rk
!
!!-----------
!! Diffusion
!!-----------
!
      _SET_ODE_BEN_(self%id_K4n, - jMN(1))
      _SET_BOTTOM_EXCHANGE_(self%id_N4n, + jMN(1))
      !_SET_ODE_BEN_(self%id_K1p, - jMN(2))
      !_SET_BOTTOM_EXCHANGE_(self%id_N1p, + jMN(2))
      !_SET_ODE_BEN_(self%id_K5s, - jMN(3))
      !_SET_BOTTOM_EXCHANGE_(self%id_N5s, + jMN(3))
      _SET_ODE_BEN_(self%id_G2o, - jMN(4))
      _SET_BOTTOM_EXCHANGE_(self%id_O2o, + jMN(4))
      !_SET_ODE_BEN_(self%id_G3c, - jMN(5))
      !_SET_BOTTOM_EXCHANGE_(self%id_O3c, + jMN(5))
      _SET_ODE_BEN_(self%id_K3n, - jMN(6))
      _SET_BOTTOM_EXCHANGE_(self%id_N3n,jMN(6))
      !_SET_ODE_BEN_(self%id_G4n, - jMN(7))

      _HORIZONTAL_LOOP_END_

   end subroutine

   subroutine compute_parabola(H,D,y0,sms_int,flux_bot,y_bot,y_int)
      real(rk),intent(in)  :: H,D,y0,sms_int,flux_bot
      real(rk),intent(out) :: y_bot,y_int
      real(rk) :: a,b,c
      ! ----------------------------------------------------------------------------------------------------
      ! Determine equilibrium distribution from:
      ! - layer height (m)
      ! - diffusivity (m2/d)
      ! - top concentration (#/m3)
      ! - layer-integrated source terms (#/m2/d)
      ! - bottom flux (#/m2/d)
      ! Returns equilibrium concentration at bottom interface (#/m3), depth integral of layer contents (#/m2)
      ! ----------------------------------------------------------------------------------------------------
      ! Governing equation: dy/dt = D d2y/dz2 + sms
      ! Assumption: diffusivity D and sources-minus-sinks sms are independent of z within the layer.
      ! Thus, sms can be written as layer integral divided by layer height: sms = sms_int/H
      ! At equilibrium: D d^2y/dz^2 + sms_int/H = 0
      ! Thus, d^2y/dz^2 = -sms_int/H/D
      ! Solution is a quadratic equation: y(z) = a z^2 + b z + c
      ! From second derivative: d^2y/dz^2 = 2a = -sms_int/H/D. Thus, a = -sms_int/H/D/2.
      ! ----------------------------------------------------------------------------------------------------
      ! y(z) = -sms_int/H/D/2 z^2 + b z + c
      !
      ! Lets adopt a bottom-to-top coordinate system, with the bottom interface at z=0 and the top interface at z=H.
      !
      ! Constraint 1: flux over bottom interface is known: flux_bot
      ! (typically chosen to balance demand depth-integrated sinks-sources in deeper layers)
      ! For consistency, this flux must equal that produced by diffusion at the boundary, i.e., D*dy/dz
      ! Note: gradient (bottom to top!) must be positive when deeper layers are a sink (i.e., sms_int_deep<0)
      ! D*dy/dz(0) = D*b = -flux_bot -> b = -flux_bot/D
      !
      ! Constraint 2: top concentration y(H)=y0 is known.
      ! y(H) = -sms_int/D/2 H - flux_bot/D H + c = y0
      ! -> c = y0 + (sms_int/2 + flux_bot)/D H
      ! This is also the concentration at bottom interface: y(0) = c
      !
      ! Depth-integrated layer contents is found by integrating the parabola between 0 and H:
      !    \int{a * z^2 + b z + c} = [a/3 z^3 + b/2 z^2 + c z]_0^H = a/3 H^3 + b/2 H^2 + c Z
      ! ----------------------------------------------------------------------------------------------------
      a = -sms_int/H/D/2
      b = -flux_bot/D
      c = y0 + H*(sms_int/2 + flux_bot)/D

      y_bot = c
      y_int = (a/3*H*H + b/2*H + c)*H
   end subroutine

   subroutine compute_parabola_end(D,y0,sms_int,Hmax,H,y_int)
      real(rk),intent(in)  :: D,y0,sms_int,Hmax
      real(rk),intent(out) :: H,y_int
      real(rk) :: y_bot
      ! ----------------------------------------------------------------------------------------------------
      ! Determine steady-state layer depth and equilibrium distribution from:
      ! - diffusivity (m2/d)
      ! - top concentration (#/m3)
      ! - layer-integrated source terms (#/m2/d)
      ! - concentration and flux at bottom interface must be zero
      ! ----------------------------------------------------------------------------------------------------
      ! Governing equation in 1st layer: dy/dt = D d2y/dz2 + sms
      ! Assumption: diffusivity D and sink-minus-source sms are independent of z within the layer.
      ! Thus, sms can be written as layer integral divided by layer height: sms = sms_int/H
      ! At equilibrium: D d^2y/dz^2 + sms_int/H = 0
      ! Thus, d^2y/dz^2 = -sms_int/H/D
      ! Solution is a quadratic equation: y(z) = a z^2 + b z + c
      ! From second derivative: d^2y/dz^2 = 2a = -sms_int/H/D. Thus, a = -sms_int/H/D/2.
      ! ----------------------------------------------------------------------------------------------------
      ! Task 1: find layer height
      !   Bottom-to-top coordinate system, bottom interface at z=0, top interface at z=H
      !   Bottom constraint: concentration y(0) = 0, flux dy/dz = 0 (i.e., minimum of parabola)
      !   Thus, b=0 and c=0, and y(z) = -sms_int/H/D/2 z^2
      !   Top constraint: concentration y(H) is known: y0
      !   Thus y(H) = -sms_int/D/2 H = y0 -> layer height H = -2 D y0/sms_int
      !   Note that H>=0 for sms_int<0 only: this solution is valid only if the layer is a sink.
      !   Also, if the layer is neither sink nor source (sms_int=0), H->infinity. That is,
      !   the solution is a flat profile, with surface concentration y0 extending forever downward.
      ! Task 2: compute depth-integrated layer contents:
      !   \int{-sms_int/H/D/2 z^2} = [-sms_int/H/D/6 z^3]_0^H = -sms_int/D/6 H^2
      ! Check consistency: flux at top interface is D*dy/dz = -sms_int
      !    [OK: in steady state, surface exchange compensates internal loss]
      ! ----------------------------------------------------------------------------------------------------
      if (Hmax*sms_int>-2*D*y0) then
         ! Loss rate within layer is too low (or layer experiences net production, i.e., sms_int>0).
         ! Layer would extend beyond maximum depth. Fix depth at maximum depth and
         ! use parabola with zero flux but non-zero concentration at bottom interface.
         H = Hmax
         call compute_parabola(H,D,y0,sms_int,0.0_rk,y_bot,y_int)
      else
         H = -2*D*y0/sms_int
         y_int = -sms_int/D/6*H*H
      end if
   end subroutine

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Prof_Parameter \label{sec:profParameter}
!
! !DESCRIPTION:
!  TODO - CHECK THIS
!
!  Load layer productions and volume factors into profile.
!\\
!\\
! !INTERFACE:
   SUBROUTINE Prof_Parameter(prof,s1,s2,s3,v1,v2,v3)
!
! !INPUT PARAMETERS:
! ! TODO - document these.
  real(rk),intent(in) :: s1,s2,s3,v1,v2,v3
!
! !OUTPUT PARAMETERS:
! ! TODO - document these
  real(rk),intent(out) :: prof(15)
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
!  code................................................

      prof(1) = s1  !layer fluxes
      prof(2) = s2
      prof(3) = s3
      prof(8) = v1 !layer volume factors
      prof(9) = v2
      prof(10) = v3

      RETURN
      END SUBROUTINE Prof_Parameter
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: EquProfile \label{sec:equProfile}
!
! !DESCRIPTION:
!  TODO - CHECK THIS.
!
!  computes diffusion-production equilibrium profiles
!    from layer productions, layer depths and benthic/pelagic interface 
!    concentration.
! 
!    a\_water: concentration in water box above the seafloor (per $m^3$)
!    mass: conent in benthic section (per $m^2$)
!    p: working array with
!       1: first layer production
!       2: second layer production
!       3: third layer production
!       4: concentration at seafloor (benthic/pelagic interface)
!       5: concentration at 1st/2nd layer interface
!       6: concentration at 2nd/3rd layer interface
!       7: concentraion at lower boundary of 3rd layer
!       8: 1st layer volume factor (porosity)
!       9: 2nd layer volume factor (porosity)
!      10: 3rd layer volume factor (porosity)
!      11: 1st layer content (per $m^2$)
!      12: 2nd layer content (per $m^2$)
!      13: 3rd layer content (per $m^2$)
!      14: computed difference in total benthic layer content
!      15: lower boundary of 3rd layer
!    f\_0: flux from benthic to pelagic layer 
!\\
!\\
! !INTERFACE:
   SUBROUTINE EquProfile(a_water,mass,p,f0,cmix,d1,d2,d3,diff1,diff2,diff3)
!
! !INPUT PARAMETERS:
   real(rk),intent(in) :: a_water,mass,cmix,d1,d2,d3,diff1,diff2,diff3
!
! !INPUT/OUTPUT PARAMETERS:
!  ! TODO - document these.
   real(rk) :: f0,p(15)
!
! !LOCAL VARIABLES:
!  ! TODO - Document these.
   real(rk) :: f1,f2
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
! fluxes from Sources:
      f0 = p(1)+ p(2)+ p(3)
      f1 =       p(2)+ p(3)
      f2 =             p(3)

! Surface concentration a0:
      p(4) = modconc(a_water,f0,cmix)

! Profile in form of three parabular pieces:
      CALL arc   (p(4),p(5), f0,f1, diff1, 0._rk,d1, p(11))
      CALL arc   (p(5),p(6), f1,f2, diff2, d1,d2, p(12))
      CALL arc   (p(6),p(7), f2,0._rk, diff3, d2,d3, p(13)) ! 0 bottom flux
      p(15) = d3

! Difference of equilibrium masses to real mass:
      p(14) = p(11)*p(8) + p(12)*p(9) + p(13)*p(10) - mass

   END SUBROUTINE EquProfile

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: modconc \label{sec:modconc}
!
! !DESCRIPTION:
!  TODO - CHECK THIS
!
!  Modification of surface concentration by the diffusion flux
!\\
!\\
! !INTERFACE:
   real(rk) FUNCTION modconc(conc,flux,cmix)
!
! !LOCAL VARIABLES:
!  ! TODO - document these.
   real(rk), intent(in) :: conc,flux,cmix
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
! conc mg/m3, flux mg/m2day, cmix day/m

      IF (flux.GE.0._rk) THEN
         modconc = conc + cmix*flux
      ELSE
         modconc = conc*conc/(conc - cmix*flux)
      ENDIF

      RETURN

   END FUNCTION modconc
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: arc \label{sec:arc}
!
! !DESCRIPTION:
!  TODO - CHECK THIS
!
!  Computes for a layer the lower face concentration and the layer content 
!  from upper face concentration and lower and upper face fluxes, assuming
!  equilibrium between inner production and diffusion.
!  Solves: \verb+a''=(f1-f0)/(diff*(d1-d0)+ 
!     with \verb+diff*a'=f0 at face 0 and diff*a'=f1\verb+ at face 1
!    \verb+=> a(d,t) = (f1-f0)/(diff*(d1-d0))*.5*(d-d0)^2+f0/diff*(d-d0)+a0+
!
!  a0: upper face concentration (per $m^3$)
!  a1: lower face concentration (per $m^3$)
!  f0: upper flace flux (outwards)
!  f1: lower face flux (inwards)
!  d0: depth level upper face
!  d1: depth level lower face
!  m1: total benthic content (per $m^2$)
!\\
!\\
! !INTERFACE:
   SUBROUTINE arc(a0,a1,f0,f1,diff,d0,d1,m1)
!
! !INPUT PARAMETERS:
! ! TODO - document these
  real(rk),intent(in) :: a0,d0,d1,f0,f1,diff
!
! !INPUT/OUTPUT PARAMETERS:
! ! TODO - document these.
  real(rk),intent(out) :: a1,m1
!
! !LOCAL VARIABLES:
! ! TODO - document these.
  real(rk) :: g0,g1
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
!  code................................................

      g0 = f0/diff
      g1 = f1/diff

      a1 = a0 + (g0+g1)*(d1-d0)/2._rk
      m1 = (a0 + (g0+g1/2._rk)*(d1-d0)/3._rk)*(d1-d0)

      RETURN
  END SUBROUTINE arc
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: endarc \label{sec:endarc}
!
! !DESCRIPTION:
!  TODO - CHECK THIS
!
!  For profiles with negative concentrations at the bottom interface cut 
!  the layer at the intersection. Works only for f1<=0.
!  a0: upper face concentration (per $m^3$)
!  a1: lower face concentration (per $m^3$)
!  f0: upper flace flux (outwards)
!  f1: lower face flux (inwards)
!  d0: depth level upper face
!  d1: depth level lower face
!  m1: total benthic content (per $m^2$)
!\\
!\\
! !INTERFACE:
  SUBROUTINE endarc(a0,a1,f0,f1,diff,d0,d1,dt,m1)
!
! !INPUT PARAMETERS:
!     ! TODO - document these.
      real(rk),intent(in) :: diff,f0,f1,a0,d0,dt
!
! !OUTPUT PARAMETERS:
!     ! TODO - document these.
      real(rk),intent(out) :: a1,d1,m1
!
! !LOCAL VARIABLES:
!     ! TODO - document these.
      real(rk) :: df,pf,vlow
      DATA vlow / 1.e-10_rk /
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOC
!
!  code................................................

         pf = f0 + f1


         if ( abs(pf).gt.vlow) then
            !df = (f0 + 2._rk*f1)/pf
            df = 1._rk + f1/pf
         else
            pf = vlow * SIGN(1._rk,pf)
            if ((abs(f0).lt.vlow).AND.(abs(f1).lt.vlow)) then
               df = 1.5_rk
            else
               if ( abs(f0).lt.vlow ) df = 2._rk
               if ( abs(f1).lt.vlow ) df = 1._rk
            endif
         endif

         d1 = d0 - 2._rk*a0/pf*diff

         IF (d1.GE.d0.AND.d1.LE.dt) THEN ! cut profile
            m1 = a0*(d1-d0)*df/3._rk
            a1 = 0._rk
         ELSE ! standard case
            CALL arc(a0,a1,f0,f1,diff,d0,dt,m1)
            d1=dt
         ENDIF

      RETURN
   END SUBROUTINE endarc
!
!EOC
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: EndProfile \label{sec:endProfile}
!
! !DESCRIPTION:
!  TODO - CHECK THIS
!  computes diffusion-production equilibrium profiles
!    from layer productions, layer depths and benthic/pelagic interface 
!    concentration. Clips negativity.
! 
!    a\_water: concentration in water box above the seafloor (per $m^3$)
!    mass: conent in benthic section (per $m^2$)
!    p: working array with
!       1: first layer production
!       2: second layer production
!       3: third layer production
!       4: concentration at seafloor (benthic/pelagic interface)
!       5: concentration at 1st/2nd layer interface
!       6: concentration at 2nd/3rd layer interface
!       7: concentraion at lower boundary of 3rd layer
!       8: 1st layer volume factor (porosity)
!       9: 2nd layer volume factor (porosity)
!      10: 3rd layer volume factor (porosity)
!      11: 1st layer content (per $m^2$)
!      12: 2nd layer content (per $m^2$)
!      13: 3rd layer content (per $m^2$)
!      14: computed difference in total benthic layer content
!      15: lower boundary of 3rd layer
!    f\_0: flux from benthic to pelagic layer
!\\
!\\
! !INTERFACE:
   SUBROUTINE EndProfile(a_water,mass,p,f0,cmix,d1,d2,d3,diff1,diff2,diff3)
!
! !INPUT PARAMETERS:
! ! TODO - document these
  real(rk),intent(in) :: a_water,mass,cmix,d1,d2,d3,diff1,diff2,diff3
!
! !INPUT/OUTPUT PARAMETERS:
! ! TODO - document this.
  real(rk) :: p(15),f0
!
! !LOCAL VARIABLES:
! ! TODO - document these
  real(rk) :: f1,f2
  LOGICAL case1,case2,case3,case4
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
! fluxes from Sources:
      f0 = p(1)+ p(2)+ p(3)
      f1 =       p(2)+ p(3)
      f2 =             p(3)

! Surface concentration a0:
       p(4) = modconc(a_water,f0,cmix)
! Case differentiation:
! case1: oxigen
! case2: nitrate
! case4: ammonium,phosphate,silicate,co2

      case1 = ((p(1).LT.0._rk) .AND. (p(2).LE.0._rk) .AND. (p(3).EQ.0._rk))
      case2 = ((p(1).GE.0._rk) .AND. (p(2).LT.0._rk) .AND. (p(3).LE.0._rk))
      case3 = ((p(1).GE.0._rk) .AND. (p(2).GE.0._rk) .AND. (p(3).LT.0._rk))
      case4 = ((p(1).GE.0._rk) .AND. (p(2).GE.0._rk) .AND. (p(3).EQ.0._rk))

! Profile in form of up to three parabular pieces:

      IF (case1) THEN
         CALL endarc(p(4),p(5), f0,f1, diff1, 0._rk,p(15),d3,p(11))
         p(12) = 0._rk
         p(13) = 0._rk
      ELSEIF (case2) THEN
         CALL arc   (p(4),p(5), f0,f1, diff1, 0._rk,d1,      p(11))
         CALL endarc(p(5),p(6), f1,f2, diff2, d1,p(15),d3, p(12))
         p(13) = 0._rk
      ELSEIF (case3) THEN
         CALL arc   (p(4),p(5), f0,f1, diff1, 0._rk,d1,      p(11))
         CALL arc   (p(5),p(6), f1,f2, diff2, d1,d2,      p(12))
         CALL endarc(p(6),p(7), f2,0._rk, diff3, d2,p(15),d3, p(13))
      ELSEIF (case4) THEN
         CALL arc   (p(4),p(5), f0,f1, diff1, 0._rk,d1,      p(11))
         CALL arc   (p(5),p(6), f1,f2, diff2, d1,d2,      p(12))
         CALL arc   (p(6),p(7), f2,0._rk, diff3, d2,d3,      p(13))
         p(15) = d3
      ELSE
         CALL arc   (p(4),p(5), f0,f1, diff1, 0._rk,d1,      p(11))
         CALL arc   (p(5),p(6), f1,f2, diff2, d1,d2,      p(12))
         CALL arc   (p(6),p(7), f2,0._rk, diff3, d2,d3,      p(13))
         p(15) = d3
      ENDIF

      p(14) = p(11)*p(8) + p(12)*p(9) + p(13)*p(10) - mass

      RETURN

      END SUBROUTINE EndProfile
!
!EOC
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: NonEquFlux \label{sec:nonEquFlux}
!
! !DESCRIPTION:
!  TODO - CHECK THIS
!
!  Adds non-equilibrium modification to equilibrium flux.  
!              assumes homogenous distribution of flux correction
!   prof:    (1-3): layer fluxes
!            (8-10): layer volume-factors (porosity+adsorption)
!            (14): Difference of eq. content to real content
!            (15): flux to pelagic layer
!   profD:   (1-3): layer thicknesses
!            (11:13): layer equilibrium contents (per $m^2$)
!   j\_diff:  outward flux on benthic surface 
!\\
!\\
! !INTERFACE:
   SUBROUTINE NonEquFlux(prof,profD,j_diff)
!
! !OUTPUT PARAMETERS:
!  ! TODO - document this.
   real(rk),intent(out) :: j_diff
!
! !INPUT/OUTPUT PARAMETERS:
!  ! TODO - docuemnt these
   real(rk) :: prof(15),profD(15)
!
! !LOCAL VARIABLES:
!  ! TODO - docuemnt these.
   real(rk) :: md,factor
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
!  code.................................................................
!
!  projects content difference p(14) on standard parabolic with 
!    0 surface concentration and 0 bottom flux
      md = prof(8)*profD(11) + prof(9)*profD(12) +prof(10)*profD(13)
      factor = prof(14)/md 

      prof(1) = prof(1) - profD(1)*factor
      prof(2) = prof(2) - profD(2)*factor
      prof(3) = prof(3) - profD(3)*factor

      j_diff = prof(1) + prof(2) + prof(3)

      RETURN

  END SUBROUTINE NonEquFlux
!
!EOC
!-----------------------------------------------------------------------


end module