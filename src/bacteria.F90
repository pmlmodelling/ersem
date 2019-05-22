#include "fabm_driver.h"

module ersem_bacteria

   use fabm_types
   use fabm_particle

   use ersem_shared
   use ersem_pelagic_base

   implicit none

   private

   type,extends(type_ersem_pelagic_base),public :: type_ersem_bacteria
      ! Variables
      type (type_state_variable_id) :: id_O3c, id_O2o, id_TA
      type (type_state_variable_id) :: id_R1c, id_R1p, id_R1n
      type (type_state_variable_id) :: id_R2c
      type (type_state_variable_id) :: id_N1p,id_N4n,id_N7f
      type (type_dependency_id)     :: id_ETW,id_eO2mO2
      type (type_state_variable_id),allocatable,dimension(:) :: id_RPc,id_RPp,id_RPn,id_RPf
      type (type_model_id),         allocatable,dimension(:) :: id_RP
      type (type_diagnostic_variable_id) :: id_fB1O3c, id_fB1NIn, id_fB1N1p
      type (type_diagnostic_variable_id) :: id_minn,id_minp

      ! Parameters
      integer  :: nRP
      integer  :: iswBlimX
      real(rk) :: q10B1X,chdB1oX
      real(rk) :: chB1nX,chB1pX
      real(rk) :: sdB1X
      real(rk) :: sumB1X
      real(rk) :: puB1X,puB1oX,srsB1X,sR1B1X
      real(rk) :: qpB1cX,qnB1cX
      real(rk) :: urB1_O2X
      real(rk) :: R1R2X
      real(rk),allocatable :: puRP_B1X(:)

      ! Remineralization
      real(rk) :: sR1N1X,sR1N4X
      real(rk) :: fsinkX
      real(rk) :: redfieldX
      real(rk) :: rR2R1X
   contains
!     Model procedures
      procedure :: initialize
      procedure :: do
   end type

   real(rk),parameter :: onedayX = 1.0_rk

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_ersem_bacteria),intent(inout),target :: self
      integer,                    intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
      integer           :: iRP
      character(len=16) :: index
      real(rk)          :: c0
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%iswBlimX,'iswBlim', '',           'nutrient limitation (1: minimum of inorganic and organic availability, 2: additive availability)')
      call self%get_parameter(self%q10B1X,  'q10',     '-',          'Q_10 temperature coefficient')
      call self%get_parameter(self%chdB1oX, 'chdo',    '-',          'Michaelis-Menten constant for oxygen limitation')
      call self%get_parameter(self%chB1nX,  'chn',     'mmol N/m^3', 'Michaelis-Menten constant for nitrate limitation')
      call self%get_parameter(self%chB1pX,  'chp',     'mmol P/m^3', 'Michaelis-Menten constant for phosphate limitation')
      call self%get_parameter(self%sdB1X,   'sd',      '1/d',        'specific mortality at reference temperature')
      call self%get_parameter(self%sumB1X,  'sum',     '1/d',        'maximum specific uptake at reference temperature')
      call self%get_parameter(self%puB1X,   'pu',      '-',          'efficiency at high oxygen levels')
      call self%get_parameter(self%puB1oX,  'puo',     '-',          'efficiency at low oxygen levels')
      call self%get_parameter(self%srsB1X,  'srs',     '1/d',        'specific rest respiration at reference temperature')
      call self%get_parameter(self%sR1B1X,  'sR1',     '1/d',        'maximum turn-over rate of DOM', default=1.0_rk)
      call self%get_parameter(self%qpB1cX,  'qpc',     'mmol P/mg C','maximum phosphorus to carbon ratio')
      call self%get_parameter(self%qnB1cX,  'qnc',     'mmol N/mg C','maximum nitrogen to carbon ratio')
      call self%get_parameter(self%urB1_O2X,'ur_O2',   'mmol O_2/mg C','oxygen consumed per carbon respired')

      ! Remineralization parameters
      call self%get_parameter(self%sR1N1X,   'sR1N1',   '1/d',    'mineralisation rate of labile dissolved organic phosphorus')
      call self%get_parameter(self%sR1N4X,   'sR1N4',   '1/d',    'mineralisation rate of labile dissolved organic nitrogen')
      call self%get_parameter(self%fsinkX,   'fsink',   '1/d',    'scavenging rate for iron')

      call self%get_parameter(self%redfieldX,'redfield','mol/mol','Redfield carbon to nitrogen ratio')
      call self%get_parameter(self%rR2R1X,   'rR2R1', '1/d', 'specific rate of breakdown of semi-labile to labile DOC')

      call self%get_parameter(c0,'c0','mg C/m^3','background carbon concentration')

      ! Allow ERSEM base model to declare our own state variables.
      call self%initialize_ersem_base(sedimentation=.false.)
      call self%add_constituent('c',1.e-4_rk,   c0)
      call self%add_constituent('n',1.26e-6_rk, self%qnB1cX*c0)
      call self%add_constituent('p',4.288e-8_rk,self%qpB1cX*c0)

      ! Register links to nutrient pools.
      call self%register_state_dependency(self%id_N1p,'N1p','mmol P/m^3','phosphate')
      call self%register_state_dependency(self%id_N4n,'N4n','mmol N/m^3','ammonium')
      if (use_iron) call self%register_state_dependency(self%id_N7f,'N7f','umol Fe/m^3','inorganic iron')

      ! Register links to labile dissolved organic matter pools.
      call self%register_state_dependency(self%id_R1c,'R1c','mg C/m^3',  'labile dissolved organic carbon')
      call self%register_state_dependency(self%id_R1p,'R1p','mmol P/m^3','labile dissolved organic phosphorus')
      call self%register_state_dependency(self%id_R1n,'R1n','mmol N/m^3','labile dissolved organic nitrogen')

      ! Register links to semi-labile dissolved organic matter pools.
      call self%register_state_dependency(self%id_R2c,'R2c','mg C/m^3','semi-labile dissolved organic carbon')

      ! Register links to particulate organic matter pools.
      call self%get_parameter(self%nRP,'nRP','','number of substrates',default=0)
      allocate(self%id_RP(self%nRP))
      allocate(self%id_RPc(self%nRP))
      allocate(self%id_RPn(self%nRP))
      allocate(self%id_RPp(self%nRP))
      allocate(self%id_RPf(self%nRP))
      do iRP=1,self%nRP
         write (index,'(i0)') iRP
         call self%register_state_dependency(self%id_RPc(iRP),'RP'//trim(index)//'c','mg C/m^3',   'carbon in substrate '//trim(index))
         call self%register_state_dependency(self%id_RPn(iRP),'RP'//trim(index)//'n','mmol N/m^3', 'nitrogen in substrate '//trim(index))
         call self%register_state_dependency(self%id_RPp(iRP),'RP'//trim(index)//'p','mmol P/m^3', 'phosphorus in substrate '//trim(index))
         call self%register_model_dependency(self%id_RP(iRP),'RP'//trim(index))
         call self%request_coupling_to_model(self%id_RPc(iRP),self%id_RP(iRP),'c')  ! For now link to hardcoded "c" to get a direct link to state mg C/m3 (and not a diagnostic for mmol C/m3)
         call self%request_coupling_to_model(self%id_RPn(iRP),self%id_RP(iRP),standard_variables%total_nitrogen)
         call self%request_coupling_to_model(self%id_RPp(iRP),self%id_RP(iRP),standard_variables%total_phosphorus)
         if (use_iron) then
            call self%register_state_dependency(self%id_RPf(iRP),'RP'//trim(index)//'f','umol Fe/m^3','iron in substrate '//trim(index))
            call self%request_coupling_to_model(self%id_RPf(iRP),self%id_RP(iRP),standard_variables%total_iron)
         end if
      end do

      allocate(self%puRP_B1X(self%nRP))
      do iRP=1,self%nRP
         write (index,'(i0)') iRP
         call self%get_parameter(self%puRP_B1X(iRP),'puRP'//trim(index),'-','fraction of substrate '//trim(index)//' available to bacteria')
      end do

      call self%get_parameter(self%R1R2X,'R1R2','-','labile fraction of produced dissolved organic carbon')

      ! Register links to external total dissolved inorganic carbon, dissolved oxygen pools
      call self%register_state_dependency(self%id_O2o,'O2o','mmol O_2/m^3','oxygen')
      call self%register_state_dependency(self%id_O3c,'O3c','mmol C/m^3','carbon dioxide')
      call self%register_state_dependency(self%id_TA,standard_variables%alkalinity_expressed_as_mole_equivalent)

      ! Register environmental dependencies (temperature, suspendend sediment, pH, oxygen saturation)
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      call self%register_dependency(self%id_eO2mO2,standard_variables%fractional_saturation_of_oxygen)

      ! Register diagnostics.
      call self%register_diagnostic_variable(self%id_fB1O3c,'fB1O3c','mg C/m^3/d','respiration')
      call self%register_diagnostic_variable(self%id_fB1NIn,'fB1NIn','mmol N/m^3/d','release of DIN')
      call self%register_diagnostic_variable(self%id_fB1N1p,'fB1N1p','mmol P/m^3/d','release of DIP')
      call self%register_diagnostic_variable(self%id_minn,'minn','mmol N/m^3/d','mineralisation of N')
      call self%register_diagnostic_variable(self%id_minp,'minp','mmol P/m^3/d','mineralisation of P')

   end subroutine

   subroutine do(self,_ARGUMENTS_DO_)

      class (type_ersem_bacteria),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

   ! !LOCAL VARIABLES:
      real(rk) :: ETW,eO2mO2
      real(rk) :: B1c,B1n,B1p
      real(rk) :: B1cP,B1nP,B1pP
      real(rk) :: N1pP,N4nP,R1c,R1cP,R1pP,R1nP,R2c
      real(rk) :: qpB1c,qnB1c
      real(rk) :: etB1,eO2B1
      real(rk) :: sB1RD,sutB1,rumB1,sugB1,rugB1,rraB1,fB1O3c
      real(rk) :: fB1RDc
      real(rk) :: netb1,bge
      real(rk) :: fB1N1p,fR1B1p,fB1RDp
      real(rk) :: fB1NIn,fR1B1n,fB1RDn
      real(rk) :: Nlim,Plim
      real(rk) :: CORROX

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         _GET_(self%id_ETW,ETW)
         _GET_(self%id_eO2mO2,eO2mO2)
         eO2mO2 = min(1.0_rk,eO2mO2)

         _GET_WITH_BACKGROUND_(self%id_c,B1c)
         _GET_WITH_BACKGROUND_(self%id_p,B1p)
         _GET_WITH_BACKGROUND_(self%id_n,B1n)
         _GET_(self%id_c,B1cP)
         _GET_(self%id_p,B1pP)
         _GET_(self%id_n,B1nP)

         _GET_(self%id_N1p,N1pP)
         _GET_(self%id_N4n,N4nP)
         _GET_WITH_BACKGROUND_(self%id_R1c,R1c)
         _GET_WITH_BACKGROUND_(self%id_R2c,R2c)
         _GET_(self%id_R1c,R1cP)
         _GET_(self%id_R1p,R1pP)
         _GET_(self%id_R1n,R1nP)

         qpB1c = B1p/B1c
         qnB1c = B1n/B1c
!..Temperature effect on pelagic bacteria:

         etB1 = self%q10B1X**((ETW-10._rk)/10._rk) - self%q10B1X**((ETW-32._rk)/3._rk)

!..Prevailing Oxygen limitation for bacteria:

         CORROX = 1._rk + self%chdB1oX
         eO2B1 = min(1._rk,CORROX*(eO2mO2/( eO2mO2 + self%chdB1oX )))

!..Potential nutrient limitation on bacterial growth after Anderson

         IF (self%iswBlimX .EQ. 1 ) THEN
           Nlim = min(N4nP/(N4nP+self%chB1nX),R1nP/(R1nP+self%chB1nX))
           Plim = min(N1pP/(N1pP+self%chB1pX),R1pP/(R1pP+self%chB1pX))
         ELSE IF (self%iswBlimX .EQ. 2 ) THEN
           Nlim = (N4nP+R1nP)/(N4nP+R1nP+self%chB1nX)
           Plim = (N1pP+R1pP)/(N1pP+R1pP+self%chB1pX)
         END IF

!..bacterial mortality

         sB1RD=self%sdB1X*etB1
         fB1RDc = sB1RD*B1cP

!..Total amount of substrate available

         sutB1 = self%sR1B1X  ! DOM-specific uptake rate
         ! rutB1 = sutB1*R1cP(I)

!..Potential uptake :
         rumB1 = self%sumB1X*etB1*eO2B1*min(Nlim,Plim)*B1c

!..Actual uptake

      ! rugB1 = MIN(rumB1,rutB1)
      ! specific in substrate concentration:
         sugB1 = rumB1/max(rumB1/sutB1,R1c)
             ! = MIN(rumB1,rutB1)=MIN(rumB1/R1cP,sutB1) avoid pot. div. by 0
         rugB1 = sugB1*R1cP

!..Respiration :

! activity respiration
          rraB1 = rugB1 * ( 1._rk - self%puB1X*eO2mO2  - self%puB1oX*( 1._rk - eO2mO2 ) )
! total respiration:
          fB1O3c = rraB1 + self%srsB1X * B1cP * etB1
!      fB1O3c(I) = ( 1._fp8 - puB1X*eO2mO2  - puB1oX*( 1._fp8 - eO2mO2 ) )&
!                * rugB1  + srsB1X * B1cP(I) * etB1
          _SET_DIAGNOSTIC_(self%id_fB1O3c,fB1O3c)

!..net bacterial production

         netb1 = rugB1 - fB1o3c - fB1RDc
         IF (netB1.gt.0._rk) THEN
            BGE=netB1/rugB1
         ELSE
            BGE=0._rk
         ENDIF

!..Source equations


         _SET_ODE_(self%id_c,netb1)
         _SET_ODE_(self%id_R1c,- rugB1 + (fB1RDc * self%R1R2X))
         _SET_ODE_(self%id_R2c,+ (fB1RDc * (1._rk - self%R1R2X)))

         _SET_ODE_(self%id_O3c,+ fB1O3c/CMass)
         _SET_ODE_(self%id_O2o,- fB1O3c*self%urB1_O2X)

!..Phosphorus dynamics in bacteria........................................

         IF ((qpB1c - self%qpB1cX).gt.0._rk) THEN
           fB1N1p = (qpB1c - self%qpB1cX ) * B1cP /onedayX ! sink
         ELSE
           fB1N1p = (qpB1c - self%qpB1cX ) * B1c *&
                           N1pP/(N1pP+self%chB1pX) /onedayX ! source
         ENDIF

!..uptake of DOP

         fR1B1p = sugB1*R1pP

!..flux of DOP from B1

         fB1RDp = sB1RD*B1pP


!..Source equations

         _SET_ODE_(self%id_p, + fR1B1p - fB1N1p - fB1RDp)
         _SET_ODE_(self%id_N1p, + fB1N1p)
         _SET_ODE_(self%id_R1p, + fB1RDp - fR1B1p)
         _SET_ODE_(self%id_TA,  - fB1N1p)   ! Contribution to alkalinity: -1 for phosphate

!..Set diagnostics
         _SET_DIAGNOSTIC_(self%id_fB1N1p,fB1N1p)

!..Nitrogen dynamics in bacteria........................................

         IF ((qnB1c - self%qnB1cX).gt.0._rk) THEN
           fB1NIn = (qnB1c - self%qnB1cX ) * B1cP /onedayX ! sink
         ELSE
           fB1NIn = (qnB1c - self%qnB1cX ) * B1c *    &
                           N4nP/(N4nP+self%chB1nX) /onedayX ! source
         ENDIF

!..uptake of DON

         fR1B1n = sugB1*R1nP

!..flux of DON from B1

         fB1RDn = sB1RD*B1nP

!..Source equations

         _SET_ODE_(self%id_N4n, + fB1NIn)
         _SET_ODE_(self%id_n,   + fR1B1n - fB1NIn - fB1RDn)
         _SET_ODE_(self%id_R1n, + fB1RDn   - fR1B1n)
         _SET_ODE_(self%id_TA,  + fB1NIn)   ! Contribution to alkalinity: +1 for ammonium

!..Set diagnostics
         _SET_DIAGNOSTIC_(self%id_fB1NIn,fB1NIn)

      ! Leave spatial loops (if any)
      _LOOP_END_

      call remineralization(self,_ARGUMENTS_DO_)

   end subroutine do

   subroutine remineralization(self,_ARGUMENTS_DO_)

      class (type_ersem_bacteria),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      integer :: iRP
      real(rk) :: R2cP
      real(rk),dimension(self%nRP) :: RPc,RPn,RPcP,RPpP,RPnP
      real(rk),dimension(self%nRP) :: qnRPc
#ifdef NOBAC
      real(rk) :: R1cP
      real(rk) :: fB1O3c
      real(rk) :: ETW,etRemin
#endif
      real(rk) :: R1pP,R1nP
      real(rk) :: N7fP
      real(rk) :: fR1N1p,fR1NIn
      real(rk) :: sRPr1(self%nRP)
      real(rk) :: RPfP(self%nRP),n7fsink

      _LOOP_BEGIN_

         ! -----------------------------------------
         ! POM breakdown
         ! -----------------------------------------

         ! Get particulate organic matter concentrations. Also retrieve POC and PON including background,
         ! so N:C can be computed without risk of division by 0.
         do iRP=1,self%nRP
            _GET_(self%id_RPc(iRP),RPcP(iRP))
            _GET_(self%id_RPp(iRP),RPpP(iRP))
            _GET_(self%id_RPn(iRP),RPnP(iRP))
            _GET_WITH_BACKGROUND_(self%id_RPc(iRP),RPc(iRP))
            _GET_WITH_BACKGROUND_(self%id_RPn(iRP),RPn(iRP))
         end do

         ! Specific rate of particulate matter conversion (1/d) - substrate-specific and proportional to N:C
         qnRPc = RPn/RPc
         sRPr1  = self%redfieldX*CMass * qnRPc * self%puRP_B1X

#ifdef NOBAC
         ! Apply temperature limitation factor to particulate organic matter conversion rates.
         ! NB limitation factor is also used further down for mineralization of dissolved organic matter.
         _GET_(self%id_ETW,ETW)
         etRemin = max(0.0_rk,self%q10B1X**((ETW-10._rk)/10._rk) - self%q10B1X**((ETW-32._rk)/3._rk))
         sRPr1 = sRPr1*etRemin
#endif

         ! Decrease in particulate organic matter due to conversion to dissolved organic matter.
         do iRP=1,self%nRP
            _SET_ODE_(self%id_RPc(iRP), -sRPr1(iRP)*RPcP(iRP))
            _SET_ODE_(self%id_RPp(iRP), -sRPr1(iRP)*RPpP(iRP))
            _SET_ODE_(self%id_RPn(iRP), -sRPr1(iRP)*RPnP(iRP))
         end do

         ! Increase in dissolved organic matter due to conversion from particulate organic matter.
         _SET_ODE_(self%id_R1c, + sum(sRPr1*RPcP))
         _SET_ODE_(self%id_R1p, + sum(sRPr1*RPpP))
         _SET_ODE_(self%id_R1n, + sum(sRPr1*RPnP))

         ! Conversion of particulate iron to inorganic Fe
         if (use_iron) then
            do iRP=1,self%nRP
               _GET_(self%id_RPf(iRP),RPfP(iRP))
               if (RPfP(iRP)/=0.0_rk) _SET_ODE_(self%id_RPf(iRP),-sRPr1(iRP)*RPfP(iRP))
            end do
            _SET_ODE_(self%id_N7f, + sum(sRPr1*RPfP))
         end if

         ! -----------------------------------------
         ! DOM mineralization
         ! -----------------------------------------

#ifdef NOBAC
         ! Mineralisation of dissolved organic carbon to CO2
         _GET_(self%id_R1c,R1cP)
         fB1O3c = self%sR1N1X * R1cP*etRemin
         _SET_ODE_(self%id_R1C, - fB1O3c)
         _SET_ODE_(self%id_O3c, + fB1O3c/CMass)
         _SET_ODE_(self%id_O2o, - fB1O3c*self%urB1_O2X)
         _SET_DIAGNOSTIC_(self%id_fB1O3c,fB1O3c)
#endif

         ! Mineralisation of dissolved organic phosphorus to PO4
         _GET_(self%id_R1p,R1pP)
         fR1N1p = self%sR1N1X * R1pP

         ! Mineralisation of dissolved organic nitrogen to NH4
         _GET_(self%id_R1n,R1nP)
         fR1NIn = self%sR1N4X * R1nP

#ifdef NOBAC
         ! Apply temperature limitation factor to DOP, DON mineralization.
         fR1N1p = fR1N1p*etRemin
         fR1NIn = fR1NIn*etRemin
#endif

         ! Fluxes due to mineralization of DOP and DON to PO4 and NH4.
         _SET_ODE_(self%id_R1p, - fR1N1p)
         _SET_ODE_(self%id_R1n, - fR1NIn)
         _SET_ODE_(self%id_N1p, + fR1N1p)
         _SET_ODE_(self%id_N4n, + fR1NIn)
         _SET_ODE_(self%id_TA,  - fR1N1p + fR1NIn)   ! Contributions to alkalinity: -1 for phosphate, +1 for ammonium

         !.. set Diagnostics
         _SET_DIAGNOSTIC_(self%id_minn,fR1NIn)
         _SET_DIAGNOSTIC_(self%id_minp,fR1N1p)

         ! -----------------------------------------
         ! Iron scavenging
         ! -----------------------------------------

         if (use_iron) then
            ! This term takes into account the scavenging due to hydroxide precipitation and it is supposed to be
            ! regulated by a threshold concentration (0.6 nM). See Aumont et al., 2003 (GBC) and Vichi et al., 2007 (JMS) for references.
            ! (Luca, 12/08)
            _GET_(self%id_N7f,N7fP)
            n7fsink = self%fsinkX*max(0._rk,N7fP-0.6_rk)
            _SET_ODE_(self%id_N7f,-n7fsink)
         end if

         ! -----------------------------------------
         ! Conversion from semi-labile to labile DOM
         ! -----------------------------------------

         _GET_(self%id_R2c,R2cP)
         _SET_ODE_(self%id_R1c, + self%rR2R1X * R2cP)
         _SET_ODE_(self%id_R2c, - self%rR2R1X * R2cP)

      ! Leave spatial loops (if any)
      _LOOP_END_

   end subroutine remineralization

end module
