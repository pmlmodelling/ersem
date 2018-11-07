#include "fabm_driver.h"

module ersem_bacteria_docdyn

   use fabm_types
   use fabm_particle

   use ersem_shared
   use ersem_pelagic_base

   implicit none

   private

   type,extends(type_ersem_pelagic_base),public :: type_ersem_bacteria_docdyn
      ! Variables
      type (type_state_variable_id) :: id_O3c, id_O2o, id_TA
      type (type_state_variable_id) :: id_R1c, id_R2c, id_R3c
      type (type_state_variable_id) :: id_R1p
      type (type_state_variable_id) :: id_R1n
      type (type_state_variable_id) :: id_N1p,id_N4n,id_N7f
      type (type_dependency_id)     :: id_ETW,id_eO2mO2
      type (type_state_variable_id),allocatable,dimension(:) :: id_RPc,id_RPp,id_RPn,id_RPf
      type (type_model_id),         allocatable,dimension(:) :: id_RP
      type (type_diagnostic_variable_id) :: id_fB1O3c, id_fB1NIn, id_fB1N1p
      type (type_diagnostic_variable_id) :: id_fR1B1c, id_fR2B1c, id_fR3B1c,id_fB1R1c, id_fB1R2c, id_fB1R3c
      type (type_diagnostic_variable_id) :: id_fR1B1n,id_fB1R1n,id_fR1B1p,id_fB1R1p
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
      real(rk) :: rR2B1X,rR3B1X
      real(rk),allocatable :: sRPR1(:)
      real(rk) :: frB1R3

      ! Remineralization
      real(rk) :: sR1N1X,sR1N4X
      real(rk) :: fsinkX
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
      class (type_ersem_bacteria_docdyn),intent(inout),target :: self
      integer,                        intent(in)           :: configunit
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

      ! Register links to semi-refractory dissolved organic matter pool.
      call self%register_state_dependency(self%id_R3c,'R3c','mg C/m^3','semi-refractory DOC')

      allocate(self%sRPR1(self%nRP))
      do iRP=1,self%nRP
         write (index,'(i0)') iRP
         call self%get_parameter(self%sRPR1(iRP),'sRP'//trim(index)//'R1','1/d','remineralisation of substrate '//trim(index)//' to DOM')
      end do

      call self%get_parameter(self%rR2B1X,'rR2','-','fraction of semi-labile DOC available to bacteria')
      call self%get_parameter(self%rR3B1X,'rR3','-','fraction of semi-refractory DOC available to bacteria')
      call self%get_parameter(self%frB1R3,'frR3','-','fraction of activity respiration converted to semi-refractory DOC')

      ! Register links to external total dissolved inorganic carbon, dissolved oxygen pools
      call self%register_state_dependency(self%id_O3c,'O3c','mmol C/m^3','carbon dioxide')
      call self%register_state_dependency(self%id_O2o,'O2o','mmol O_2/m^3','oxygen')
      call self%register_state_dependency(self%id_TA,standard_variables%alkalinity_expressed_as_mole_equivalent)

      ! Register environmental dependencies (temperature, suspendend sediment, pH, oxygen saturation)
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      call self%register_dependency(self%id_eO2mO2,standard_variables%fractional_saturation_of_oxygen)

      ! Register diagnostics.
      call self%register_diagnostic_variable(self%id_fB1O3c,'fB1O3c','mg C/m^3/d','respiration')
      call self%register_diagnostic_variable(self%id_fB1NIn,'fB1NIn','mmol N/m^3/d','release of DIN')
      call self%register_diagnostic_variable(self%id_fB1N1p,'fB1N1p','mmol P/m^3/d','release of DIP')

      call self%register_diagnostic_variable(self%id_fB1R1c,'fB1R1c','mg C/m^3/d','release of labile DOC ')
      call self%register_diagnostic_variable(self%id_fB1R2c,'fB1R2c','mg C/m^3/d','release of semi-labile DOC ')
      call self%register_diagnostic_variable(self%id_fB1R3c,'fB1R3c','mg C/m^3/d','release of semi-refractory DOC ')
      call self%register_diagnostic_variable(self%id_fB1R1n,'fB1R1n','mmol N/m^3/d','release of DON')
      call self%register_diagnostic_variable(self%id_fB1R1p,'fB1R1p','mmol P/m^3/d','release of DOP')

      call self%register_diagnostic_variable(self%id_fR1B1c,'fR1B1c','mg C/m^3/d','uptake of labile DOC ')
      call self%register_diagnostic_variable(self%id_fR2B1c,'fR2B1c','mg C/m^3/d','uptake of semi-labile DOC ')
      call self%register_diagnostic_variable(self%id_fR3B1c,'fR3B1c','mg C/m^3/d','uptake of semi-refractory DOC ')
      call self%register_diagnostic_variable(self%id_fR1B1n,'fR1B1n','mmol N/m^3/d','uptake of DON')
      call self%register_diagnostic_variable(self%id_fR1B1p,'fR1B1p','mmol P/m^3/d','uptake of DOP')

      call self%register_diagnostic_variable(self%id_minn,'minn','mmol N/m^3/d','mineralisation of DON to DIN')
      call self%register_diagnostic_variable(self%id_minp,'minp','mmol P/m^3/d','mineralisation of DOP to DIP')
   end subroutine

   subroutine do(self,_ARGUMENTS_DO_)

      class (type_ersem_bacteria_docdyn),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

   ! !LOCAL VARIABLES:
      real(rk) :: ETW,eO2mO2
      real(rk) :: B1c,B1n,B1p
      real(rk) :: B1cP,B1nP,B1pP
      real(rk) :: N1pP,N4nP,R1c,R1cP,R1pP,R1nP,R2c
      real(rk) :: qpB1c,qnB1c
      real(rk) :: etB1,eO2B1
      real(rk) :: sB1RD,sutB1,rumB1,sugB1,rugB1,rraB1,fB1O3c
      real(rk) :: sB1R2,fB1R2c,fB1R3c,fB1RDc
      real(rk) :: netb1,bge
      real(rk) :: fB1N1p,fR1B1p,fB1RDp
      real(rk) :: fB1NIn,fR1B1n,fB1RDn
      real(rk) :: R3c,R2cP,R3cP
      real(rk) :: fB1R1c
      real(rk) :: totsubst
      real(rk) :: CORROX
      integer  :: iRP
      real(rk),dimension(self%nRP) :: RPc,RPcP,RPnP,RPpP
      real(rk),dimension(self%nRP) :: fRPB1c,fRPB1p,fRPB1n

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

         _GET_WITH_BACKGROUND_(self%id_R3c,R3c)
         _GET_(self%id_R2c,R2cP)
         _GET_(self%id_R3c,R3cP)
         do iRP=1,self%nRP
            _GET_WITH_BACKGROUND_(self%id_RPc(iRP),RPc(iRP))
            _GET_(self%id_RPc(iRP),RPcP(iRP))
            _GET_(self%id_RPn(iRP),RPnP(iRP))
            _GET_(self%id_RPp(iRP),RPpP(iRP))
         end do

         qpB1c = B1p/B1c
         qnB1c = B1n/B1c
!..Temperature effect on pelagic bacteria:

         etB1 = max(0.0_rk,self%q10B1X**((ETW-10._rk)/10._rk) - self%q10B1X**((ETW-32._rk)/3._rk))

!..Prevailing Oxygen limitation for bacteria:

         CORROX = 1._rk + self%chdB1oX
         eO2B1 = min(1._rk,CORROX*(eO2mO2/( eO2mO2 + self%chdB1oX )))

!..bacterial mortality

         sB1RD=self%sdB1X*etB1
         fB1R1c = sB1RD*B1cP

!..Total amount of substrate available

         sutB1 = self%sR1B1X  ! DOM-specific uptake rate
         ! rutB1 = sutB1*R1cP(I)

!..Potential uptake :
         rumB1 = self%sumB1X*etB1*eO2B1*B1c

!..Actual uptake

      ! rugB1 = MIN(rumB1,rutB1)
      ! specific in substrate concentration:

      totsubst = R1cP+R2cP*self%rR2B1X+R3cP*self%rR3B1X+sum(RPcP*self%sRPR1/sutB1)
      ! Jorn: check whether total substrate>0 to prevent NaNs
      if (totsubst>0.0_rk) then
         sugB1 = rumB1/max(rumB1/sutB1,totsubst)
      else
         sugB1 = 0.0_rk
      end if
            ! = MIN(rumB1,rutB1)=MIN(rumB1/(R1cP+R2cP*rR2B1X,sutB1) avoid pot. div. by 0
      fRPB1c = sugB1*RPcP*self%sRPR1/sutB1
      rugB1 = sugB1*(R1cP+R2cP*self%rR2B1X+R3cP*self%rR3B1X)+sum(fRPB1c)

!..Respiration :

! activity respiration
          rraB1 = rugB1 * ( 1._rk - self%puB1X*eO2mO2  - self%puB1oX*( 1._rk - eO2mO2 ) )
! total respiration:
          fB1O3c = rraB1 + self%srsB1X * B1cP * etB1
!      fB1O3c(I) = ( 1._fp8 - puB1X*eO2mO2  - puB1oX*( 1._fp8 - eO2mO2 ) )&
!                * rugB1  + srsB1X * B1cP(I) * etB1
          _SET_DIAGNOSTIC_(self%id_fB1O3c,fB1O3c)

! specific release of semilabile DOC
! fudge factor 1. as used in Polimene et al. AME 2006
         sB1R2=max(0._rk,max(1._rk-(qpB1c/self%qpB1cX),1._rk-(qnB1c/self%qnB1cX)))*1._rk
         fB1R2c=sB1R2*B1cP
         fB1R3c=self%frB1R3*rraB1
         fB1RDc = fB1R1c + fB1R2c + fB1R3c

!..net bacterial production

         netb1 = rugB1 - fB1o3c - fB1RDc
         IF (netB1.gt.0._rk) THEN
            BGE=netB1/rugB1
         ELSE
            BGE=0._rk
         ENDIF

!..Source equations


         _SET_ODE_(self%id_c,netb1)
         _SET_ODE_(self%id_R1c,+ fB1R1c - sugB1*R1cP)
         _SET_ODE_(self%id_R2c,+ fB1R2c - sugB1*R2cP*self%rR2B1X)
         _SET_ODE_(self%id_R3c,+ fB1R3c - sugB1*R3cP*self%rR3B1X)

         _SET_DIAGNOSTIC_(self%id_fB1R1c, fB1R1c)
         _SET_DIAGNOSTIC_(self%id_fB1R2c, fB1R2c)
         _SET_DIAGNOSTIC_(self%id_fB1R3c, fB1R3c)
         _SET_DIAGNOSTIC_(self%id_fR1B1c, sugB1*R1cP)
         _SET_DIAGNOSTIC_(self%id_fR2B1c, sugB1*R2cP*self%rR2B1X)
         _SET_DIAGNOSTIC_(self%id_fR3B1c, sugB1*R3cP*self%rR3B1X)

         do iRP=1,self%nRP
            _SET_ODE_(self%id_RPc(iRP), -fRPB1c(iRP))
         end do

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

         fRPB1p = sugB1*RPpP*self%sRPR1/sutB1
         _SET_ODE_(self%id_p, sum(fRPB1p))
         do iRP=1,self%nRP
            _SET_ODE_(self%id_RPp(iRP), - fRPB1p(iRP))
         end do

         _SET_ODE_(self%id_p, + fR1B1p - fB1N1p - fB1RDp)
         _SET_ODE_(self%id_N1p, + fB1N1p)
         _SET_ODE_(self%id_R1p, + fB1RDp - fR1B1p)
         _SET_ODE_(self%id_TA,  - fB1N1p)   ! Contribution to alkalinity: -1 for phosphate

!..Set diagnostics
         _SET_DIAGNOSTIC_(self%id_fB1N1p,fB1N1p)
         _SET_DIAGNOSTIC_(self%id_fB1R1p, fB1RDp)
         _SET_DIAGNOSTIC_(self%id_fR1B1p, fR1B1p)

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

         fRPB1n = sugB1*RPnP*self%sRPR1/sutB1
         _SET_ODE_(self%id_n, sum(fRPB1n))
         do iRP=1,self%nRP
            _SET_ODE_(self%id_RPn(iRP), - fRPB1n(iRP))
         end do

         _SET_ODE_(self%id_N4n, + fB1NIn)
         _SET_ODE_(self%id_n,   + fR1B1n - fB1NIn - fB1RDn)
         _SET_ODE_(self%id_R1n, + fB1RDn   - fR1B1n)
         _SET_ODE_(self%id_TA,  + fB1NIn)   ! Contribution to alkalinity: +1 for ammonium

!..Set diagnostics
         _SET_DIAGNOSTIC_(self%id_fB1NIn,fB1NIn)
         _SET_DIAGNOSTIC_(self%id_fB1R1n, fB1RDn)
         _SET_DIAGNOSTIC_(self%id_fR1B1n, fR1B1n)

      ! Leave spatial loops (if any)
      _LOOP_END_

      call remineralization(self,_ARGUMENTS_DO_)

   end subroutine do


   subroutine remineralization(self,_ARGUMENTS_DO_)

      class (type_ersem_bacteria_docdyn),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      integer :: iRP
      real(rk) :: R1pP,R1nP
      real(rk) :: N7fP
      real(rk) :: fR1N1p,fR1NIn
      real(rk) :: RPfP(self%nRP),fRPN7f(self%nRP),n7fsink

      _LOOP_BEGIN_

         ! Mineralisation of DOP to PO4
         _GET_(self%id_R1p,R1pP)
         fR1N1p = self%sR1N1X * R1pP

         ! Mineralisation of DON to NH4
         _GET_(self%id_R1n,R1nP)
         fR1NIn = self%sR1N4X * R1nP

         !..Source equations
         _SET_ODE_(self%id_R1p, - fR1N1p)
         _SET_ODE_(self%id_R1n, - fR1NIn)

         _SET_ODE_(self%id_N1p, + fR1N1p)
         _SET_ODE_(self%id_N4n, + fR1NIn)
         _SET_ODE_(self%id_TA,  - fR1N1p + fR1NIn)   ! Contributions to alkalinity: -1 for phosphate, +1 for ammonium

         !.. set Diagnostics
         _SET_DIAGNOSTIC_(self%id_minn,fR1NIn)
         _SET_DIAGNOSTIC_(self%id_minp,fR1N1p)

         if (use_iron) then
            ! remineralization of particulate iron to Fe
            do iRP=1,self%nRP
               _GET_(self%id_RPf(iRP),RPfP(iRP))
            end do
            fRPN7f=self%sRPr1*RPfP

            ! sink of Fe

            ! This term takes into account the scavenging due to hydroxide precipitation and it is supposed to be
            ! regulated by a threshold concentration (0.6 nM). See Aumont et al., 2003 (GBC) and Vichi et al., 2007 (JMS) for references.
            ! (Luca, 12/08)
            _GET_(self%id_N7f,N7fP)
            n7fsink=self%fsinkX*max(0._rk,N7fP-0.6_rk)

            do iRP=1,self%nRP
               if (fRPn7f(iRP)/=0.0_rk) _SET_ODE_(self%id_RPf(iRP),-fRPn7f(iRP))
            end do
            _SET_ODE_(self%id_N7f,+sum(fRPn7f)-n7fsink)
         end if

      ! Leave spatial loops (if any)
      _LOOP_END_

   end subroutine remineralization

end module
