#include "fabm_driver.h"

!#define DOCDYN
!#define NOBAC
!#define IRON

module pml_ersem_bacteria

   use fabm_types
   use fabm_particle

   use pml_ersem_shared
   use pml_ersem_pelagic_base

   implicit none

   private

   type,extends(type_ersem_pelagic_base_model),public :: type_pml_ersem_bacteria
      ! Variables
      type (type_state_variable_id) :: id_O3c, id_O2o
      type (type_state_variable_id) :: id_R1c, id_R2c, id_R3c
      type (type_state_variable_id) :: id_R1p
      type (type_state_variable_id) :: id_R1n
      type (type_state_variable_id) :: id_N1p,id_N3n,id_N4n,id_N7f
      type (type_dependency_id)     :: id_ETW,id_ESS,id_phx,id_eO2mO2
      type (type_state_variable_id),allocatable,dimension(:) :: id_RPc,id_RPp,id_RPn,id_RPf
      type (type_model_id),         allocatable,dimension(:) :: id_RP

      ! Parameters
      integer  :: nRP
      integer  :: iswBlimX
      real(rk) :: q10B1X,chdB1oX
      real(rk) :: chB1nX,chB1pX
      real(rk) :: sdB1X
      real(rk) :: sumB1X
      real(rk) :: puB1X,puB1oX,srsB1X
      real(rk) :: qpB1cX,qnB1cX
      real(rk) :: urB1_O2X
#ifdef DOCDYN
      real(rk) :: rR2B1X,rR3B1X
      real(rk),allocatable :: sRPR1(:)
      real(rk) :: frB1R3
#else
      real(rk) :: R1R2X
      real(rk),allocatable :: puRP_B1X(:)
#endif

      ! Remineralization
      real(rk) :: redfieldX
      real(rk) :: sR1N1X,sR1N4X
      real(rk) :: fsinkX
      real(rk) :: sN4N3X,cessX
      integer  :: ISWphx
#ifndef DOCDYN
      real(rk) :: rR2R1X
#endif
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
      class (type_pml_ersem_bacteria),intent(inout),target :: self
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
      call self%get_parameter(self%iswBlimX,'iswBlim', '',           'nutrient limitation (0=minimum of inorganic and organic availability,1=additive availability)')
      call self%get_parameter(self%q10B1X,  'q10',     '-',          'Q_10 temperature coefficient')
      call self%get_parameter(self%chdB1oX, 'chdo',    '-',          'Michaelis-Menten constant for oxygen limitation')
      call self%get_parameter(self%chB1nX,  'chn',     'mmol N/m^3', 'Michaelis-Menten constant for nitrate limitation')
      call self%get_parameter(self%chB1pX,  'chp',     'mmol P/m^3', 'Michaelis-Menten constant for phosphate limitation')
      call self%get_parameter(self%sdB1X,   'sd',      '1/d',        'specific mortality at reference temperature')
      call self%get_parameter(self%sumB1X,  'sum',     '1/d',        'maximum specific uptake at reference temperature')
      call self%get_parameter(self%puB1X,   'pu',      '-',          'efficiency at high oxygen levels')
      call self%get_parameter(self%puB1oX,  'puo',     '-',          'efficiency at low oxygen levels')
      call self%get_parameter(self%srsB1X,  'srs',     '1/d',        'specific rest respiration at reference temperature')
      call self%get_parameter(self%qpB1cX,  'qpc',     'mmol P/mg C','maximum phophorus to carbon ratio')
      call self%get_parameter(self%qnB1cX,  'qnc',     'mmol N/mg C','maximum nitrogen to carbon ratio')
      call self%get_parameter(self%urB1_O2X,'ur_O2',   'mmol O_2/mg C','oxygen consumed per carbon respired')

      ! Remineralization parameters
      call self%get_parameter(self%redfieldX,'redfield','mol/mol','Redfield carbon to nitrogen ratio')
      call self%get_parameter(self%sR1N1X,   'sR1N1',   '1/d',    'mineralisation rate of labile dissolved organic phosphorus')
      call self%get_parameter(self%sR1N4X,   'sR1N4',   '1/d',    'mineralisation rate of labile dissolved organic nitrogen')
      call self%get_parameter(self%fsinkX,   'fsink',   '1/d',    'scavenging rate for iron')
      call self%get_parameter(self%ISWphx,   'ISWph',   '',       'switch for pH impact on nitrification')
      call self%get_parameter(self%sN4N3X,   'sN4N3',   '1/d',    'specific nitrification rate')
      call self%get_parameter(self%cessX,    'cess',    'mg/m^3', 'silt concentration at which relative rate of nitrification is 1')
      
#ifndef DOCDYN
      call self%get_parameter(self%rR2R1X,   'rR2R1', '1/d', 'specific rate of breakdown of semi-labile to labile DOC')
#endif

      call self%get_parameter(c0,'c0','mg C/m^3','background carbon concentration')

      ! Allow ERSEM base model to declare our own state variables.
      call self%initialize_ersem_base(sedimentation=.false.)
      call self%add_constituent('c',1.e-4_rk,   c0)
      call self%add_constituent('n',1.26e-6_rk, self%qnB1cX*c0)
      call self%add_constituent('p',4.288e-8_rk,self%qpB1cX*c0)

      ! Register links to nutrient pools.
      call self%register_state_dependency(self%id_N1p,'N1p','mmol P/m^3', 'Phosphate')
      call self%register_state_dependency(self%id_N3n,'N3n','mmol N/m^3', 'Nitrate')
      call self%register_state_dependency(self%id_N4n,'N4n','mmol N/m^3', 'Ammonium')
#ifdef IRON
      call self%register_state_dependency(self%id_N7f,'N7f','umol Fe/m^3','Inorganic Iron')
#endif

      ! Register links to labile dissolved organic matter pools.
      call self%register_state_dependency(self%id_R1c,'R1c','mg C/m^3',  'DOC')
      call self%register_state_dependency(self%id_R1p,'R1p','mmol P/m^3','DOP')    
      call self%register_state_dependency(self%id_R1n,'R1n','mmol N/m^3','DON')    

      ! Register links to semi-labile dissolved organic matter pools.
      call self%register_state_dependency(self%id_R2c,'R2c','mg C/m^3','Semi-labile DOC')    

      ! Register links to particulate organic matter pools.
      call self%get_parameter(self%nRP,'nRP','','number of substrates',default=0)
      allocate(self%id_RP(self%nRP))
      allocate(self%id_RPc(self%nRP))
      allocate(self%id_RPn(self%nRP))
      allocate(self%id_RPp(self%nRP))
      allocate(self%id_RPf(self%nRP))
      do iRP=1,self%nRP
         write (index,'(i0)') iRP
         call self%register_state_dependency(self%id_RPc(iRP),'RP'//trim(index)//'c','mg C/m^3',   'POC '//trim(index))    
         call self%register_state_dependency(self%id_RPn(iRP),'RP'//trim(index)//'n','mmol N/m^3', 'PON '//trim(index))    
         call self%register_state_dependency(self%id_RPp(iRP),'RP'//trim(index)//'p','mmol P/m^3', 'POP '//trim(index))    
         call self%register_model_dependency(self%id_RP(iRP),'RP'//trim(index))
         call self%request_coupling_to_model(self%id_RPc(iRP),self%id_RP(iRP),'c')
         call self%request_coupling_to_model(self%id_RPn(iRP),self%id_RP(iRP),'n')
         call self%request_coupling_to_model(self%id_RPp(iRP),self%id_RP(iRP),'p')
#ifdef IRON
         call self%register_state_dependency(self%id_RPf(iRP),'RP'//trim(index)//'f','umol Fe/m^3','POFe '//trim(index),required=.false.)    
         call self%request_coupling_to_model(self%id_RPf(iRP),self%id_RP(iRP),'f')
#endif
      end do
      
#ifdef DOCDYN
      ! Register links to semi-refractory dissolved organic matter pool.
      call self%register_state_dependency(self%id_R3c,'R3c','mg C/m^3','Semi-refractory DOC')

      allocate(self%sRPR1(self%nRP))
      do iRP=1,self%nRP
         write (index,'(i0)') iRP
         call self%get_parameter(self%sRPR1(iRP),'sRP'//trim(index)//'R1','1/d','remineralisation of substrate '//trim(index)//' to DOM')
      end do

      call self%get_parameter(self%rR2B1X,'rR2','-','fraction of semi-labile DOC available to bacteria')
      call self%get_parameter(self%rR3B1X,'rR3','-','fraction of semi-refractory DOC available to bacteria')
      call self%get_parameter(self%frB1R3,'frR3','-','fraction of activity respiration converted to semi-refractory DOC')
#else
      allocate(self%puRP_B1X(self%nRP))
      do iRP=1,self%nRP
         write (index,'(i0)') iRP
         call self%get_parameter(self%puRP_B1X(iRP),'puRP'//trim(index),'-','fraction of substrate '//trim(index)//' available to bacteria')
      end do

      call self%get_parameter(self%R1R2X,'R1R2','-','labile fraction of produced DOM')
#endif

      ! Register links to external total dissolved inorganic carbon, dissolved oxygen pools
      call self%register_state_dependency(self%id_O3c,'O3c','mmol C/m^3','Carbon Dioxide')
      call self%register_state_dependency(self%id_O2o,'O2o','mmol O/m^3','Oxygen')

      ! Register environmental dependencies (temperature, suspendend sediment, pH, oxygen saturation)
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      call self%register_dependency(self%id_ESS,type_bulk_standard_variable(name='mass_concentration_of_silt'))
      if (self%ISWphx==1) call self%register_dependency(self%id_phx,standard_variables%ph_reported_on_total_scale)
      call self%register_dependency(self%id_eO2mO2,standard_variables%fractional_saturation_of_oxygen)

   end subroutine
   
   subroutine do(self,_ARGUMENTS_DO_)

      class (type_pml_ersem_bacteria),intent(in) :: self
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
#ifdef DOCDYN
      real(rk) :: R3c,R2cP,R3cP
      real(rk) :: fB1R1c
      real(rk) :: totsubst
      integer  :: iRP
      real(rk),dimension(self%nRP) :: RPc,RPcP,RPnP,RPpP
      real(rk),dimension(self%nRP) :: fRPB1p,fRPB1n
#else
      real(rk) :: Nlim,Plim
#endif

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
#ifdef DOCDYN
         _GET_WITH_BACKGROUND_(self%id_R3c,R3c)
         _GET_(self%id_R2c,R2cP)
         _GET_(self%id_R3c,R3cP)
         do iRP=1,self%nRP
            _GET_WITH_BACKGROUND_(self%id_RPc(iRP),RPc(iRP))
            _GET_(self%id_RPc(iRP),RPcP(iRP))
            _GET_(self%id_RPn(iRP),RPnP(iRP))
            _GET_(self%id_RPp(iRP),RPpP(iRP))
         end do
#endif
         qpB1c = B1p/B1c
         qnB1c = B1n/B1c
!..Temperature effect on pelagic bacteria:

         etB1 = self%q10B1X**((ETW-10._rk)/10._rk) - self%q10B1X**((ETW-32._rk)/3._rk)

!..Prevailing Oxygen limitation for bacteria:

         eO2B1 = eO2mO2/( eO2mO2 + self%chdB1oX )

!..Potential nutrient limitation on bacterial growth after Anderson

#ifndef DOCDYN
         IF (self%iswBlimX .EQ. 1 ) THEN
           Nlim = min(N4nP/(N4nP+self%chB1nX),R1nP/(R1nP+self%chB1nX))
           Plim = min(N1pP/(N1pP+self%chB1pX),R1pP/(R1pP+self%chB1pX))
         ELSE IF (self%iswBlimX .EQ. 2 ) THEN
           Nlim = (N4nP+R1nP)/(N4nP+R1nP+self%chB1nX)
           Plim = (N1pP+R1pP)/(N1pP+R1pP+self%chB1pX)
         END IF
#endif

!..bacterial mortality

         sB1RD=self%sdB1X*etB1
#ifdef DOCDYN
         fB1R1c = sB1RD*B1cP
#else
         fB1RDc = sB1RD*B1cP
#endif

!..Total amount of substrate available

         sutB1 = 1._rk  ! DOM-specific uptake rate
         ! rutB1 = sutB1*R1cP(I)

!..Potential uptake :
#ifdef DOCDYN
         rumB1 = self%sumB1X*etB1*eO2B1*B1c
#else
         rumB1 = self%sumB1X*etB1*eO2B1*min(Nlim,Plim)*B1c
#endif
!..Actual uptake

      ! rugB1 = MIN(rumB1,rutB1)
      ! specific in substrate concentration:
#ifdef DOCDYN
         totsubst = R1c+R2c*self%rR2B1X+R3c*self%rR3B1X+sum(RPc*self%sRPR1/sutB1)
         ! Jorn: check whether total substrate>0 to prevent NaNs
         if (totsubst>0.0_rk) then
            sugB1 = rumB1/max(rumB1/sutB1,totsubst)
         else
            sugB1 = 0.0_rk
         end if
             ! = MIN(rumB1,rutB1)=MIN(rumB1/(R1cP+R2cP*rR2B1X,sutB1) avoid pot. div. by 0
         rugB1 = sugB1*(R1cP+R2cP*self%rR2B1X+R3c*self%rR3B1X+sum(RPc*self%sRPR1/sutB1))
#else
         sugB1 = rumB1/max(rumB1/sutB1,R1c)
             ! = MIN(rumB1,rutB1)=MIN(rumB1/R1cP,sutB1) avoid pot. div. by 0
         rugB1 = sugB1*R1cP
#endif

!..Respiration :

! activity respiration
          rraB1 = rugB1 * ( 1._rk - self%puB1X*eO2mO2  - self%puB1oX*( 1._rk - eO2mO2 ) )
! total respiration:
          fB1O3c = rraB1 + self%srsB1X * B1cP * etB1
!      fB1O3c(I) = ( 1._fp8 - puB1X*eO2mO2  - puB1oX*( 1._fp8 - eO2mO2 ) )&
!                * rugB1  + srsB1X * B1cP(I) * etB1

#ifdef DOCDYN
! specific release of semilabile DOC 
! fudge factor 1. as used in Polimene et al. AME 2006
         sB1R2=max(0._rk,max(1._rk-(qpB1c/self%qpB1cX),1._rk-(qnB1c/self%qnB1cX)))*1._rk
         fB1R2c=sB1R2*B1cP
         fB1R3c=self%frB1R3*rraB1
         fB1RDc = fB1R1c + fB1R2c + fB1R3c
#endif

!..net bacterial production

         netb1 = rugB1 - fB1o3c - fB1RDc
         IF (netB1.gt.0._rk) THEN
            BGE=netB1/rugB1
         ELSE
            BGE=0._rk
         ENDIF

!..Source equations


         _SET_ODE_(self%id_c,netb1)
#ifdef DOCDYN
         _SET_ODE_(self%id_R1c,+ fB1R1c - sugB1*R1cP)
         _SET_ODE_(self%id_R2c,+ fB1R2c - sugB1*R2cP*self%rR2B1X)
         _SET_ODE_(self%id_R3c,+ fB1R3c - sugB1*R3cP*self%rR3B1X)
         do iRP=1,self%nRP
            _SET_ODE_(self%id_RPc(iRP),- sugB1*RPcP(iRP)*self%sRPR1(iRP))
         end do
#else
         _SET_ODE_(self%id_R1c,- rugB1 + (fB1RDc * self%R1R2X))
         _SET_ODE_(self%id_R2c,+ (fB1RDc * (1._rk - self%R1R2X)))
#endif

         _SET_ODE_(self%id_O3c,+ fB1O3c/CMass)
         _SET_ODE_(self%id_O2o,- fB1O3c*self%urB1_O2X)

!..Phosporous dynamics in bacteria........................................

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

#ifdef DOCDYN
         fRPB1p = sugB1*RPpP*self%sRPR1
         _SET_ODE_(self%id_p, sum(fRPB1p))
         do iRP=1,self%nRP
            _SET_ODE_(self%id_RPp(iRP), - fRPB1p(iRP))
         end do
#endif
         _SET_ODE_(self%id_p, + fR1B1p - fB1N1p - fB1RDp)
         _SET_ODE_(self%id_N1p, + fB1N1p)
         _SET_ODE_(self%id_R1p, + fB1RDp - fR1B1p)

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

#ifdef DOCDYN
         fRPB1n = sugB1*RPnP*self%sRPR1
         _SET_ODE_(self%id_n, sum(fRPB1n))
         do iRP=1,self%nRP
            _SET_ODE_(self%id_RPn(iRP), - fRPB1n(iRP))
         end do
#endif
         _SET_ODE_(self%id_N4n, + fB1NIn)
         _SET_ODE_(self%id_n,   + fR1B1n - fB1NIn - fB1RDn)
         _SET_ODE_(self%id_R1n, + fB1RDn   - fR1B1n)

      ! Leave spatial loops (if any)
      _LOOP_END_

      call remineralization(self,_ARGUMENTS_DO_)

   end subroutine do

   
   subroutine remineralization(self,_ARGUMENTS_DO_)

      class (type_pml_ersem_bacteria),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      
      integer :: iRP
      real(rk) :: ETW,ESS,phx
      real(rk) :: etB1
#ifndef DOCDYN
      real(rk) :: R2cP
      real(rk),dimension(self%nRP) :: RPc,RPn,RPp,RPcP,RPpP,RPnP
      real(rk),dimension(self%nRP) :: qnRPc,fRPr1
      real(rk) :: ruRPRD
#endif
#ifdef NOBAC
      real(rk) :: R1cP
      real(rk) :: fB1O3c
#endif
      real(rk) :: R1pP,R1nP,N4nP
#ifdef IRON
      real(rk) :: N7fP
#endif
      real(rk) :: fR1N1p,fR1NIn
      real(rk) :: sRPr1(self%nRP)
      real(rk) :: RPfP(self%nRP),fRPN7f(self%nRP),n7fsink
      
      real(rk) :: Fph,fN4N3n

      _LOOP_BEGIN_

         _GET_(self%id_ETW,ETW)
         _GET_(self%id_ESS,ESS)

#ifdef NOBAC
         _GET_(self%id_R1c,R1cP)
#endif
         _GET_(self%id_R1p,R1pP)
         _GET_(self%id_R1n,R1nP)
         RPfP = 0.0_rk
         do iRP=1,self%nRP
            if (_AVAILABLE_(self%id_RPf(iRP))) _GET_(self%id_RPf(iRP),RPfP(iRP))
         end do

         _GET_(self%id_N4n,N4nP)
#ifdef IRON
         _GET_(self%id_N7f,N7fP)
#endif

         etB1 = self%q10B1X**((ETW-10._rk)/10._rk) - self%q10B1X**((ETW-32._rk)/3._rk)

!..fluxes from particulate to dissolved

#ifdef DOCDYN
         sRPr1 = self%sRPr1
#else
         _GET_(self%id_R2c,R2cP)
         do iRP=1,self%nRP
            _GET_WITH_BACKGROUND_(self%id_RPc(iRP),RPc(iRP))
            _GET_WITH_BACKGROUND_(self%id_RPn(iRP),RPn(iRP))
            _GET_WITH_BACKGROUND_(self%id_RPp(iRP),RPp(iRP))
            _GET_(self%id_RPc(iRP),RPcP(iRP))
            _GET_(self%id_RPp(iRP),RPpP(iRP))
            _GET_(self%id_RPn(iRP),RPnP(iRP))
         end do
         qnRPc = RPn/RPc
         sRPr1  = self%redfieldX*CMass * qnRPc * self%puRP_B1X 
         fRPr1  = sRPr1 * RPcP
         ruRPRD = sum(fRPr1)
#ifdef NOBAC
         ruRPRD = ruRPRD*etB1
#endif
#ifdef SAVEFLX
         fRPRDc(I) = ruRPRD
#endif
#endif

#ifdef NOBAC
!..mineralisation of DOC to CO2

         fB1O3c = self%sR1N1X * R1cP*etB1

!..mineralisation of DOP to PO4

         fR1N1p = self%sR1N1X * R1pP*etB1

!..mineralisation of DON to Nh4

         fR1NIn = self%sR1N4X * R1nP*etB1

#else 
!..mineralisation of DOP to PO4

         fR1N1p = self%sR1N1X * R1pP

!..mineralisation of DON to Nh4

         fR1NIn = self%sR1N4X * R1nP
#endif

#ifdef IRON
! remineralization of particulate iron to Fe

         fRPN7f=sRPr1*RPfP

! sink of Fe

! This term takes into account the scavenging due to hydroxide precipitation and it is supposed to be 
! regulated by a threshold concentration (0.6 nM). See Aumont et al., 2003 (GBC) and Vichi et al., 2007 (JMS) for references.
! (Luca, 12/08)
         n7fsink=self%fsinkX*max(0._rk,N7fP-0.6_rk)
#endif

    !..Source equations

#ifdef DOCDYN
         _SET_ODE_(self%id_R1p, - fR1N1p)
         _SET_ODE_(self%id_R1n, - fR1NIn)
#else
#ifdef NOBAC
         _SET_ODE_(self%id_R1C,  + ruRPRD - fB1O3c)
         _SET_ODE_(self%id_O3c, + fB1O3c/CMass)
         _SET_ODE_(self%id_O2o, - fB1O3c*self%urB1_O2X)
#else
         _SET_ODE_(self%id_R1C, + ruRPRD)
#endif
         do iRP=1,self%nRP
            _SET_ODE_(self%id_RPc(iRP), - fRPr1(iRP))
            _SET_ODE_(self%id_RPp(iRP), - sRPr1(iRP)*RPpP(iRP))
            _SET_ODE_(self%id_RPn(iRP), - sRPr1(iRP)*RPnP(iRP))
         end do
         _SET_ODE_(self%id_R1p, + sum(sRPr1*RPpP) - fR1N1p)
         _SET_ODE_(self%id_R1n, + sum(sRPr1*RPnP) - fR1NIn)
#endif

         _SET_ODE_(self%id_N1p, + fR1N1p)
         _SET_ODE_(self%id_N4n, + fR1NIn)

#ifdef IRON
         do iRP=1,self%nRP
            if (fRPn7f(iRP)/=0.0_rk) _SET_ODE_(self%id_RPf(iRP),-fRPn7f(iRP))
         end do
         _SET_ODE_(self%id_N7f,+sum(fRPn7f)-n7fsink)
#endif

!..Nitrification..
!Ph influence on nitrification - empirical equation
! use(1) or not(2)
         if(self%ISWphx.eq.1)then
            _GET_(self%id_phx,phx)
            Fph = MIN(2._rk,MAX(0._rk,0.6111_rk*phx-3.8889_rk))
            fN4N3n = self%sN4N3X  * Fph * N4nP * etB1 * ESS / self%cessX
         else
            fN4N3n = self%sN4N3X  * N4nP * etB1 * ESS / self%cessX
         endif
         _SET_ODE_(self%id_N3n, + fN4N3n)
         _SET_ODE_(self%id_N4n, - fN4N3n)

!..Breakdown of semi labile to labile DOM

#ifndef DOCDYN
         _SET_ODE_(self%id_R1c, + self%rR2R1X * R2cP)
         _SET_ODE_(self%id_R2c, - self%rR2R1X * R2cP)
#endif

      ! Leave spatial loops (if any)
      _LOOP_END_

   end subroutine remineralization
   
end module