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
! Iron use is controlled by use_iron defined in ersem/shared.F90.
!
! CO2-enhanced photosynthesis is controlled by runtime switch cenh.
! -----------------------------------------------------------------------------

module ersem_primary_producer

   use fabm_types
   use fabm_particle

   use ersem_shared
   use ersem_pelagic_base

   implicit none

   private

   type,extends(type_ersem_pelagic_base),public :: type_ersem_primary_producer
      ! NB: own state variables (c,n,p,s,f,chl) are added implicitly by deriving from type_ersem_pelagic_base!

      ! Identifiers for state variables of other models
      type (type_state_variable_id) :: id_O3c,id_O2o,id_TA                  ! dissolved inorganic carbon, oxygen, total alkalinity
      type (type_state_variable_id) :: id_N1p,id_N3n,id_N4n,id_N5s,id_N7f   ! nutrients: phosphate, nitrate, ammonium, silicate, iron
      type (type_state_variable_id) :: id_R1c,id_R1p,id_R1n,id_R2c          ! dissolved organic carbon (R1: labile, R2: semi-labile)
      type (type_state_variable_id) :: id_RPc,id_RPp,id_RPn,id_RPs,id_RPf   ! particulate organic carbon
      type (type_state_variable_id) :: id_L2c                               ! Free calcite (liths) - used by calcifiers only

      ! Environmental dependencies
      type (type_dependency_id)            :: id_parEIR,id_ETW   ! PAR and temperature
      type (type_dependency_id)            :: id_RainR           ! Rain ratio - used by calcifiers only
      type (type_horizontal_dependency_id) :: id_pco2a3          ! Atmospheric pCO2 - used only if cenh is active

      ! Identifiers for diagnostic variables
      type (type_diagnostic_variable_id) :: id_fO3PIc  ! Gross primary production rate
      type (type_diagnostic_variable_id) :: id_fPIO3c  ! Respiration rate
      type (type_diagnostic_variable_id) :: id_fN3PIn,id_fN4PIn,id_fN1PIp,id_fN5PIs  ! nutrient uptake
      type (type_diagnostic_variable_id) :: id_netPI   ! Net primary production rate
      type (type_diagnostic_variable_id) :: id_lD      ! Cell-bound calcite - used by calcifiers only
      type (type_diagnostic_variable_id) :: id_O3L2c   ! Calcification
      type (type_diagnostic_variable_id) :: id_fPIRPc,id_fPIRPn,id_fPIRPp,id_fPIRPs  ! Total loss to Particulate detritus
      type (type_diagnostic_variable_id) :: id_fPIR1c,id_fPIR1n,id_fPIR1p ! Total loss to labile dissovled detritus
      type (type_diagnostic_variable_id) :: id_fPIR2c  ! Total loss to non-labile dissovled detritus
      type (type_diagnostic_variable_id) :: id_iNI

      ! Parameters (described in subroutine initialize, below)
      real(rk) :: sum
      real(rk) :: q10,srs,pu_ea,pu_ra,chs,qnlc,qplc,xqcp,exulim
      real(rk) :: xqcn,xqn,xqp,qun3,qun4,qurp,qsc,esni,snplux,xqcnx,xqcpx
      real(rk) :: resm,sdo
      real(rk) :: alpha,beta,phim
      real(rk) :: R1R2,uB1c_O2,urB1_O2
      real(rk) :: qflc,qfRc,qurf
      integer :: Limnut
      logical :: use_Si, calcify, docdyn, cenh

   contains

      ! Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: get_vertical_movement
      procedure :: get_sinking_rate

   end type type_ersem_primary_producer

   ! Constants
   real(rk),parameter :: ChlCmin = 0.002_rk

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_ersem_primary_producer),intent(inout),target :: self
      integer,                            intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
      real(rk) :: c0,EPS,iopABS,iopBBS
!EOP
!-----------------------------------------------------------------------
!BOC
      ! Obtain the values of all model parameters from FABM.
      ! Specify the long name and units of the parameters, which could be used by FABM (or its host)
      ! to present parameters to the user for configuration (e.g., through a GUI)
      call self%get_parameter(self%sum,   'sum',  '1/d',        'maximum specific productivity at reference temperature')
      call self%get_parameter(self%q10,   'q10',  '-',          'Q_10 temperature coefficient')
      call self%get_parameter(self%srs,   'srs',  '1/d',        'specific rest respiration at reference temperature')
      call self%get_parameter(self%pu_ea, 'pu_ea','-',          'excreted fraction of primary production')
      call self%get_parameter(self%exulim, 'exulim','-',          'excreted fraction of primary production due to nutrient stress',maximum=1._rk-self%pu_ea)
      call self%get_parameter(self%pu_ra, 'pu_ra','-',          'respired fraction of primary production')
      call self%get_parameter(self%qnlc,  'qnlc', 'mmol N/mg C','minimum nitrogen to carbon ratio')
      call self%get_parameter(self%qplc,  'qplc', 'mmol P/mg C','minimum phosphorus to carbon ratio')
      call self%get_parameter(self%xqcp,  'xqcp', '-',          'threshold for phosphorus limitation (relative to Redfield ratio)')
      call self%get_parameter(self%xqcn,  'xqcn', '-',          'threshold for nitrogen limitation (relative to Redfield ratio)')
      call self%get_parameter(self%xqcpx,  'xqcpx', '-','threshold for phosphorus limitation, relative to minimum ratio',default=2.0_rk )
      call self%get_parameter(self%xqcnx,  'xqcnx', '-','threshold for nitrogen limitation, relative to minimum ratio', default=2.0_rk )
      call self%get_parameter(self%xqp,   'xqp',  '-',          'maximum phosphorus to carbon ratio (relative to Redfield ratio)')
      call self%get_parameter(self%xqn,   'xqn',  '-',          'maximum nitrogen to carbon ratio (relative to Redfield ratio)')
      call self%get_parameter(self%qun3,  'qun3', 'm^3/mg C/d', 'nitrate affinity')
      call self%get_parameter(self%qun4,  'qun4', 'm^3/mg C/d', 'ammonium affinity')
      call self%get_parameter(self%qurp,  'qurp', 'm^3/mg C/d', 'phosphate affinity')
      call self%get_parameter(self%snplux,'snplux','1/d',       'specific tendency of luxury uptake of nutrients towards maximum quota',default=1.0_rk)
      call self%get_parameter(self%use_Si,'use_Si','',          'use silicate',default=.false.)
      if (self%use_Si) then
         call self%get_parameter(self%qsc,'qsc', 'mmol Si/mg C','maximum silicate to carbon ratio')
         call self%get_parameter(self%chs,'chs', 'mmol/m^3',    'Michaelis-Menten constant for silicate limitation')
      end if
      call self%get_parameter(self%sdo,   'sdo',  '1/d',        '1.1 of minimal specific lysis rate')
      call self%get_parameter(self%alpha, 'alpha','mg C m^2/mg Chl/W/d', 'initial slope of PI-curve')
      call self%get_parameter(self%beta,  'beta', 'mg C m^2/mg Chl/W/d','photoinhibition parameter')
      call self%get_parameter(self%phim,  'phim', 'mg Chl/mg C','maximum effective chlorophyll to carbon photosynthesis ratio')
      call self%get_parameter(self%Limnut,'Limnut','',          'nitrogen-phosphorus colimitation formulation (0: geometric mean, 1: minimum, 2: harmonic mean)',minimum=0,maximum=2)
      call self%get_parameter(self%docdyn,'docdyn','','use dynamic ratio of labile to semi-labile DOM production', default=.false.)
      if (.not.self%docdyn) call self%get_parameter(self%R1R2,'R1R2','-','labile fraction of produced dissolved organic carbon')
      call self%get_parameter(self%uB1c_O2,'uB1c_O2','mmol O_2/mg C','oxygen produced per unit of carbon fixed')
      call self%get_parameter(self%urB1_O2,'urB1_O2','mmol O_2/mg C','oxygen consumed per unit of carbon respired')
      if (use_iron) then
         call self%get_parameter(self%qflc,'qflc','umol Fe/mg C','minimum iron to carbon ratio')
         call self%get_parameter(self%qfRc,'qfRc','umol Fe/mg C','maximum/optimal iron to carbon ratio')
         call self%get_parameter(self%qurf,'qurf','m^3/mg C/d',  'specific affinity for iron')
      end if
      call self%get_parameter(EPS,     'EPS',    'm^2/mg C','specific shortwave attenuation', default=4.E-4_rk)
      call self%get_parameter(iopABS,  'iopABS', 'm^2/mg Chl','specific shortwave absorption', default=8.E-3_rk)
      call self%get_parameter(iopBBS,  'iopBBS', 'm^2/mg Chl','specific shortwave backscatter', default=3.E-3_rk)
      call self%get_parameter(c0,          'c0',     'mg C/m^3','background carbon concentration', default=0.0_rk)
      call self%get_parameter(self%calcify,'calcify','',        'calcify',                         default=.false.)
      call self%get_parameter(self%rm,     'rm',     'm/d',     'background sinking velocity',     default=0.0_rk)
      call self%get_parameter(self%resm,   'resm',   'm/d',     'maximum nutrient-limitation-induced sinking velocity', default=0.0_rk)
      call self%get_parameter(self%esni,   'esni',   '-',       'level of nutrient limitation below which sinking commences')
      call self%get_parameter(self%cenh,   'cenh',   '',        'enable atmospheric CO2 influence on photosynthesis', default=.false.)

      ! Register state variables (handled by type_ersem_pelagic_base)
      call self%initialize_ersem_base(sedimentation=.true.)
      call self%add_constituent('c',1.e-4_rk,   c0)
      call self%add_constituent('n',1.26e-6_rk, c0*qnrpicX)
      call self%add_constituent('p',4.288e-8_rk,c0*qprpicX)
      call self%add_constituent('f',5.e-6_rk,   0.0_rk)  ! NB this does nothing if iron support is disabled.
      call self%add_constituent('chl',3.e-6_rk, c0*self%phim)
      if (self%use_Si) call self%add_constituent('s',1.e-6_rk,c0*self%qsc)

      ! Register links to external nutrient pools.
      call self%register_state_dependency(self%id_N1p,'N1p','mmol P/m^3','phosphate')
      call self%register_state_dependency(self%id_N3n,'N3n','mmol N/m^3','nitrate')
      call self%register_state_dependency(self%id_N4n,'N4n','mmol N/m^3','ammonium')
      if (self%use_Si) call self%register_state_dependency(self%id_N5s,'N5s','mmol Si/m^3','silicate')
      if (use_iron)    call self%register_state_dependency(self%id_N7f,'N7f','umol Fe/m^3','inorganic iron')

      ! Register links to external labile dissolved organic matter pools (sink for excretion and lysis).
      call self%register_state_dependency(self%id_R1c,'R1c','mg C/m^3',  'dissolved organic carbon')
      call self%register_state_dependency(self%id_R1p,'R1p','mmol P/m^3','dissolved organic phosphorus')
      call self%register_state_dependency(self%id_R1n,'R1n','mmol N/m^3','dissolved organic nitrogen')

      ! Register links to external semi-labile dissolved organic matter pool (sink for excretion and lysis).
      call self%register_state_dependency(self%id_R2c,'R2c','mg C/m^3','semi-labile dissolved organic carbon')

      ! Register links to external particulate organic matter pools (sink for dead phytoplankton).
      ! At run-time, these can be coupled to any available pools (e.g., R4, R6, R8 in ERSEM)
      call self%register_state_dependency(self%id_RPc,'RPc','mg C/m^3',   'particulate organic carbon')
      call self%register_state_dependency(self%id_RPp,'RPp','mmol P/m^3', 'particulate organic phosphorus')
      call self%register_state_dependency(self%id_RPn,'RPn','mmol N/m^3', 'particulate organic nitrogen')
      if (self%use_Si) call self%register_state_dependency(self%id_RPs,'RPs','mmol Si/m^3','particulate organic silicate')
      if (use_iron)    call self%register_state_dependency(self%id_RPf,'RPf','umol Fe/m^3','particulate organic iron')

      ! Automatically hook up all components of external particulate organic matter,
      ! by obtaining them from a single named model "RP". This takes away the need to couple each RP?
      ! constituent individually.
      call self%request_coupling_to_model(self%id_RPc,'RP','c')
      call self%request_coupling_to_model(self%id_RPn,'RP','n')
      call self%request_coupling_to_model(self%id_RPp,'RP','p')
      if (_VARIABLE_REGISTERED_(self%id_RPs)) call self%request_coupling_to_model(self%id_RPs,'RP','s')
      if (_VARIABLE_REGISTERED_(self%id_RPf)) call self%request_coupling_to_model(self%id_RPf,'RP','f')

      ! Register links to external total dissolved inorganic carbon, dissolved oxygen pools
      call self%register_state_dependency(self%id_O2o,'O2o','mmol O_2/m^3','oxygen')
      call self%register_state_dependency(self%id_O3c,'O3c','mmol C/m^3','carbon dioxide')
      call self%register_state_dependency(self%id_TA,standard_variables%alkalinity_expressed_as_mole_equivalent)

      ! Register diagnostic variables (i.e., model outputs)
      call self%register_diagnostic_variable(self%id_netPI, 'netPI', 'mg C/m^3/d','net primary production')
      call self%register_diagnostic_variable(self%id_fO3PIc,'fO3PIc','mg C/m^3/d','gross primary production')
      call self%register_diagnostic_variable(self%id_fPIO3c,'fPIO3c','mg C/m^3/d','respiration')
      call self%register_diagnostic_variable(self%id_fN3PIn,'fN3PIn','mmol N/m^3/d','uptake of oxidised N')
      call self%register_diagnostic_variable(self%id_fN4PIn,'fN4PIn','mmol N/m^3/d','uptake of reduced N')
      call self%register_diagnostic_variable(self%id_fN1PIp,'fN1PIp','mmol P/m^3/d','uptake of P')
      if (self%use_Si) call self%register_diagnostic_variable(self%id_fN5PIs,'fN5PIs','mmol Si/m^3/d','uptake of Si')
      call self%register_diagnostic_variable(self%id_fPIR1c,'fPIR1c','mg C/m^3/d','loss to labile DOC')
      call self%register_diagnostic_variable(self%id_fPIR1n,'fPIR1n','mmol N/m^3/d','loss to labile DON')
      call self%register_diagnostic_variable(self%id_fPIR1p,'fPIR1p','mmol P/m^3/d','loss to labile DOP')
      call self%register_diagnostic_variable(self%id_fPIR2c,'fPIR2c','mg C/m^3/d','loss to semi-labile DOC')
      call self%register_diagnostic_variable(self%id_fPIRPc,'fPIRPc','mg C/m^3/d','loss to POC')
      call self%register_diagnostic_variable(self%id_fPIRPn,'fPIRPn','mmol N/m^3/d','loss to PON')
      call self%register_diagnostic_variable(self%id_fPIRPp,'fPIRPp','mmol P/m^3/d','loss to POP')
      call self%register_diagnostic_variable(self%id_iNI,'iNI','-','nutrient limitation factor')

      if (self%use_Si) call self%register_diagnostic_variable(self%id_fPIRPs,'fPIRPs','mmol Si/m^3/d','loss to POSi')

      ! Register environmental dependencies (temperature, shortwave radiation)
      call self%register_dependency(self%id_parEIR,standard_variables%downwelling_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      if (self%calcify) then
         ! Link to rain ratio (set by calcification module; excludes effect of nutrient limitation, which is computed below)
         call self%register_dependency(self%id_RainR,'RainR','-','rain ratio (PIC : POC)')

         ! Link to external detached liths (sink for calcite of dying phytoplankton)
         call self%register_state_dependency(self%id_L2c,'L2c','mg C/m^3','free calcite')

         ! Create diagnostic variable for cell-bound calcite, and register its contribution to the known quantity "total_calcite_in_biota".
         ! This quantity can be used by other models (e.g., predators) to determine how much calcite is released when phytoplankton is broken down.
         call self%register_diagnostic_variable(self%id_lD,'l','mg C/m^3','bound calcite',missing_value=0._rk,output=output_none)
         call self%register_diagnostic_variable(self%id_O3L2c,'calcification','mg C/m^3/d','calcification rate')
         call self%add_to_aggregate_variable(total_calcite_in_biota,self%id_lD)
      end if

      ! Link to atmospheric CO2 (only if using CO2-enhanced primary production).
      if (self%cenh) call self%register_dependency(self%id_pco2a3,standard_variables%mole_fraction_of_carbon_dioxide_in_air)

      ! Register contribution to light extinction
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &
         self%id_c,scale_factor=EPS,include_background=.true.)
      call self%add_to_aggregate_variable(particulate_organic_absorption_coefficient, &
         self%id_chl,scale_factor=iopABS,include_background=.true.)
      call self%add_to_aggregate_variable(particulate_organic_backscatter_coefficient, &
         self%id_chl,scale_factor=iopBBS,include_background=.true.)

   end subroutine

   subroutine do(self,_ARGUMENTS_DO_)

      class (type_ersem_primary_producer),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

   ! !LOCAL VARIABLES:
      real(rk) :: ETW,parEIR
      real(rk) :: c, p, n, Chl
      real(rk) :: cP,pP,nP,sP,ChlP
      real(rk) :: N5s,N1pP,N3nP,N4nP
      real(rk) :: iNn,iNp,iNs,iNf,iNI, iNnx, iNpx, iNIx
      real(rk) :: qpc,qnc

      real(rk) :: srs
      real(rk) :: fPIRPc,fPIRDc
      real(rk) :: fPIO3c,fO3PIc
      real(rk) :: fPIRPp,fPIRDp,fN1PIp
      real(rk) :: fPIRPn,fPIRDn
      real(rk) :: fNIPIn,fN3PIn,fN4PIn
      real(rk) :: fPIRPs,fPIN5s,fN5PIs

      real(rk) :: sdo,sum,sug,seo,sea,sra,sun,run,sPIRP
      real(rk) :: runp,misp,rump
      real(rk) :: runn,misn,rumn,rumn3,rumn4
      real(rk) :: et,pe_RP
      real(rk) :: rho,Chl_inc,Chl_loss
      real(rk) :: ChlCpp
      real(rk) :: N7fP,f,fP,qfc
      real(rk) :: runf,rumf,misf
      real(rk) :: fN7PIf,fPIRPf
      real(rk) :: fPIR1c,fPIR2c
      real(rk) :: pco2a3,cenh
      real(rk) :: RainR, t

      real(rk),parameter :: onedayX = 1.0_rk

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve local biomass (carbon, phosphorus, nitrogen, chlorophyll concentrations).

         ! Concentrations including background (used in source terms)
         _GET_WITH_BACKGROUND_(self%id_c,c)
         _GET_WITH_BACKGROUND_(self%id_p,p)
         _GET_WITH_BACKGROUND_(self%id_n,n)
         _GET_WITH_BACKGROUND_(self%id_chl,Chl)

         ! Concentrations excluding background (used in sink terms)
         _GET_(self%id_c,cP)
         _GET_(self%id_p,pP)
         _GET_(self%id_n,nP)
         _GET_(self%id_chl,ChlP)

         ! Retrieve ambient nutrient concentrations
         _GET_(self%id_N1p,N1pP)
         _GET_(self%id_N3n,N3nP)
         _GET_(self%id_N4n,N4nP)

         ! Retrieve environmental dependencies (water temperature, photosynthetically active radation)
         _GET_(self%id_ETW,ETW)
         _GET_(self%id_parEIR,parEIR)

         ! Ratios: phosphorus to carbon, nitrogen to carbon, chlorophyll to carbon.
         ! Note: these are protected against division by zero because c includes the background concentration.
         qpc = p/c
         qnc = n/c
         ChlCpp = Chl/c

         ! Regulation factors...................................................

         ! Nitrogen and phosphorus limitation factors based on internal quota
         iNp = MIN(1._rk,  &
                 MAX(0._rk, (qpc-self%qplc) / (self%xqcp*qpRPIcX-self%qplc) ))
         iNn = MIN(1._rk,  &
                 MAX(0._rk, (qnc-self%qnlc) / (self%xqcn*qnRPIcX-self%qnlc) ))

         if (use_iron) then
            ! Limitation factor based on internal iron quota
            _GET_WITH_BACKGROUND_(self%id_f,f)
            qfc = f/c
            iNf = MIN(1._rk,  &
                    MAX(ZeroX, (qfc-self%qflc) / (self%qfRc-self%qflc) ))
         else
            ! No iron limitation
            iNf = 1.0_rk
         end if

         if (self%use_Si) then
            ! Limitation factor based on ambient silicate
            if (legacy_ersem_compatibility) then
               ! Legacy ERSEM includes background value, but this is inappropriate as it is used in a sink term.
               _GET_WITH_BACKGROUND_(self%id_N5s,N5s)
            else
               _GET_(self%id_N5s,N5s)
            end if
            iNs = MIN(1._rk, N5s/(N5s+self%chs))
         else
            ! No silicate limitation
            iNs = 1.0_rk
         end if

         ! Co-limitation of inorganic nitrogen and phosphorus
         if (self%Limnut==0) then  ! Jorn: select case would be cleaner but makes vectorization impossible for ifort 14
            ! Geometric mean
            iNI = sqrt(iNp * iNn)
         elseif (self%Limnut==1) then
            ! Minimum
            iNI = min(iNp, iNn)
         else
            ! Harmonic mean
            iNI = 2.0_rk / (1._rk/iNp + 1._rk/iNn)
         end if


         _SET_DIAGNOSTIC_(self%id_iNI,iNI)


! Limitation factor for C uptake - as approach minimum C:P/C:N
         iNpx = MIN(1._rk,  &
                 MAX(0._rk, (qpc-self%qplc) / (self%xqcpx*self%qplc-self%qplc) ))
         iNnx = MIN(1._rk,  &
                 MAX(0._rk, (qnc-self%qnlc) / (self%xqcnx*self%qnlc-self%qnlc) ))

         iNIx= min(iNpx, iNnx)

         ! Temperature response
         et = max(0.0_rk,self%q10**((ETW-10._rk)/10._rk) - self%q10**((ETW-32._rk)/3._rk))

         ! Production...........................................................

         ! Gross photosynthetic activity (1/d), limited by availability of silicate and iron,
         ! but not nitrogen and phosphorus.
         sum = self%sum*et*iNs*iNf

         if (parEIR>zeroX) then
            sum = iNIx * sum * (1._rk-exp(-self%alpha*parEIR*ChlCpp/sum)) * exp(-self%beta*parEIR*ChlCpp/sum)
            rho = (self%phim - ChlCmin) * (sum/(self%alpha*parEIR*ChlCpp)) + ChlCmin
         else
            sum = 0._rk
            rho = ChlCmin
         end if

         if (self%cenh) then
            ! Retrieve atmospheric pCO2 if using CO2-enhanced primary production.
            _GET_HORIZONTAL_(self%id_pco2a3,pco2a3)

            ! Enhancement factor (from MEECE D1.5) 379.48 = pco2a @ 2005
            cenh=1.0_rk+(pco2a3-379.48_rk)*0.0005_rk
         else
            ! No influence of atmospheric CO2 on photosynthesis
            cenh = 1.0_rk
         end if

         ! Nutrient-stress lysis rate (1/d)
         sdo = (1._rk/(MIN(iNs, iNI)+0.1_rk))*self%sdo

         ! Excretion rate, as regulated by nutrient-stress (1/d)
     !    seo = sum*(1._rk-iNI)*(1._rk-self%pu_ea) - old ERSEM
         seo = sum*(1._rk-iNI)*self%exulim


         ! Activity-dependent excretion (1/d)
         sea = sum*self%pu_ea

         ! Productivity after subtracting excretion fluxes (1/d)
         sug = sum-seo-sea

         ! Fraction of lysis flux that is particulate (remainder is dissolved)
         pe_RP = MIN(self%qplc/(qpc+ZeroX), self%qnlc/(qnc+ZeroX), 1.0_rk)

         ! Specific loss (1/d) and carbon loss (mg C/m3/d) to particulate matter due to lysis
         sPIRP = pe_RP*sdo
         fPIRPc = sPIRP*cP

         if (self%docdyn) then
            ! Lysis produces labile DOM, excretion produces semi-labile DOM.
            fPIR1c = (1._rk-pe_RP)*sdo*cP
            fPIR2c = (seo + sea)*cP
            fPIRDc = fPIR1c+fPIR2c
         else
            ! Fixed ratio between production of labile and semi-labile DOM.
            fPIRDc = (1._rk-pe_RP)*sdo*cP + (seo + sea)*cP
            fPIR1c = fPIRDc*self%R1R2
            fPIR2c = fPIRDc*(1._rk-self%R1R2)
         end if

         ! Calcified matter:
         if (self%calcify) then
            _GET_(self%id_RainR,RainR)

            ! Calculate nutrient limitation impact on rain ratio:
            t=max(0._rk,ETW)  ! this is to avoid funny values of rain ratio when ETW ~ -2 degrees
            RainR = RainR * min((1._rk-iNp),iNn) * (t/(2._rk+t)) !* max(1.,P2c(I)/2.) removd as P2 is a broad class not just calicifiers
            RainR = max(RainR,0.005_rk)

            ! Compute virtual calcite attached to live cells. It is virtual in the sense that it has not been subtracted
            ! from the DIC budget yet. Thus, its carbon is not yet allocated. When consumed by predators, this is (partially)
            ! transformed into free liths. At that time, carbon used to form liths is subtracted from the DIC pool.
            _SET_DIAGNOSTIC_(self%id_lD,RainR*cP)

            ! For dying cells: convert virtual cell-attached calcite to actual calcite in free liths.
            ! Update DIC and lith balances accordingly (only now take away the DIC needed to form the calcite)
            _SET_ODE_(self%id_L2c,   fPIRPc*RainR)
            _SET_DIAGNOSTIC_(self%id_O3L2c,fPIRPc*RainR)
            _SET_ODE_(self%id_O3c,  -fPIRPc*RainR/Cmass)
            _SET_ODE_(self%id_TA, -2*fPIRPc*RainR/Cmass)   ! CaCO3 formation decreases alkalinity by 2 units
         end if

         ! Respiration..........................................................

         ! Rest respiration rate (1/d)
         srs = et*self%srs

         ! Activity respiration rate (1/d)
         sra = sug*self%pu_ra

         ! Total respiration = production of CO2 (mg C/m3/d)
         fPIO3c = srs*cP + sra*c*cenh

         ! Gross production = uptake of CO2 (mg C/m3/d)
         fO3PIc = sum*c*cenh

         ! Productivity and production (as used in nutrient uptake equations)
         sun = sum-(seo+sea+sra)  ! net productivity (1/d) - excludes excretion and activity respition, but not rest respiration
         run = sun*c-srs*cP       ! net production (mg C/m3/d) - excludes excretion and respiration

         ! Save net production (equals run if cenh=1)
         _SET_DIAGNOSTIC_(self%id_netPI,(sum*cenh-seo-sea-sra*cenh)*c-srs*cP)

         ! Chl changes (note that Chl is a component of PXc and not involved
         ! in mass balance)
         if (use_iron) then
           Chl_inc = min(iNf,iNI)*rho*(sum-sra-seo-sea)*c
         else
           Chl_inc = iNI*rho*(sum-sra-seo-sea)*c
         endif
         Chl_loss = (sdo+srs)*ChlP

         _SET_ODE_(self%id_c,(fO3PIc-fPIO3c-fPIRPc-fPIRDc))

         _SET_ODE_(self%id_R1c,fPIR1c)
         _SET_ODE_(self%id_R2c,fPIR2c)
         _SET_ODE_(self%id_RPc,fPIRPc)
         _SET_ODE_(self%id_chl,(Chl_inc - Chl_loss))

         _SET_ODE_(self%id_O3c,(fPIO3c - fO3PIc)/CMass)
         _SET_ODE_(self%id_O2o,(fO3PIc*self%uB1c_O2 - fPIO3c*self%urB1_O2))

         ! Save rates of photosynthesis (a.k.a., gross primary production) and respiration
         _SET_DIAGNOSTIC_(self%id_fPIO3c,fPIO3c)
         _SET_DIAGNOSTIC_(self%id_fO3PIc,fO3PIc)

         _SET_DIAGNOSTIC_(self%id_fPIR1c,fPIR1c)
         _SET_DIAGNOSTIC_(self%id_fPIR2c,fPIR2c)
         _SET_DIAGNOSTIC_(self%id_fPIRPc,fPIRPc)
         ! Phosphorus flux...........................................

         ! Lysis loss of phosphorus (mmol P m-3 d-1)
         fPIRPp = sPIRP * min(self%qplc*cP,pP)
         fPIRDp = sdo * pP - fPIRPp

         ! Net phosphorus uptake
         ! Maximum achievable uptake (mmol P m-3 d-1)
         rump = self%qurp * N1pP * c
 
         !Regulation term relaxing internal quota towards maximum quota (using nutrient luxury uptake)  (mmol P m-3 d-1)
         misp = self%snplux*(self%xqp * qpRPIcX*cP - pP)

         !Assimilation demand at maximum quota and compensating for rest respiration (mmol P m-3 d-1)
         runp = sun*c * qpRPIcX*self%xqp - srs*pP

         ! Uptake capped at maximum achievable uptake (mmol P m-3 d-1)
         fN1PIp = MIN(rump, runp+misp)

         ! Source equations
         _SET_ODE_(self%id_p,(fN1PIp-fPIRDp-fPIRPp))
         _SET_ODE_(self%id_N1p,-fN1PIp)
         _SET_ODE_(self%id_TA,  fN1PIp) ! Alkalinity contributions: -1 for PO4
         _SET_ODE_(self%id_RPp,fPIRPp)
         _SET_ODE_(self%id_R1p,fPIRDp)

         ! Set Diagnostic
         _SET_DIAGNOSTIC_(self%id_fN1PIp, fN1PIp)
         _SET_DIAGNOSTIC_(self%id_fPIR1p,fPIRDp)
         _SET_DIAGNOSTIC_(self%id_fPIRPp,fPIRPp)

         ! Nitrogen flux.............................................

         ! Nitrogen loss by lysis
         fPIRPn = sPIRP * min(self%qnlc*cP,nP)
         fPIRDn = sdo * nP - fPIRPn

         ! Net nitrogen uptake

         ! maximum acheivable uptake of nitrate  (mmol N m-3 d-1)
         rumn3 = self%qun3 * N3nP * c

         ! Maximum achievable uptake of ammonium (mmol N m-3 d-1)
         rumn4 = self%qun4 * N4nP * c

         !Total maximum achievable uptake of nitrogen (mmol N m-3 d-1)
         rumn = rumn3 + rumn4

         !Regulation term relaxing internal quota towards maximum quota (using nutrient luxury uptake) (mmol N m-3 d-1)
         misn = self%snplux * (self%xqn * qnRPIcX*cP - nP)

         !Assimilation demand at maximum quota and compensating for rest respiration (mmol N m-3 d-1)
         runn = sun*c * qnRPIcX*self%xqn - srs*nP

         ! Uptake capped at maximum achievable uptake of nitrogen (mmol N m-3 d-1)
         fNIPIn = MIN(rumn, runn + misn)

         ! Partitioning over NH4 and NO3 uptake
         IF (fNIPIn .gt. 0._rk) THEN
            fN3PIn = fNIPIn * rumn3 / rumn
            fN4PIn = fNIPIn * rumn4 / rumn
         ELSE
            fN3PIn = 0._rk
            fN4PIn = fNIPIn
         ENDIF

         ! Source equations
         _SET_ODE_(self%id_n,(fN4PIn+fN3PIn-fPIRDn-fPIRPn))
         _SET_ODE_(self%id_N3n,-fN3PIn)
         _SET_ODE_(self%id_N4n,-fN4PIn)
         _SET_ODE_(self%id_TA,  fN3PIn-fN4PIn) ! Alkalinity contributions: -1 for NO3, +1 for NH4
         _SET_ODE_(self%id_RPn,fPIRPn)
         _SET_ODE_(self%id_R1n,fPIRDn)

         ! Set Diagnostic
         _SET_DIAGNOSTIC_(self%id_fN3PIn, fN3PIn)
         _SET_DIAGNOSTIC_(self%id_fN4PIn, fN4PIn)
         _SET_DIAGNOSTIC_(self%id_fPIR1n,fPIRDn)
         _SET_DIAGNOSTIC_(self%id_fPIRPn,fPIRPn)

         if (self%use_Si) then
            ! Silicate flux.............................................
            _GET_(self%id_s,sP)

            ! Loss of silicate due to lysis (mmol Si m-3 d-1)
            fPIRPs = sdo * sP

            ! Loss of excess silicate (qsP1c > qsc)  (mmol Si m-3 d-1)
            fPIN5s = MAX ( 0._rk, sP-self%qsc * cP) / onedayX

            ! Net silicate uptake  (mmol Si m-3)
            fN5PIs = MAX ( 0._rk, self%qsc*run) - fPIN5s
            _SET_DIAGNOSTIC_(self%id_fN5PIs, fN5PIs)

            ! Source equations  
            _SET_ODE_(self%id_s,(fN5PIs - fPIRPs))
            _SET_ODE_(self%id_N5s,-fN5PIs)
            _SET_ODE_(self%id_RPs, fPIRPs)

            _SET_DIAGNOSTIC_(self%id_fPIRPs,fPIRPs)
         end if

         if (use_iron) then
            ! Iron flux................................................

            ! Obtain internal iron concentration (umol/m^3, excludes background), and external biologically available iron.
            _GET_(self%id_f,fP)
            _GET_(self%id_N7f,N7fP)

            ! Iron loss by lysis
            !  Because its high affinity with particles all the iron lost from phytoplankton by lysis is supposed to be
            !  associated to organic particulate detritus. (luca)

            fPIRPf = sdo * fP

            ! Net iron uptake
            rumf = self%qurf * N7fP * c
            misf = self%qfRc*cP - fP
            runf = sun*c * self%qfRc - srs*fP
            fN7PIf = MIN(rumf, runf+misf)

            ! Source equations
            _SET_ODE_(self%id_f,(fN7PIf-fPIRPf))
            _SET_ODE_(self%id_N7f,-fN7PIf)
            _SET_ODE_(self%id_RPf,fPIRPf)
         end if

      ! Leave spatial loops (if any)
      _LOOP_END_

   end subroutine do

   function get_sinking_rate(self,_ARGUMENTS_LOCAL_) result(SD)
      class (type_ersem_primary_producer),intent(in) :: self
      _DECLARE_ARGUMENTS_LOCAL_

      real(rk) :: c,p,n
      real(rk) :: qpc,qnc
      real(rk) :: SD
      real(rk) :: iNp,iNn,iNI

      _GET_WITH_BACKGROUND_(self%id_c,c)
      _GET_WITH_BACKGROUND_(self%id_p,p)
      _GET_WITH_BACKGROUND_(self%id_n,n)

      qpc = p/c
      qnc = n/c

      iNp = MIN(1._rk,  &
               MAX(0._rk, (qpc-self%qplc) / (self%xqcp*qpRPIcX-self%qplc) ))
      iNn = MIN(1._rk,  &
               MAX(0._rk, (qnc-self%qnlc) / (self%xqcn*qnRPIcX-self%qnlc) ))

      if (self%Limnut==0) then  ! NB select case would be cleaner but makes vectorization impossible for ifort 14
         iNI = sqrt(iNp * iNn)
      elseif (self%Limnut==1) then
         iNI = min(iNp, iNn)
      else
         iNI = 2.0_rk / (1._rk/iNp + 1._rk/iNn)
      end if

      ! Sedimentation and resting stages.....................................
      SD = self%resm * MAX(0._rk, self%esni - iNI) + self%rm
   end function

   subroutine get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)
      class (type_ersem_primary_producer),intent(in) :: self
      _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

      real(rk) :: SD

      _LOOP_BEGIN_

         ! Retrieve local state
         SD = -get_sinking_rate(self,_ARGUMENTS_LOCAL_)
         _SET_VERTICAL_MOVEMENT_(self%id_c,SD)
         _SET_VERTICAL_MOVEMENT_(self%id_p,SD)
         _SET_VERTICAL_MOVEMENT_(self%id_n,SD)
         if (self%use_Si) _SET_VERTICAL_MOVEMENT_(self%id_s,SD)
         _SET_VERTICAL_MOVEMENT_(self%id_chl,SD)
         if (use_iron) _SET_VERTICAL_MOVEMENT_(self%id_f,SD)

      _LOOP_END_

   end subroutine get_vertical_movement

end module
