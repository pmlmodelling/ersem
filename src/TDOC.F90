#include "fabm_driver.h"

module ersem_TDOC

   use fabm_types
   use fabm_builtin_models
   use fabm_particle

   use ersem_shared
   use ersem_pelagic_base

   implicit none

   private

   type,extends(type_ersem_pelagic_base),public :: type_ersem_TDOC
      ! Variables
      type (type_state_variable_id) :: id_RPc,id_T2c, id_O3c
      type (type_state_variable_id) :: id_RPn, id_RPp, id_T2n, id_T2p, id_N1p, id_N4n
      type (type_state_variable_id) :: id_T1_older_c,id_T1_older_n,id_T1_older_p
      type (type_model_id) :: id_T2, id_T1_older
      type (type_dependency_id)                     :: id_ETW, id_chemEIR
      type (type_dependency_id)            :: id_X1X             ! Salinity
!      type (type_horizontal_dependency_id)          :: id_R1c,id_R1d
      type (type_horizontal_diagnostic_variable_id) :: id_surface_photolysis, id_surface_photo_aging
      type (type_diagnostic_variable_id) :: id_photolysis,id_flocc, id_photo_aging

   ! Parameters

      real(rk) :: suva, iref,phyref,surf_phyref,phyt,floc,qp,qn,scx,sbx,chemEIR_scaling
      real(rk) :: age, photoaging,bioaging
      
   contains
!     Model procedures
      procedure :: initialize
      procedure :: do_surface
      procedure :: do
   end type
   
   ! Submodel that is seen by bacteria for microbial degradation
   ! this is currently needed because the microbial degradation flux (calculated by the bacteria module) is needed here for the aging of tDOC
   type,extends(type_particle_model),public :: type_ersem_shadow_tDOC
      type (type_diagnostic_variable_id) :: id_c_shadow   ! shadow state variable
      type (type_state_variable_id) :: id_parent_c, id_T1_older_parent_c ! ,id_parent_n, id_T1_older_parent_n,id_parent_p, id_T1_older_parent_p
      type (type_model_id) :: id_parent, id_T1_older_parent
      
      real (rk)  :: qn, qp, bioaging

   contains
      procedure :: initialize => shadow_initialize
      procedure :: do  => shadow_do
   end type


contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_ersem_TDOC),intent(inout),target :: self
      integer,                         intent(in)           :: configunit

! !LOCAL VARIABLES:
      real(rk) :: c0,EPS
      class (type_ersem_shadow_tDOC), pointer :: shadow      
!
!EOP
!-----------------------------------------------------------------------
!BOC
!      call self%get_parameter(self%kuw,'kuw', 'm-1',              'extinction of UV in the water')
!      call self%get_parameter(self%suva,'suva','m2*mmol C-1',            'specific UV absorption at 350 nm')
      call self%get_parameter(self%iref,'iref','W m-2',            'reference irradiance')
      call self%get_parameter(self%phyref,'phyref','d-1 ',            'reference photooxidation rate')
      call self%get_parameter(self%surf_phyref,'surf_phyref','d-1 ',            'reference surface_photooxidation rate',default=0._rk)
      call self%get_parameter(self%phyt,'phyt','adim',            'photooxidated fraction of T1 going into T2) ')
      call self%get_parameter(self%floc,'floc','(mmol C-3)-1*d-1 ',            'reference photooxidation rate')
      call self%get_parameter(self%qp,    'qp',    'mmol P/mg C','phosphorus to carbon ratio')
      call self%get_parameter(self%qn,    'qn',    'mmol N/mg C','nitrogen to carbon ratio')
      call self%get_parameter(EPS,     'EPS',    'm^2/mg C','specific shortwave attenuation', default=4.E-4_rk)
      call self%get_parameter(c0,'c0','mg C/m^3','background carbon concentration')
      call self%get_parameter(self%sbx,   'sbx',   'psu',        'optimal salinity')
      call self%get_parameter(self%scx,   'scx',   '',        'salinity function parameter')
      call self%get_parameter(self%chemEIR_scaling, 'chemEIR_scaling', '', 'scaling factor for incoming radiation activating photochemical reaction', default=1._rk)
      !call self%get_parameter(self%age, 'age', 'd', 'Characteristic age of the tDOC pool', default=0._rk))  ! for now the parameter is calculated a priori
      call self%get_parameter(self%photoaging, 'photoaging', '-', 'aging due to photo oxidation')
      call self%get_parameter(self%bioaging, 'bioaging', '-', 'aging due to microbial oxidation')
      
! Allow ERSEM base model to declare our own state variables.
      call self%initialize_ersem_base(sedimentation=.false.)
      call self%add_constituent('c',1.e-4_rk,c0,qn=self%qn,qp=self%qp)
!      call self%add_constituent('c',0.0_rk)

! Register links to nutrient pools.
!      call self%register_dependency(self%id_EIR,'EIR','W/m^2','downwelling_shortwave_flux', &
!              standard_variable=standard_variables%downwelling_shortwave_flux,source=source_do_column)
      !call self%register_dependency(self%id_EIR,standard_variables%downwelling_shortwave_flux)
      call self%register_dependency(self%id_chemEIR,'chemEIR','W/m^2','incoming radiation activating photochemical reaction')
      call self%register_dependency(self%id_X1X,standard_variables%practical_salinity)
      call self%register_state_dependency(self%id_RPc,'RPc','mg C/m^3',   'particulate organic carbon')
      call self%register_state_dependency(self%id_RPp,'RPp','mmol P/m^3', 'particulate organic phosphorus')
      call self%register_state_dependency(self%id_RPn,'RPn','mmol N/m^3', 'particulate organic nitrogen')
      call self%register_model_dependency(self%id_T1_older,'T1_older')
      call self%register_state_dependency(self%id_T1_older_c,'T1_older_c','mg C/m^3','non photolabile terrigenous DOC')
      call self%register_state_dependency(self%id_T1_older_n,'T1_older_n','mmol N/m^3','non photolabile terrigenous DON')
      call self%register_state_dependency(self%id_T1_older_p,'T1_older_p','mmol P/m^3','non photolabile terrigenous DOP')
      call self%request_coupling_to_model(self%id_T1_older_c,self%id_T1_older,'c')
      call self%request_coupling_to_model(self%id_T1_older_n,self%id_T1_older,standard_variables%total_nitrogen)
      call self%request_coupling_to_model(self%id_T1_older_p,self%id_T1_older,standard_variables%total_phosphorus)

      call self%register_model_dependency(self%id_T2,'T2')
      call self%register_state_dependency(self%id_T2c,'T2c','mg C/m^3','non photolabile terrigenous DOC')
      call self%register_state_dependency(self%id_T2n,'T2n','mmol N/m^3','non photolabile terrigenous DON')
      call self%register_state_dependency(self%id_T2p,'T2p','mmol P/m^3','non photolabile terrigenous DOP')
      call self%request_coupling_to_model(self%id_T2c,self%id_T2,'c')
      call self%request_coupling_to_model(self%id_T2n,self%id_T2,standard_variables%total_nitrogen)
      call self%request_coupling_to_model(self%id_T2p,self%id_T2,standard_variables%total_phosphorus)

      call self%register_state_dependency(self%id_O3c,'O3c','mmol C/m^3','carbon dioxide sink')
      ! Register links to nutrient pools.
      call self%register_state_dependency(self%id_N1p,'N1p','mmol P/m^3','phosphate')
      call self%register_state_dependency(self%id_N4n,'N4n','mmol N/m^3','ammonium')
      call self%register_diagnostic_variable(self%id_photolysis,'photolysis','mgC/m^3/d','photolysis')
      call self%register_diagnostic_variable(self%id_photo_aging,'photoaging','mgC/m^3/d','aging due to 3D photoloysis')
      call self%register_horizontal_diagnostic_variable(self%id_surface_photolysis,'surface_photolysis','mgC/m^2/d','surface photolysis')
      call self%register_horizontal_diagnostic_variable(self%id_surface_photo_aging,'surface_photoaging','mgC/m^2/d','aging due to surface photolysis')
      call self%register_diagnostic_variable(self%id_flocc,'flocc','mgC/m^3/d','flocculation')

      ! Register contribution to light extinction
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &
         self%id_c,scale_factor=EPS,include_background=.true.)
         

      allocate(shadow)
      shadow%qn=self%qn
      shadow%qp=self%qp
      shadow%bioaging=self%bioaging
      call shadow%couplings%set_string('parent',self%name)
      call shadow%couplings%set_string('T1_older_parent','T1_older')
      call self%add_child(shadow,'shadow',configunit=configunit)      

   end subroutine initialize
   
   
   subroutine do(self,_ARGUMENTS_DO_)

      class (type_ersem_TDOC),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

   ! !LOCAL VARIABLES:
    real(rk) :: flocc,photolysis,R8cP,R8c,T2c,T2n,T2p,O3c,O3cP,T2cP,chemEIR,c,Xp,Xn,XXn,T1T2,px,nx
    real(rk) :: X1X,sal

! Enter spatial loops (if any)
      _LOOP_BEGIN_

      _GET_(self%id_c,c)
      _GET_(self%id_T2c,T2cP)
!      _GET_(self%id_R8c,R8cP)
      _GET_(self%id_O3c,O3cP)
      _GET_WITH_BACKGROUND_(self%id_T2c,T2c)
      _GET_WITH_BACKGROUND_(self%id_T2n,T2n)
      _GET_WITH_BACKGROUND_(self%id_T2p,T2p)
!      _GET_WITH_BACKGROUND_(self%id_R8c,R8c)
      _GET_(self%id_chemEIR,chemEIR)
      _GET_(self%id_X1X,X1X)

    IF(T2c.gt.0._rk) then
    Xn=self%qn/(T2n/T2c)
    Xp=self%qp/(T2p/T2c)
    nx=self%qn-(T2n/T2c)
    px=self%qp-(T2p/T2c)
    else
    Xn=0._rk
    Xp=0._rk
    nx=1._rk
    px=1._rk
    endif

    XXn=min(Xn,Xp)

!Flocculation is assumed to be (log-normal) function of salinity(sal)

     X1X=max(X1X,0.01_rk) !to avoid log(0) in sal!!
     sal=exp(-((log(X1X)-self%sbx)**2._rk)/(2._rk*self%scx**2._rk))
     flocc=self%floc*sal*c**2

! Photolysis
    photolysis=self%phyref*(chemEIR*self%chemEIR_scaling/self%iref)*c

! in order to ensure mass conservation, the fraction of carbon going from T1 to T2 after photolysis is down regulated 
! if T1 has less  N and/or P than T2. Any excess of nutrients is redirected into the dissolved pool (N1p and N4n). (Luca, July 2018)    
    T1T2=self%phyt*min(1._rk,XXn)

    _SET_ODE_(self%id_c,-(1._rk+self%photoaging)*photolysis-flocc)
    
    _SET_ODE_(self%id_T1_older_c,self%photoaging*photolysis)
    _SET_ODE_(self%id_T1_older_n,self%photoaging*photolysis*self%qn)
    _SET_ODE_(self%id_T1_older_p,self%photoaging*photolysis*self%qp)
    
    _SET_ODE_(self%id_RPc,+flocc)
    _SET_ODE_(self%id_RPn,+flocc*self%qn)
    _SET_ODE_(self%id_RPp,+flocc*self%qp)
    _SET_ODE_(self%id_T2c,+photolysis*T1T2)
    _SET_ODE_(self%id_O3c,+photolysis*(1._rk-T1T2)/CMass)
    _SET_ODE_(self%id_N1p,+photolysis*(1._rk-T1T2)*self%qp+photolysis*T1T2*px)
    _SET_ODE_(self%id_N4n,+photolysis*(1._rk-T1T2)*self%qn+photolysis*T1T2*nx)
    _SET_DIAGNOSTIC_(self%id_photolysis, photolysis)
    _SET_DIAGNOSTIC_(self%id_flocc, flocc)

    _LOOP_END_

   end subroutine do

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_ersem_TDOC), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

         real (rk) :: chemEIR,photolysis,c,T1T2,T2c,T2n,T2p
         real (rk) :: Xn,Xp,XXn,nx,px

! Enter horizontal loops (if any)
      _HORIZONTAL_LOOP_BEGIN_

         _GET_(self%id_chemEIR,chemEIR)
         _GET_(self%id_c,c)
         _GET_WITH_BACKGROUND_(self%id_T2c,T2c)
         _GET_WITH_BACKGROUND_(self%id_T2n,T2n)
         _GET_WITH_BACKGROUND_(self%id_T2p,T2p)

         IF(T2c.gt.0._rk) then
            Xn=self%qn/(T2n/T2c)
            Xp=self%qp/(T2p/T2c)
            nx=self%qn-(T2n/T2c)
            px=self%qp-(T2p/T2c)
         else
            Xn=0._rk
            Xp=0._rk
            nx=1._rk
            px=1._rk
         endif

          XXn=min(Xn,Xp)
         ! Photolysis
         photolysis=self%surf_phyref*(chemEIR*self%chemEIR_scaling/self%iref)*c
         ! in order to ensure mass conservation, the fraction of carbon going from T1 to T2 after photolysis is down regulated 
         ! if T1 has less  N and/or P than T2. Any excess of nutrients is redirected into the dissolved pool (N1p and N4n). (Luca, July 2018)    
         T1T2=self%phyt*min(1._rk,XXn)
         

         _SET_SURFACE_EXCHANGE_(self%id_c,-photolysis)
         _SET_SURFACE_EXCHANGE_(self%id_T2c,+photolysis*T1T2)
         _SET_SURFACE_EXCHANGE_(self%id_O3c,+photolysis*(1._rk-T1T2)/CMass)
         _SET_SURFACE_EXCHANGE_(self%id_N1p,+photolysis*(1._rk-T1T2)*self%qp+photolysis*T1T2*px)
         _SET_SURFACE_EXCHANGE_(self%id_N4n,+photolysis*(1._rk-T1T2)*self%qn+photolysis*T1T2*nx)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_surface_photolysis, photolysis)

      _HORIZONTAL_LOOP_END_    
    end subroutine

   subroutine shadow_initialize(self,configunit)
   
   ! !INPUT PARAMETERS:
      class (type_ersem_shadow_TDOC),intent(inout),target :: self
      integer,                       intent(in)           :: configunit
      
      self%dt = 3600._rk*24._rk

      call self%register_diagnostic_variable(self%id_c_shadow,'c','mgC/m^3','shadow carbon')
      call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_c_shadow)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_c_shadow,scale_factor=self%qp)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_c_shadow,scale_factor=self%qn)
      
      call self%register_model_dependency(self%id_parent,'parent')
      call self%register_state_dependency(self%id_parent_c,'parent_c','mg C/m^3','non photolabile terrigenous DOC')


      call self%request_coupling_to_model(self%id_parent_c,self%id_parent,'c')


      
      call self%register_model_dependency(self%id_T1_older_parent,'T1_older')
      call self%register_state_dependency(self%id_T1_older_parent_c,'T1_older_parent_c','mg C/m^3','non photolabile terrigenous DOC')
      call self%request_coupling_to_model(self%id_T1_older_parent_c,self%id_T1_older_parent,'c')

      call copy_fluxes(self,self%id_c_shadow,self%id_parent_c,(1._rk+self%bioaging))
      call copy_fluxes(self,self%id_c_shadow,self%id_T1_older_parent_c,self%bioaging)
      
   end subroutine shadow_initialize
   
   subroutine shadow_do(self,_ARGUMENTS_DO_)

      class (type_ersem_shadow_TDOC),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

         real (rk) :: c_shadow

! Enter horizontal loops (if any)
      _LOOP_BEGIN_
      
      _GET_(self%id_parent_c,c_shadow)
      _SET_DIAGNOSTIC_(self%id_c_shadow,c_shadow)
      
      _LOOP_END_    
    end subroutine shadow_do
      
    
end module
