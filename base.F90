#include "fabm_driver.h"
#define IRON
module pml_ersem_base

   use fabm_types

   implicit none

!  default: all is private.
   private

   type,extends(type_base_model),public :: type_ersem_pelagic_base_model
      type (type_state_variable_id)      :: id_c,id_n,id_p,id_f,id_s,id_chl
      type (type_diagnostic_variable_id) :: id_cD,id_nD,id_pD,id_fD,id_sD,id_chlD
      type (type_conserved_quantity_id)  :: id_totc,id_totn,id_totp,id_tots,id_totf
      type (type_horizontal_dependency_id) :: id_bedstress,id_wdepth
      type (type_dependency_id)            :: id_dens

      ! Target variables for sedimentation
      type (type_model_id)                      :: id_Q1,id_Q6,id_Q7
      type (type_bottom_state_variable_id)      :: id_Q6c,id_Q6n,id_Q6p,id_Q6s,id_Q6f
      type (type_bottom_state_variable_id)      :: id_Q1c,id_Q1n,id_Q1p
      type (type_bottom_state_variable_id)      :: id_Q7c,id_Q7n,id_Q7p
      
      logical  :: sedimentation
      real(rk) :: qQ7c,xR7nX,xR7pX
      real(rk) :: qQ1c,xR1nX,xR1pX
      real(rk) :: w
   contains
      procedure :: initialize_ersem_base
      procedure :: do_bottom
      procedure :: get_sinking_rate
   end type

   type,extends(type_base_model),public :: type_ersem_benthic_base_model
      type (type_bottom_state_variable_id) :: id_c,id_n,id_p,id_f,id_s,id_chl
      type (type_conserved_quantity_id)    :: id_totc,id_totn,id_totp,id_tots,id_totf
   contains
      procedure :: initialize_ersem_benthic_base
   end type

   type type_component
      character(len=attribute_length) :: name   = ''
      real(rk)                        :: weight = 1._rk
      type (type_dependency_id)       :: id
      type (type_component),pointer   :: next   => null()
   end type

   type,extends(type_base_model) :: type_weighted_sum
      character(len=attribute_length) :: output_name      = ''
      character(len=attribute_length) :: output_long_name = ''
      character(len=attribute_length) :: output_units     = ''
      type (type_diagnostic_variable_id) :: id_output
      type (type_component),pointer   :: first => null()
   contains
      procedure :: initialize    => weighted_sum_initialize
      procedure :: do            => weighted_sum_do
      procedure :: add_component => weighted_sum_add_component
   end type
      
contains
      
   subroutine initialize_ersem_base(self,c_ini,n_ini,p_ini,s_ini,f_ini,chl_ini,qn,qp,w,sedimentation)
      class (type_ersem_pelagic_base_model), intent(inout), target :: self
      real(rk),optional,             intent(in)            :: c_ini,n_ini,p_ini,s_ini,f_ini,chl_ini,w
      real(rk),optional,             intent(in)            :: qn,qp
      logical,optional,              intent(in)            :: sedimentation
      
      type (type_weighted_sum),pointer :: child

      self%dt = 86400._rk

      ! Sedimentation
      self%sedimentation = .false.
      if (present(sedimentation)) self%sedimentation = sedimentation
      if (self%sedimentation) then
         call self%register_dependency(self%id_bedstress,standard_variables%bottom_stress)
         call self%register_dependency(self%id_dens,     standard_variables%density)
         call self%get_parameter(self%xR1nX,'xR1nX')
         call self%get_parameter(self%xR1pX,'xR1pX')
         call self%get_parameter(self%xR7nX,'xR7nX')
         call self%get_parameter(self%xR7pX,'xR7pX')
         call self%get_parameter(self%qQ1c, 'qQ1c')
         call self%get_parameter(self%qQ7c, 'qQ7c')

         ! Links to external benthic state variables
         call self%register_state_dependency(self%id_Q1c,'Q1c','mg C m-2',  'Q1c')
         call self%register_state_dependency(self%id_Q1n,'Q1n','mmol N m-2','Q1n')
         call self%register_state_dependency(self%id_Q1p,'Q1p','mmol P m-2','Q1p')
         call self%register_state_dependency(self%id_Q6c,'Q6c','mg C m-2',  'Q6c')
         call self%register_state_dependency(self%id_Q6n,'Q6n','mmol N m-2','Q6n')
         call self%register_state_dependency(self%id_Q6p,'Q6p','mmol P m-2','Q6p')
         call self%register_state_dependency(self%id_Q6s,'Q6s','mmol Si m-2','Q6s')
         call self%register_state_dependency(self%id_Q6f,'Q6f','mmol Fe m-2','Q6f')
         call self%register_state_dependency(self%id_Q7c,'Q7c','mg C m-2',  'Q7c')
         call self%register_state_dependency(self%id_Q7n,'Q7n','mmol N m-2','Q7n')
         call self%register_state_dependency(self%id_Q7p,'Q7p','mmol P m-2','Q7p')

         ! Retrieve external benthic states from individual benthic models
         call self%register_model_dependency(self%id_Q1,'Q1')
         call self%request_coupling(self%id_Q1c,'c',source=self%id_Q1)
         call self%request_coupling(self%id_Q1n,'n',source=self%id_Q1)
         call self%request_coupling(self%id_Q1p,'p',source=self%id_Q1)

         call self%register_model_dependency(self%id_Q6,'Q6')
         call self%request_coupling(self%id_Q6c,'c',source=self%id_Q6)
         call self%request_coupling(self%id_Q6n,'n',source=self%id_Q6)
         call self%request_coupling(self%id_Q6p,'p',source=self%id_Q6)
         call self%request_coupling(self%id_Q6s,'s',source=self%id_Q6)
         call self%request_coupling(self%id_Q6f,'f',source=self%id_Q6)

         call self%register_model_dependency(self%id_Q7,'Q7')
         call self%request_coupling(self%id_Q7c,'c',source=self%id_Q7)
         call self%request_coupling(self%id_Q7n,'n',source=self%id_Q7)
         call self%request_coupling(self%id_Q7p,'p',source=self%id_Q7)
      end if

      ! Conserved quantities
      call self%register_conserved_quantity(self%id_totc,standard_variables%total_carbon)
      call self%register_conserved_quantity(self%id_totn,standard_variables%total_nitrogen)
      call self%register_conserved_quantity(self%id_totp,standard_variables%total_phosphorus)
      call self%register_conserved_quantity(self%id_tots,standard_variables%total_silicate)
      call self%register_conserved_quantity(self%id_totf,standard_variables%total_iron)

      ! Vertical velocity (positive: upwards, negative: downwards)
      self%w = 0.0_rk
      if (present(w)) self%w = w

      ! Carbon
      if (present(c_ini)) then
         call self%register_state_variable(self%id_c,'c','mg m-3','carbon',c_ini,minimum=0._rk,vertical_movement=self%w/self%dt)
         call self%add_conserved_quantity_component(self%id_totc,self%id_c,scale_factor=1._rk/12._rk)
      else
         call self%register_diagnostic_variable(self%id_cD,'c','mg m-3','carbon',missing_value=0._rk)
      end if

      ! Nitrogen
      if (present(n_ini)) then
         ! Nitrogen as state variable
         call self%register_state_variable(self%id_n, 'n', 'mmol N m-3', 'nitrogen', n_ini, minimum=0._rk,vertical_movement=self%w/self%dt)
         call self%add_conserved_quantity_component(self%id_totn,self%id_n)
      elseif (present(qn)) then
         ! Nitrogen derived from carbon (constant C:P)
         ! Create an expression that computes N content, so external models (e.g., predators) can access it.
         allocate(child)
         child%output_name = 'n'
         child%output_units = 'mmol N m-3'
         call child%add_component('c',qn)
         call self%add_child(child,'dynamic_n',configunit=-1)
         call self%add_conserved_quantity_component(self%id_totn,self%id_c,scale_factor=qn)
      else
         ! No nitrogen
         call self%register_diagnostic_variable(self%id_nD, 'n', 'mmol N m-3', 'nitrogen', missing_value=0._rk)
      end if

      ! Phosphorus
      if (present(p_ini)) then
         ! Phosphorous as state variable
         call self%register_state_variable(self%id_p, 'p', 'mmol P m-3', 'phosphorus', p_ini, minimum=0._rk,vertical_movement=self%w/self%dt)
         call self%add_conserved_quantity_component(self%id_totp,self%id_p)
      elseif (present(qp)) then
         ! Phosphorous derived from carbon (constant C:P)
         ! Create an expression that computes P content, so external models (e.g., predators) can access it.
         allocate(child)
         child%output_name = 'p'
         child%output_units = 'mmol P m-3'
         call child%add_component('c',qp)
         call self%add_child(child,'dynamic_p',configunit=-1)
         call self%add_conserved_quantity_component(self%id_totp,self%id_c,scale_factor=qp)
      else
         ! No phosphorous
         call self%register_diagnostic_variable(self%id_pD, 'p', 'mmol P m-3', 'phosphorus', missing_value=0._rk)
      end if

      ! Silicate
      if (present(s_ini)) then
         call self%register_state_variable(self%id_s,'s','mmol Si m-3','silicate',s_ini,minimum=0._rk,vertical_movement=self%w/self%dt)
         call self%add_conserved_quantity_component(self%id_tots,self%id_s)
      else
         call self%register_diagnostic_variable(self%id_sD,'s','mmol Si m-3','silicate',missing_value=0._rk)
      end if

      ! Chlorophyll
      if (present(chl_ini)) then
         call self%register_state_variable(self%id_chl,'Chl','mg C m-3','chlorophyll-a',chl_ini,minimum=0._rk,vertical_movement=self%w/self%dt)
      else
         call self%register_diagnostic_variable(self%id_chlD,'Chl','mg C m-3','chlorophyll-a',missing_value=0._rk)
      end if

#ifdef IRON   
      ! Iron
      if (present(f_ini)) then
         ! Iron as state variable
         call self%register_state_variable(self%id_f,'f','umol Fe m-3','iron',f_ini,minimum=0._rk,vertical_movement=self%w/self%dt)
         call self%add_conserved_quantity_component(self%id_totf,self%id_f)
      else
         ! No iron
         call self%register_diagnostic_variable(self%id_fD,'f','umol Fe m-3','iron',missing_value=0._rk)
      end if
#endif

   end subroutine
   
   function get_sinking_rate(self,_ARGUMENTS_LOCAL_) result(w)
      class (type_ersem_pelagic_base_model),intent(in) :: self
      _DECLARE_ARGUMENTS_LOCAL_
      real(rk) :: w

      w = self%w
   end function

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_pelagic_base_model), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: tbed,density
      real(rk) :: tdep,fac,sdrate
      real(rk) :: Pc,Pn,Pp,Ps
      real(rk) :: fsd,fsdc,fsdn,fsdp,fsds

      if (.not.self%sedimentation) return

      _HORIZONTAL_LOOP_BEGIN_

         _GET_HORIZONTAL_(self%id_bedstress,tbed)
         _GET_(self%id_dens,density)
         
         Pc = 0.0_rk
         Pn = 0.0_rk
         Pp = 0.0_rk
         Ps = 0.0_rk
         if (_AVAILABLE_(self%id_c)) _GET_(self%id_p,Pc)
         if (_AVAILABLE_(self%id_n)) _GET_(self%id_n,Pn)
         if (_AVAILABLE_(self%id_p)) _GET_(self%id_p,Pp)
         if (_AVAILABLE_(self%id_s)) _GET_(self%id_s,Ps)

         fsd = -self%get_sinking_rate(_ARGUMENTS_LOCAL_)

         ! From actual stress (Pa) to shear velocity (m/s)
         tbed = sqrt(tbed/density)
         tbed = 0.0_rk
!
!     Bed characteristics - from Puls and Sundermann 1990
!     Critical stress for deposition (m/s)
!
      tdep=0.01_rk**2
!
     if(tbed .lt. tdep) then
       fac=1._rk-tbed/tdep
     else
       fac=0._rk
     endif
!
      !sdrate = min(fsd*fac,pdepth(I)/timestep) ! Jorn: CFL criterion disabled because FABM does not provide timestep
      sdrate = fsd*fac
!
      fsdc = sdrate*Pc
      fsdn = sdrate*Pn
      fsdp = sdrate*Pp
      fsds = sdrate*Ps
      
      _SET_ODE_BEN_(self%id_Q6c, fsdc * (1._rk - self%qQ1c - self%qQ7c))
      _SET_ODE_BEN_(self%id_Q6n, fsdn * (1._rk - MIN(1._rk, self%qQ1c * self%xR1nX) - MIN(1._rk, self%qQ7c * self%xR7nX)))
      _SET_ODE_BEN_(self%id_Q6p, fsdp * (1._rk - MIN(1._rk, self%qQ1c * self%xR1pX) - MIN(1._rk, self%qQ7c * self%xR7pX)))
      _SET_ODE_BEN_(self%id_Q6s, fsds)
      
      _SET_ODE_BEN_(self%id_Q1c, fsdc * self%qQ1c)
      _SET_ODE_BEN_(self%id_Q1n, fsdn * MIN(1._rk, self%qQ1c * self%xR1nX))
      _SET_ODE_BEN_(self%id_Q1p, fsdp * MIN(1._rk, self%qQ1c * self%xR1pX))
      
      _SET_ODE_BEN_(self%id_Q7c, fsdc * self%qQ7c)
      _SET_ODE_BEN_(self%id_Q7n, fsdn * MIN(1._rk, self%qQ7c * self%xR7nX))
      _SET_ODE_BEN_(self%id_Q7p, fsdp * MIN(1._rk, self%qQ7c * self%xR7pX))
      
      _SET_BOTTOM_EXCHANGE_(self%id_c,-fsdc)
      _SET_BOTTOM_EXCHANGE_(self%id_n,-fsdn)
      _SET_BOTTOM_EXCHANGE_(self%id_p,-fsdp)
      if (_AVAILABLE_(self%id_s)) _SET_BOTTOM_EXCHANGE_(self%id_s,-fsds)
      
      _HORIZONTAL_LOOP_END_
   end subroutine
   
   subroutine initialize_ersem_benthic_base(self,c_ini,n_ini,p_ini,s_ini,f_ini)
      class (type_ersem_benthic_base_model), intent(inout), target :: self
      real(rk),optional,                     intent(in)            :: c_ini,n_ini,p_ini,s_ini,f_ini

      self%dt = 86400._rk

      call self%register_conserved_quantity(self%id_totc,standard_variables%total_carbon)
      call self%register_conserved_quantity(self%id_totn,standard_variables%total_nitrogen)
      call self%register_conserved_quantity(self%id_totp,standard_variables%total_phosphorus)
      call self%register_conserved_quantity(self%id_tots,standard_variables%total_silicate)
      call self%register_conserved_quantity(self%id_totf,standard_variables%total_iron)

      if (present(c_ini)) then
         call self%register_state_variable(self%id_c,'c','mg m-2','carbon',c_ini,minimum=0._rk)
         call self%add_conserved_quantity_component(self%id_totc,self%id_c,scale_factor=1._rk/12._rk)
      end if

      if (present(s_ini)) then
         call self%register_state_variable(self%id_s,'s','mmol Si m-2','silicate',s_ini,minimum=0._rk)
         call self%add_conserved_quantity_component(self%id_tots,self%id_s)
      end if

#ifdef IRON   
      if (present(f_ini)) then
         ! Iron as state variable
         call self%register_state_variable(self%id_f,'f','umol Fe m-2','iron',f_ini,minimum=0._rk)
         call self%add_conserved_quantity_component(self%id_totf,self%id_f)
      end if
#endif

      if (present(p_ini)) then
         ! Phosphorous as state variable
         call self%register_state_variable(self%id_p, 'p', 'mmol P m-2', 'phosphorus', p_ini, minimum=0._rk)
         call self%add_conserved_quantity_component(self%id_totp,self%id_p)
      end if

      if (present(n_ini)) then
         ! Nitrogen as state variable
         call self%register_state_variable(self%id_n, 'n', 'mmol N m-2', 'nitrogen', n_ini, minimum=0._rk)
         call self%add_conserved_quantity_component(self%id_totn,self%id_n)
      end if

   end subroutine

   subroutine weighted_sum_initialize(self,configunit)
      class (type_weighted_sum),intent(inout),target :: self
      integer,                  intent(in)           :: configunit

      type (type_component),pointer :: component
      component => self%first
      do while (associated(component))
         call self%register_dependency(component%id,trim(component%name))
         component => component%next
      end do
      if (self%output_long_name=='') self%output_long_name = self%output_name
      call self%register_diagnostic_variable(self%id_output,self%output_name,self%output_units,self%output_long_name)
      call self%parent%add_alias(self%id_output,trim(self%output_name))
   end subroutine

   subroutine weighted_sum_add_component(self,name,weight)
      class (type_weighted_sum),intent(inout) :: self
      character(len=*),         intent(in)    :: name
      real(rk),optional,        intent(in)    :: weight
      
      type (type_component),pointer :: component

      if (.not.associated(self%first)) then
         allocate(self%first)
         component => self%first
      else
         component => self%first
         do while (associated(component%next))
            component => component%next
         end do
         allocate(component%next)
         component => component%next
      end if
      component%name = name
      if (present(weight)) component%weight = weight
   end subroutine

   subroutine weighted_sum_do(self,_ARGUMENTS_DO_)
      class (type_weighted_sum),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      
      type (type_component),pointer :: component
      real(rk)                      :: sum,value
      
      _LOOP_BEGIN_
         sum = 0._rk
         component => self%first
         do while (associated(component))
            _GET_(component%id,value)
            sum = sum + component%weight*value
            component => component%next
         end do
         _SET_DIAGNOSTIC_(self%id_output,sum)
      _LOOP_END_
   end subroutine
end module