#include "fabm_driver.h"
!#define IRON
module pml_ersem_base

   use fabm_types
   use fabm_expressions

   implicit none

!  default: all is private.
   private

   real(rk),parameter,public :: CMass = 12.011_rk
   real(rk),parameter,public :: qnRPIcX = 1.26E-02_rk
   real(rk),parameter,public :: qpRPIcX = 7.86E-04_rk
   real(rk),parameter,public :: qsRPIcX=15._rk/106._rk/CMass

   type,extends(type_base_model),public :: type_ersem_pelagic_base_model
      type (type_state_variable_id)      :: id_c,id_n,id_p,id_f,id_s,id_chl
      type (type_diagnostic_variable_id) :: id_cD,id_nD,id_pD,id_fD,id_sD,id_chlD,id_lD
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
      type (type_bottom_state_variable_id) :: id_c,id_n,id_p,id_f,id_s,id_l
      type (type_conserved_quantity_id)    :: id_totc,id_totn,id_totp,id_tots,id_totf
   contains
      procedure :: initialize_ersem_benthic_base
      procedure :: add_constituent => benthic_base_add_constituent
   end type
      
contains
      
   subroutine initialize_ersem_base(self,c_ini,n_ini,p_ini,s_ini,f_ini,chl_ini,qn,qp,w,sedimentation,c0,n0,p0,s0,f0,chl0)
      class (type_ersem_pelagic_base_model), intent(inout), target :: self
      real(rk),optional,             intent(in)            :: c_ini,n_ini,p_ini,s_ini,f_ini,chl_ini,w
      real(rk),optional,             intent(in)            :: c0,n0,p0,s0,f0,chl0
      real(rk),optional,             intent(in)            :: qn,qp
      logical,optional,              intent(in)            :: sedimentation
      
      type (type_weighted_sum),pointer :: child

      self%dt = 86400._rk

      ! Sedimentation
      self%sedimentation = .false.
      if (present(sedimentation)) self%sedimentation = sedimentation
      if (self%sedimentation) call self%get_parameter(self%sedimentation,'sedimentation',default=.true.)
      if (self%sedimentation) then
         call self%register_dependency(self%id_bedstress,standard_variables%bottom_stress)
         call self%register_dependency(self%id_dens,     standard_variables%density)
         call self%get_parameter(self%xR1nX,'xR1n')
         call self%get_parameter(self%xR1pX,'xR1p')
         call self%get_parameter(self%xR7nX,'xR7n')
         call self%get_parameter(self%xR7pX,'xR7p')
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
#ifdef IRON   
         call self%register_state_dependency(self%id_Q6f,'Q6f','mmol Fe m-2','Q6f')
#endif
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
#ifdef IRON   
         call self%request_coupling(self%id_Q6f,'f',source=self%id_Q6)
#endif

         call self%register_model_dependency(self%id_Q7,'Q7')
         call self%request_coupling(self%id_Q7c,'c',source=self%id_Q7)
         call self%request_coupling(self%id_Q7n,'n',source=self%id_Q7)
         call self%request_coupling(self%id_Q7p,'p',source=self%id_Q7)
      end if

      ! Vertical velocity (positive: upwards, negative: downwards)
      self%w = 0.0_rk
      if (present(w)) self%w = w

      ! Carbon
      if (present(c_ini)) then
         call self%register_state_variable(self%id_c,'c','mg m-3','carbon',c_ini,minimum=0._rk,vertical_movement=self%w/self%dt,background_value=c0)
         call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_c,scale_factor=1._rk/CMass)
      else
         call self%register_diagnostic_variable(self%id_cD,'c','mg m-3','carbon',missing_value=0._rk,output=output_none)
      end if

      ! Nitrogen
      if (present(n_ini)) then
         ! Nitrogen as state variable
         call self%register_state_variable(self%id_n, 'n', 'mmol N m-3', 'nitrogen', n_ini, minimum=0._rk,vertical_movement=self%w/self%dt,background_value=n0)
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_n)
      elseif (present(qn)) then
         ! Nitrogen derived from carbon (constant C:P)
         ! Create an expression that computes N content, so external models (e.g., predators) can access it.
         allocate(child)
         child%output_name = 'n'
         child%output_units = 'mmol N m-3'
         call child%add_component('c',qn)
         call self%add_child(child,'dynamic_n',configunit=-1)
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_c,scale_factor=qn)
      else
         ! No nitrogen
         call self%register_diagnostic_variable(self%id_nD, 'n', 'mmol N m-3', 'nitrogen', missing_value=0._rk,output=output_none)
      end if

      ! Phosphorus
      if (present(p_ini)) then
         ! Phosphorous as state variable
         call self%register_state_variable(self%id_p, 'p', 'mmol P m-3', 'phosphorus', p_ini, minimum=0._rk,vertical_movement=self%w/self%dt,background_value=p0)
         call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_p)
      elseif (present(qp)) then
         ! Phosphorous derived from carbon (constant C:P)
         ! Create an expression that computes P content, so external models (e.g., predators) can access it.
         allocate(child)
         child%output_name = 'p'
         child%output_units = 'mmol P m-3'
         call child%add_component('c',qp)
         call self%add_child(child,'dynamic_p',configunit=-1)
         call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_c,scale_factor=qp)
      else
         ! No phosphorous
         call self%register_diagnostic_variable(self%id_pD, 'p', 'mmol P m-3', 'phosphorus', missing_value=0._rk,output=output_none)
      end if

      ! Silicate
      if (present(s_ini)) then
         call self%register_state_variable(self%id_s,'s','mmol Si m-3','silicate',s_ini,minimum=0._rk,vertical_movement=self%w/self%dt,background_value=s0)
         call self%add_to_aggregate_variable(standard_variables%total_silicate,self%id_s)
      else
         call self%register_diagnostic_variable(self%id_sD,'s','mmol Si m-3','silicate',missing_value=0._rk,output=output_none)
      end if

      ! Chlorophyll
      if (present(chl_ini)) then
         call self%register_state_variable(self%id_chl,'Chl','mg C m-3','chlorophyll-a',chl_ini,minimum=0._rk,vertical_movement=self%w/self%dt,background_value=chl0)
      else
         call self%register_diagnostic_variable(self%id_chlD,'Chl','mg C m-3','chlorophyll-a',missing_value=0._rk,output=output_none)
      end if

#ifdef IRON   
      ! Iron
      if (present(f_ini)) then
         ! Iron as state variable
         call self%register_state_variable(self%id_f,'f','umol Fe m-3','iron',f_ini,minimum=0._rk,vertical_movement=self%w/self%dt,background_value=f0)
         call self%add_to_aggregate_variable(standard_variables%total_iron,self%id_f)
      else
         ! No iron
         call self%register_diagnostic_variable(self%id_fD,'f','umol Fe m-3','iron',missing_value=0._rk,output=output_none)
      end if
#endif

      call self%register_diagnostic_variable(self%id_lD,'l','mg C m-3','calcite',missing_value=0._rk,output=output_none)

   end subroutine

   function get_sinking_rate(self,_ARGUMENTS_LOCAL_) result(w)
      ! Returns sinking rate in m/d, positive for downward movement [sinking]
      class (type_ersem_pelagic_base_model),intent(in) :: self
      _DECLARE_ARGUMENTS_LOCAL_
      real(rk) :: w

      w = -self%w
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
         if (_AVAILABLE_(self%id_c)) _GET_(self%id_c,Pc)
         if (_AVAILABLE_(self%id_n)) _GET_(self%id_n,Pn)
         if (_AVAILABLE_(self%id_p)) _GET_(self%id_p,Pp)
         if (_AVAILABLE_(self%id_s)) _GET_(self%id_s,Ps)

         fsd = self%get_sinking_rate(_ARGUMENTS_LOCAL_)

         ! Divide actual stress (Pa) by density (kg m-3) to obtain square of bed shear velocity.
         tbed = tbed/density
!
!     Bed characteristics - from Puls and Sundermann 1990
!     Critical bed shear velocity = 0.01 m/s
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
      if (present(c_ini)) call self%add_constituent('c',c_ini)
      if (present(n_ini)) call self%add_constituent('n',n_ini)
      if (present(p_ini)) call self%add_constituent('p',p_ini)
      if (present(s_ini)) call self%add_constituent('s',s_ini)
      if (present(f_ini)) call self%add_constituent('f',f_ini)
   end subroutine

   subroutine benthic_base_add_constituent(self,name,initial_value)
      class (type_ersem_benthic_base_model), intent(inout), target :: self
      character(len=*),                      intent(in)            :: name
      real(rk),                              intent(in)            :: initial_value

      select case (name)
         case ('c')
            call self%register_state_variable(self%id_c,'c','mg C m-2','carbon',initial_value,minimum=0._rk)
            call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_c,scale_factor=1._rk/CMass)
         case ('n')
            call self%register_state_variable(self%id_n, 'n', 'mmol N m-2', 'nitrogen', initial_value, minimum=0._rk)
            call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_n)
         case ('p')
            call self%register_state_variable(self%id_p, 'p', 'mmol P m-2', 'phosphorus', initial_value, minimum=0._rk)
            call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_p)
         case ('s')
            call self%register_state_variable(self%id_s,'s','mmol Si m-2','silicate',initial_value,minimum=0._rk)
            call self%add_to_aggregate_variable(standard_variables%total_silicate,self%id_s)
         case ('f')
#ifdef IRON   
            call self%register_state_variable(self%id_f,'f','umol Fe m-2','iron',initial_value,minimum=0._rk)
            call self%add_to_aggregate_variable(standard_variables%total_iron,self%id_f)
#endif
         case ('l')
            call self%register_state_variable(self%id_l,'l','mg C m-2','calcite',initial_value,minimum=0._rk)
            call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_l,scale_factor=1._rk/CMass)
         case default
            call self%fatal_error('benthic_base_add_constituent','Unknown element "'//trim(name)//'".')
      end select
   end subroutine benthic_base_add_constituent
   
end module