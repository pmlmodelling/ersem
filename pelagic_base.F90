#include "fabm_driver.h"

module ersem_pelagic_base

   use fabm_types
   use fabm_particle

   use ersem_shared

   implicit none

!  default: all is private.
   private

   type,extends(type_particle_model),public :: type_ersem_pelagic_base
      type (type_state_variable_id)        :: id_c,id_n,id_p,id_f,id_s,id_chl
      type (type_horizontal_dependency_id) :: id_bedstress,id_wdepth
      type (type_dependency_id)            :: id_dens

      ! Target variables for sedimentation
      type (type_model_id)                      :: id_Q1,id_Q6,id_Q7
      type (type_bottom_state_variable_id)      :: id_Q6c,id_Q6n,id_Q6p,id_Q6s,id_Q6f
      type (type_bottom_state_variable_id)      :: id_Q1c,id_Q1n,id_Q1p
      type (type_bottom_state_variable_id)      :: id_Q7c,id_Q7n,id_Q7p
      
      logical  :: sedimentation = .false.
      real(rk) :: qQ7c,xR7n,xR7p
      real(rk) :: qQ1c,xR1n,xR1p
      real(rk) :: rm = 0.0_rk
   contains
      procedure :: initialize
      procedure :: do_bottom

      procedure :: initialize_ersem_base
      procedure :: get_sinking_rate
      procedure :: add_constituent
   end type
      
contains

   subroutine initialize(self,configunit)
      class (type_ersem_pelagic_base), intent(inout), target :: self
      integer,                               intent(in)            :: configunit

      character(len=10) :: composition
      real(rk)          :: c0,s0,rRPmX,EPS

      call self%get_parameter(composition,'composition','',        'elemental composition')
      call self%get_parameter(EPS,        'EPS',        'm^2/mg C','specific shortwave attenuation',default=0.0_rk)
      call self%get_parameter(rRPmX,      'rm',         'm/d',     'sinking velocity',              default=0.0_rk)

      call self%initialize_ersem_base(rm=rRPmX)

      call self%get_parameter(c0,'c0','mg C/m^3','background carbon concentration',default=0.0_rk)
      if (index(composition,'c')/=0) call self%add_constituent('c',0.0_rk,c0)
      if (index(composition,'n')/=0) call self%add_constituent('n',0.0_rk,qnRPIcX*c0)
      if (index(composition,'p')/=0) call self%add_constituent('p',0.0_rk,qpRPIcX*c0)
      if (index(composition,'s')/=0) then
         call self%get_parameter(s0,'s0','mmol Si/m^3','background silicon concentration',default=qsRPIcX*c0)
         call self%add_constituent('s',0.0_rk,s0)
      end if
      if (index(composition,'f')/=0) call self%add_constituent('f',0.0_rk)

      if (EPS/=0.0_rk) call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &
         self%id_c,scale_factor=EPS,include_background=.true.)

   end subroutine

   subroutine initialize_ersem_base(self,rm,sedimentation)
      class (type_ersem_pelagic_base), intent(inout), target :: self
      real(rk),optional,               intent(in)            :: rm
      logical, optional,               intent(in)            :: sedimentation

      ! Set time unit to d-1
      ! This implies that all rates (sink/source terms, vertical velocities) are given in d-1.
      self%dt = 86400._rk

      ! Sedimentation
      if (present(sedimentation)) then
         ! Use specified sedimentation flag
         ! If active, allow the user to override this (disable sedmentation) through run-time configuration
         if (sedimentation) call self%get_parameter(self%sedimentation,'sedimentation','','enable sedimentation',default=.true.)
      else
         ! No sedimentation flag provided; get from run-time configuration.
         call self%get_parameter(self%sedimentation,'sedimentation','','enable sedimentation',default=.false.)
      end if
      if (self%sedimentation) then
         call self%register_dependency(self%id_bedstress,standard_variables%bottom_stress)
         call self%register_dependency(self%id_dens,     standard_variables%density)
         call self%get_parameter(self%qQ1c,'qQ1c','-','fraction of sedimented matter decomposing into DOM')
         call self%get_parameter(self%qQ7c,'qQ7c','-','fraction of sedimented matter decomposing into refractory matter')
         call self%get_parameter(self%xR1n,'xR1n','-','transfer of sedimented nitrogen to DOM, relative to POM')
         call self%get_parameter(self%xR1p,'xR1p','-','transfer of sedimented phosphorus to DOM, relative to POM')
         call self%get_parameter(self%xR7n,'xR7n','-','transfer of sedimented nitrogen to refractory matter, relative to POM')
         call self%get_parameter(self%xR7p,'xR7p','-','transfer of sedimented phosphorus to refractory matter, relative to POM')

         ! Links to external benthic state variables
         call self%register_state_dependency(self%id_Q1c,'Q1c','mg C/m^2',  'Q1c')
         call self%register_state_dependency(self%id_Q1n,'Q1n','mmol N/m^2','Q1n')
         call self%register_state_dependency(self%id_Q1p,'Q1p','mmol P/m^2','Q1p')
         call self%register_state_dependency(self%id_Q6c,'Q6c','mg C/m^2',  'Q6c')
         call self%register_state_dependency(self%id_Q6n,'Q6n','mmol N/m^2','Q6n')
         call self%register_state_dependency(self%id_Q6p,'Q6p','mmol P/m^2','Q6p')
         call self%register_state_dependency(self%id_Q6s,'Q6s','mmol Si/m^2','Q6s')
         if (use_iron) call self%register_state_dependency(self%id_Q6f,'Q6f','mmol Fe/m^2','Q6f')
         call self%register_state_dependency(self%id_Q7c,'Q7c','mg C/m^2',  'Q7c')
         call self%register_state_dependency(self%id_Q7n,'Q7n','mmol N/m^2','Q7n')
         call self%register_state_dependency(self%id_Q7p,'Q7p','mmol P/m^2','Q7p')

         ! Retrieve external benthic states from individual benthic models
         call self%register_model_dependency(self%id_Q1,'Q1')
         call self%request_coupling_to_model(self%id_Q1c,self%id_Q1,'c')
         call self%request_coupling_to_model(self%id_Q1n,self%id_Q1,'n')
         call self%request_coupling_to_model(self%id_Q1p,self%id_Q1,'p')

         call self%register_model_dependency(self%id_Q6,'Q6')
         call self%request_coupling_to_model(self%id_Q6c,self%id_Q6,'c')
         call self%request_coupling_to_model(self%id_Q6n,self%id_Q6,'n')
         call self%request_coupling_to_model(self%id_Q6p,self%id_Q6,'p')
         call self%request_coupling_to_model(self%id_Q6s,self%id_Q6,'s')
         if (use_iron) call self%request_coupling_to_model(self%id_Q6f,self%id_Q6,'f')

         call self%register_model_dependency(self%id_Q7,'Q7')
         call self%request_coupling_to_model(self%id_Q7c,self%id_Q7,'c')
         call self%request_coupling_to_model(self%id_Q7n,self%id_Q7,'n')
         call self%request_coupling_to_model(self%id_Q7p,self%id_Q7,'p')
      end if

      ! Vertical velocity (positive: downwards, negative: upwards)
      if (present(rm)) self%rm = rm

   end subroutine

   subroutine add_constituent(self,name,initial_value,background_value,qn,qp)
      class (type_ersem_pelagic_base), intent(inout), target :: self
      character(len=*),                      intent(in)            :: name
      real(rk),                              intent(in)            :: initial_value
      real(rk),optional,                     intent(in)            :: background_value,qn,qp

      select case (name)
         case ('c')
            call self%register_state_variable(self%id_c,'c','mg C/m^3','carbon',initial_value,minimum=0._rk,vertical_movement=-self%rm/self%dt,background_value=background_value)
            call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_c,scale_factor=1._rk/CMass)
            if (present(qn)) call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_c,scale_factor=qn)
            if (present(qp)) call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_c,scale_factor=qp)
         case ('n')
            call self%register_state_variable(self%id_n, 'n', 'mmol N/m^3', 'nitrogen', initial_value, minimum=0._rk,vertical_movement=-self%rm/self%dt,background_value=background_value)
            call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_n)
         case ('p')
            call self%register_state_variable(self%id_p, 'p', 'mmol P/m^3', 'phosphorus', initial_value, minimum=0._rk,vertical_movement=-self%rm/self%dt,background_value=background_value)
            call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_p)
         case ('s')
            call self%register_state_variable(self%id_s,'s','mmol Si/m^3','silicate',initial_value,minimum=0._rk,vertical_movement=-self%rm/self%dt,background_value=background_value)
            call self%add_to_aggregate_variable(standard_variables%total_silicate,self%id_s)
         case ('f')
            if (use_iron) then
               call self%register_state_variable(self%id_f,'f','umol Fe/m^3','iron',initial_value,minimum=0._rk,vertical_movement=-self%rm/self%dt,background_value=background_value)
               call self%add_to_aggregate_variable(standard_variables%total_iron,self%id_f)
            end if
         case ('chl')
            call self%register_state_variable(self%id_chl,'Chl','mg/m^3','chlorophyll a',initial_value,minimum=0._rk,vertical_movement=-self%rm/self%dt,background_value=background_value)
            call self%add_to_aggregate_variable(total_chlorophyll,self%id_chl)
         case default
            call self%fatal_error('add_constituent','Unknown constituent "'//trim(name)//'".')
         end select
   end subroutine

   function get_sinking_rate(self,_ARGUMENTS_LOCAL_) result(rm)
      ! Returns sinking rate in m/d, positive for downward movement [sinking]
      class (type_ersem_pelagic_base),intent(in) :: self
      _DECLARE_ARGUMENTS_LOCAL_
      real(rk) :: rm

      rm = self%rm
   end function

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_pelagic_base), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      ! Bed characteristics - from Puls and Sundermann 1990
      ! Critical bed shear velocity = 0.01 m/s
      real(rk) :: tdep = 0.01_rk**2

      real(rk) :: tbed,density
      real(rk) :: fac,sdrate
      real(rk) :: Pc,Pn,Pp,P
      real(rk) :: fsd,fsdc,fsdn,fsdp

      if (.not.self%sedimentation) return

      _HORIZONTAL_LOOP_BEGIN_

         ! Retrieve sinking rate (at centre of cell closest to the bottom)
         fsd = self%get_sinking_rate(_ARGUMENTS_LOCAL_)

         ! Retrieve bed stress and local density - needed to determine sedimentation rate from sinking rate.
         _GET_HORIZONTAL_(self%id_bedstress,tbed)
         _GET_(self%id_dens,density)

         ! Divide actual bed stress (Pa) by density (kg/m^3) to obtain square of bed shear velocity.
         tbed = tbed/density

         if(tbed<tdep) then
            ! Bed stress is low enough to allow some sedimentation.
            fac=1._rk-tbed/tdep
         else
            ! Bed stress is too high - no actual sedimentation.
            fac=0._rk
         endif

         !sdrate = min(fsd*fac,pdepth(I)/timestep) ! Jorn: CFL criterion disabled because FABM does not provide timestep
         sdrate = fsd*fac

         Pc = 0.0_rk
         Pn = 0.0_rk
         Pp = 0.0_rk
         if (_AVAILABLE_(self%id_c)) _GET_(self%id_c,Pc)
         if (_AVAILABLE_(self%id_n)) _GET_(self%id_n,Pn)
         if (_AVAILABLE_(self%id_p)) _GET_(self%id_p,Pp)

         fsdc = sdrate*Pc
         fsdn = sdrate*Pn
         fsdp = sdrate*Pp

         _SET_BOTTOM_ODE_(self%id_Q6c, fsdc * (1._rk - self%qQ1c - self%qQ7c))
         _SET_BOTTOM_ODE_(self%id_Q6n, fsdn * (1._rk - MIN(1._rk, self%qQ1c * self%xR1n) - MIN(1._rk, self%qQ7c * self%xR7n)))
         _SET_BOTTOM_ODE_(self%id_Q6p, fsdp * (1._rk - MIN(1._rk, self%qQ1c * self%xR1p) - MIN(1._rk, self%qQ7c * self%xR7p)))
      
         _SET_BOTTOM_ODE_(self%id_Q1c, fsdc * self%qQ1c)
         _SET_BOTTOM_ODE_(self%id_Q1n, fsdn * MIN(1._rk, self%qQ1c * self%xR1n))
         _SET_BOTTOM_ODE_(self%id_Q1p, fsdp * MIN(1._rk, self%qQ1c * self%xR1p))
      
         _SET_BOTTOM_ODE_(self%id_Q7c, fsdc * self%qQ7c)
         _SET_BOTTOM_ODE_(self%id_Q7n, fsdn * MIN(1._rk, self%qQ7c * self%xR7n))
         _SET_BOTTOM_ODE_(self%id_Q7p, fsdp * MIN(1._rk, self%qQ7c * self%xR7p))
      
         _SET_BOTTOM_EXCHANGE_(self%id_c,-fsdc)
         _SET_BOTTOM_EXCHANGE_(self%id_n,-fsdn)
         _SET_BOTTOM_EXCHANGE_(self%id_p,-fsdp)

         if (_AVAILABLE_(self%id_s)) then
            ! All silicate sinks into Q6f
            _GET_(self%id_s,P)
            _SET_BOTTOM_ODE_(self%id_Q6s,sdrate*P)
            _SET_BOTTOM_EXCHANGE_(self%id_s,-sdrate*P)
         end if

         if (_AVAILABLE_(self%id_f)) then
            ! All iron sinks into Q6f
            _GET_(self%id_f,P)
            _SET_BOTTOM_ODE_(self%id_Q6f,sdrate*P)
            _SET_BOTTOM_EXCHANGE_(self%id_f,-sdrate*P)
         end if

         if (_AVAILABLE_(self%id_chl)) then
            ! Chlorophyll is simply lost by sinking:
            _GET_(self%id_chl,P)
            _SET_BOTTOM_EXCHANGE_(self%id_chl,-sdrate*P)
         end if

      _HORIZONTAL_LOOP_END_
   end subroutine

end module