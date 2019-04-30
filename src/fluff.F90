#include "fabm_driver.h"

module ersem_fluff

   use fabm_types
   use fabm_particle
   use ersem_shared
   use ersem_benthic_base

   implicit none

   private

   type,extends(type_ersem_benthic_base),public :: type_ersem_fluff

      ! Target variables for sedimentation
      type (type_model_id),allocatable,dimension(:)        :: id_target
      type (type_bottom_state_variable_id),allocatable,dimension(:)      :: id_targetc,id_targetn,id_targetp,id_targets,id_targetf
      type (type_horizontal_diagnostic_variable_id),allocatable,dimension(:) :: id_inc_c(:),id_inc_n(:),id_inc_p(:),id_inc_s(:),id_inc_f(:)
      real(rk) :: sd
      integer :: ntarget
      real(rk),allocatable :: qxc(:),qxn(:),qxp(:),qxs(:),qxf(:)
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self,configunit)
      class (type_ersem_fluff), intent(inout), target :: self
      integer,                               intent(in)            :: configunit
      integer :: itarget
      character(len=16) :: num
         call self%type_ersem_benthic_base%initialize(configunit)

         call self%get_parameter(self%sd,'sd','d-1','rate of fluff incorporation into sediments')

         call self%get_parameter(self%ntarget,'ntarget','', 'number of incorporation targets',default=0)
         allocate(self%qxc(self%ntarget))
         allocate(self%qxn(self%ntarget))
         allocate(self%qxp(self%ntarget))
         allocate(self%qxs(self%ntarget))
       if (use_iron) allocate(self%qxf(self%ntarget))

       do itarget = 2,self%ntarget
         write(num,'(i0)') itarget
          call self%get_parameter(self%qxc(itarget),'qxc'//trim(num),'-','fraction of carbon incorporation into target '//trim(num))
          call self%get_parameter(self%qxn(itarget),'qxn'//trim(num),'-','fraction of nitrogen incorporation into target '//trim(num))
          call self%get_parameter(self%qxp(itarget),'qxp'//trim(num),'-','fraction of phosphorus incorporation into target '//trim(num))
          call self%get_parameter(self%qxs(itarget),'qxs'//trim(num),'-','fraction of silicate incorporation into target '//trim(num))
          if (use_iron) call self%get_parameter(self%qxf(itarget),'qxf'//trim(num),'-','fraction of iron incorporation into target '//trim(num))
      end do

          self%qxc(1)=1._rk-sum(self%qxc(2:))
          self%qxn(1)=1._rk-sum(self%qxn(2:))
          self%qxp(1)=1._rk-sum(self%qxp(2:))
          self%qxs(1)=1._rk-sum(self%qxs(2:))
       if (use_iron) self%qxf(1)=1._rk-sum(self%qxf(2:))

         allocate(self%id_target(self%ntarget))
         allocate(self%id_targetc(self%ntarget))
         allocate(self%id_targetn(self%ntarget))
         allocate(self%id_targetp(self%ntarget))
         allocate(self%id_targets(self%ntarget))
      if (use_iron) allocate(self%id_targetf(self%ntarget))

         allocate(self%id_inc_c(self%ntarget))
         allocate(self%id_inc_n(self%ntarget))
         allocate(self%id_inc_p(self%ntarget))
         allocate(self%id_inc_s(self%ntarget))
      if (use_iron) allocate(self%id_inc_f(self%ntarget))

      do itarget=1,self%ntarget
        write(num,'(i0)') itarget
        call self%register_model_dependency(self%id_target(itarget),'target'//trim(num))
        if (self%qxc(itarget)/=0) then
         call self%register_state_dependency(self%id_targetc(itarget),'target'//trim(num)//'c','mg C/m^2','target '//trim(num)//' carbon')
         call self%request_coupling_to_model(self%id_targetc(itarget),self%id_target(itarget),'c')
         call self%register_diagnostic_variable(self%id_inc_c(itarget),'into_target'//trim(num)//'c','mg C/m^2/d','rate of fluff incorporation into target '//trim(num)//' carbon',domain=domain_bottom,source=source_do_bottom)
        end if
        if (self%qxn(itarget)/=0) then
         call self%register_state_dependency(self%id_targetn(itarget),'target'//trim(num)//'n','mg N/m^2','target '//trim(num)//' nitrogen')
         call self%request_coupling_to_model(self%id_targetn(itarget),self%id_target(itarget),'n')
         call self%register_diagnostic_variable(self%id_inc_n(itarget),'into_target'//trim(num)//'n','mg N/m^2/d','rate of fluff incorporation into target '//trim(num)//' nitrogen',domain=domain_bottom,source=source_do_bottom)
        end if
        if (self%qxp(itarget)/=0) then
         call self%register_state_dependency(self%id_targetp(itarget),'target'//trim(num)//'p','mg P/m^2','target '//trim(num)//' phosphorus')
         call self%request_coupling_to_model(self%id_targetp(itarget),self%id_target(itarget),'p')
         call self%register_diagnostic_variable(self%id_inc_p(itarget),'into_target'//trim(num)//'p','mg P/m^2/d','rate of fluff incorporation into target '//trim(num)//' phosphorus',domain=domain_bottom,source=source_do_bottom)
        end if
        if (self%qxs(itarget)/=0) then
         call self%register_state_dependency(self%id_targets(itarget),'target'//trim(num)//'s','mg Si/m^2','target '//trim(num)//' silicate')
         call self%request_coupling_to_model(self%id_targets(itarget),self%id_target(itarget),'s')
         call self%register_diagnostic_variable(self%id_inc_s(itarget),'into_target'//trim(num)//'s','mg Si/m^2/d','rate of fluff incorporation into target '//trim(num)//' silicate',domain=domain_bottom,source=source_do_bottom)
        end if
        if (use_iron) then
        if (self%qxf(itarget)/=0) then
         call self%register_state_dependency(self%id_targetf(itarget),'target'//trim(num)//'f','mg Fe/m^2','target '//trim(num)//' iron')
         call self%request_coupling_to_model(self%id_targetf(itarget),self%id_target(itarget),'f')
         call self%register_diagnostic_variable(self%id_inc_f(itarget),'into_target'//trim(num)//'f','mg Fe/m^2/d','rate of fluff incorporation into target '//trim(num)//' iron',domain=domain_bottom,source=source_do_bottom)
        end if
        end if

      end do

   end subroutine initialize

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)

      class (type_ersem_fluff), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: fsdc,fsdp,fsdn,fsds,fsdf
      real(rk) :: Pc,Pn,Pp,Ps,Pf
      integer :: itarget
      call self%type_ersem_benthic_base%do_bottom(_ARGUMENTS_DO_BOTTOM_)

      _HORIZONTAL_LOOP_BEGIN_

         Pc = 0.0_rk
         Pn = 0.0_rk
         Pp = 0.0_rk
         Ps = 0.0_rk
         Pf = 0.0_rk

         if (_VARIABLE_REGISTERED_(self%id_c)) then
            _GET_HORIZONTAL_(self%id_c,Pc)
            fsdc = self%sd*Pc
            do itarget=1,self%ntarget
              if (self%qxc(itarget)/=0) then
              _SET_BOTTOM_ODE_(self%id_targetc(itarget), fsdc*self%qxc(itarget))
              _SET_HORIZONTAL_DIAGNOSTIC_(self%id_inc_c(itarget),fsdc*self%qxc(itarget))
              end if
            end do
            _SET_BOTTOM_ODE_(self%id_c,-fsdc)
         end if

         if (_VARIABLE_REGISTERED_(self%id_n)) then
            _GET_HORIZONTAL_(self%id_n,Pn)
            fsdn = self%sd*Pn
            do itarget=1,self%ntarget
              if (self%qxn(itarget)/=0) then
              _SET_BOTTOM_ODE_(self%id_targetn(itarget), fsdn*self%qxn(itarget))
              _SET_HORIZONTAL_DIAGNOSTIC_(self%id_inc_n(itarget),fsdn*self%qxn(itarget))
              end if
            end do
            _SET_BOTTOM_ODE_(self%id_n,-fsdn)
         end if

         if (_VARIABLE_REGISTERED_(self%id_p)) then
            _GET_HORIZONTAL_(self%id_p,Pp)
            fsdp = self%sd*Pp
            do itarget=1,self%ntarget
              if (self%qxp(itarget)/=0) then
              _SET_BOTTOM_ODE_(self%id_targetp(itarget), fsdp*self%qxp(itarget))
              _SET_HORIZONTAL_DIAGNOSTIC_(self%id_inc_p(itarget),fsdp*self%qxp(itarget))
              end if
            end do
            _SET_BOTTOM_ODE_(self%id_p,-fsdp)
         end if

         if (_VARIABLE_REGISTERED_(self%id_s)) then
            _GET_HORIZONTAL_(self%id_s,Ps)
            fsds = self%sd*Ps
            do itarget=1,self%ntarget
              if (self%qxs(itarget)/=0) then
              _SET_BOTTOM_ODE_(self%id_targets(itarget), fsds*self%qxs(itarget))
              _SET_HORIZONTAL_DIAGNOSTIC_(self%id_inc_s(itarget),fsds*self%qxs(itarget))
              end if
            end do
            _SET_BOTTOM_ODE_(self%id_s,-fsds)
         end if

         if (use_iron) then
         if (_VARIABLE_REGISTERED_(self%id_f)) then
            _GET_HORIZONTAL_(self%id_f,Pf)
            fsdf = self%sd*Pf
            do itarget=1,self%ntarget
              if (self%qxf(itarget)/=0) then
              _SET_BOTTOM_ODE_(self%id_targetf(itarget), fsdf*self%qxf(itarget))
              _SET_HORIZONTAL_DIAGNOSTIC_(self%id_inc_f(itarget),fsdf*self%qxf(itarget))
              end if
            end do
            _SET_BOTTOM_ODE_(self%id_f,-fsdf)
         end if
         end if

      _HORIZONTAL_LOOP_END_
   end subroutine

end module
