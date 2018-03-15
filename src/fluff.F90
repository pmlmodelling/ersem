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
      type (type_horizontal_diagnostic_variable_id) :: id_Q5inc
      real(rk) :: sdQ5
      integer :: ndeposition
      real(rk),allocatable :: qxc(:),qxn(:),qxp(:),qxs(:),qxf(:)
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self,configunit)
      class (type_ersem_fluff), intent(inout), target :: self
      integer,                               intent(in)            :: configunit
      integer :: idep
      character(len=16) :: num
         character(len=10) :: composition
         call self%type_ersem_benthic_base%initialize(configunit)

         call self%get_parameter(self%sdQ5,'sdQ5','d-1','rate of fluff incorporation into sediments')

         call self%get_parameter(self%ndeposition,'ndeposition','', 'number of deposition targets',default=0)
         allocate(self%qxc(self%ndeposition))
         allocate(self%qxn(self%ndeposition))
         allocate(self%qxp(self%ndeposition))
         allocate(self%qxs(self%ndeposition))
       if (use_iron) allocate(self%qxf(self%ndeposition))

       do idep = 2,self%ndeposition
         write(num,'(i0)') idep
          call self%get_parameter(self%qxc(idep),'qxc'//trim(num),'-','fraction of carbon sedimentation into deposition target '//trim(num))
          call self%get_parameter(self%qxn(idep),'qxn'//trim(num),'-','fraction of nitrogen sedimentation into deposition target '//trim(num))
          call self%get_parameter(self%qxp(idep),'qxp'//trim(num),'-','fraction of phosphorus sedimentation into deposition target '//trim(num))
          call self%get_parameter(self%qxs(idep),'qxs'//trim(num),'-','fraction of silicate sedimentation into deposition target '//trim(num))
          if (use_iron) call self%get_parameter(self%qxf(idep),'qxf'//trim(num),'-','fraction of iron sedimentation into deposition target '//trim(num))
      end do

          self%qxc(1)=1._rk-sum(self%qxc(2:))
          self%qxn(1)=1._rk-sum(self%qxn(2:))
          self%qxp(1)=1._rk-sum(self%qxp(2:))
          self%qxs(1)=1._rk-sum(self%qxs(2:))
       if (use_iron) self%qxf(1)=1._rk-sum(self%qxf(2:))

         allocate(self%id_target(self%ndeposition))
         allocate(self%id_targetc(self%ndeposition))
         allocate(self%id_targetn(self%ndeposition))
         allocate(self%id_targetp(self%ndeposition))
         allocate(self%id_targets(self%ndeposition))
      if (use_iron) allocate(self%id_targetf(self%ndeposition))

      do idep=1,self%ndeposition
        write(num,'(i0)') idep
        call self%register_model_dependency(self%id_target(idep),'deposition_target'//trim(num))
        if (self%qxc(idep)/=0) then
         call self%register_state_dependency(self%id_targetc(idep),'deposition_target'//trim(num)//'c','mg C/m^2','deposition_target '//trim(num)//' carbon')
         call self%request_coupling_to_model(self%id_targetc(idep),self%id_target(idep),'c')
         call self%register_diagnostic_variable(self%id_Q5inc,'Q5inc','mg C/m^2/d','Q5 to Q6',domain=domain_bottom,source=source_do_bottom)

        end if
        if (self%qxn(idep)/=0) then
         call self%register_state_dependency(self%id_targetn(idep),'deposition_target'//trim(num)//'n','mmol N/m^2','deposition_target '//trim(num)//' nitrogen')
         call self%request_coupling_to_model(self%id_targetn(idep),self%id_target(idep),'n')
        end if
        if (self%qxp(idep)/=0) then
         call self%register_state_dependency(self%id_targetp(idep),'deposition_target'//trim(num)//'p','mmol P/m^2','deposition_target '//trim(num)//' phosphorus')
         call self%request_coupling_to_model(self%id_targetp(idep),self%id_target(idep),'p')
        end if
        if (self%qxs(idep)/=0) then
         call self%register_state_dependency(self%id_targets(idep),'deposition_target'//trim(num)//'s','mmol Si/m^2','deposition_target '//trim(num)//' silicate')
         call self%request_coupling_to_model(self%id_targets(idep),self%id_target(idep),'s')
        end if
        if (use_iron) then
        if (self%qxf(idep)/=0) then
         call self%register_state_dependency(self%id_targetf(idep),'deposition_target'//trim(num)//'f','mmol Fe/m^2','deposition_target '//trim(num)//' iron')
         call self%request_coupling_to_model(self%id_targetf(idep),self%id_target(idep),'f')
        end if
        end if

      end do

         call self%get_parameter(composition,'composition','','elemental composition')

   end subroutine initialize

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)

      class (type_ersem_fluff), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: fsdc,fsdp,fsdn,fsds,fsdf
      real(rk) :: Pc,Pn,Pp,Ps,Pf
      integer :: idep
      call self%type_ersem_benthic_base%do_bottom(_ARGUMENTS_DO_BOTTOM_)

      _HORIZONTAL_LOOP_BEGIN_

         Pc = 0.0_rk
         Pn = 0.0_rk
         Pp = 0.0_rk
         Ps = 0.0_rk
         Pf = 0.0_rk

         if (_VARIABLE_REGISTERED_(self%id_c)) then
            _GET_HORIZONTAL_(self%id_c,Pc)
            fsdc = self%sdQ5*Pc
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Q5inc,fsdc)
            do idep=1,self%ndeposition
              if (self%qxc(idep)/=0) then
              _SET_BOTTOM_ODE_(self%id_targetc(idep), fsdc*self%qxc(idep))
              end if
            end do
            _SET_BOTTOM_ODE_(self%id_c,-fsdc)
         end if

         if (_VARIABLE_REGISTERED_(self%id_n)) then
            _GET_HORIZONTAL_(self%id_n,Pn)
            fsdn = self%sdQ5*Pn
            do idep=1,self%ndeposition
              if (self%qxn(idep)/=0) then
              _SET_BOTTOM_ODE_(self%id_targetn(idep), fsdn*self%qxn(idep))
              end if
            end do
            _SET_BOTTOM_ODE_(self%id_n,-fsdn)
         end if

         if (_VARIABLE_REGISTERED_(self%id_p)) then
            _GET_HORIZONTAL_(self%id_p,Pp)
            fsdp = self%sdQ5*Pp
            do idep=1,self%ndeposition
              if (self%qxp(idep)/=0) then
              _SET_BOTTOM_ODE_(self%id_targetp(idep), fsdp*self%qxp(idep))
              end if
            end do
            _SET_BOTTOM_ODE_(self%id_p,-fsdp)
         end if

         if (_VARIABLE_REGISTERED_(self%id_s)) then
            _GET_HORIZONTAL_(self%id_s,Ps)
            fsds = self%sdQ5*Ps
            do idep=1,self%ndeposition
              if (self%qxs(idep)/=0) then
              _SET_BOTTOM_ODE_(self%id_targets(idep), fsds*self%qxs(idep))
              end if
            end do
            _SET_BOTTOM_ODE_(self%id_s,-fsds)
         end if

         if (use_iron) then
         if (_VARIABLE_REGISTERED_(self%id_f)) then
            _GET_HORIZONTAL_(self%id_f,Pf)
            fsdf = self%sdQ5*Pf
            do idep=1,self%ndeposition
              if (self%qxf(idep)/=0) then
              _SET_BOTTOM_ODE_(self%id_targetf(idep), fsdf*self%qxf(idep))
              end if
            end do
            _SET_BOTTOM_ODE_(self%id_f,-fsdf)
         end if
         end if

      _HORIZONTAL_LOOP_END_
   end subroutine

end module
