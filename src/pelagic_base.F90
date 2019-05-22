#include "fabm_driver.h"

module ersem_pelagic_base

   use fabm_types
   use fabm_particle

   use ersem_shared

   implicit none

!  default: all is private.
   private

   type,extends(type_particle_model),public :: type_ersem_pelagic_base
      type (type_state_variable_id)                 :: id_c,id_n,id_p,id_f,id_s,id_chl
      type (type_horizontal_dependency_id)          :: id_bedstress,id_wdepth
      type (type_dependency_id)                     :: id_dens
      type (type_horizontal_diagnostic_variable_id) :: id_w_bot
      type (type_horizontal_diagnostic_variable_id),allocatable,dimension(:) :: id_cdep,id_ndep,id_pdep,id_sdep,id_fdep

      ! Target variables for sedimentation
      type (type_bottom_state_variable_id),allocatable,dimension(:) :: id_targetc,id_targetn,id_targetp,id_targets,id_targetf

      real(rk) :: rm = 0.0_rk
      real(rk) :: tdep
      integer :: ndeposition
      logical :: no_river_dilution = .false.
      real(rk),allocatable :: qxc(:),qxn(:),qxp(:),qxs(:),qxf(:)
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
      integer,                         intent(in)            :: configunit

      character(len=10) :: composition
      real(rk)          :: c0,s0,rRPmX,EPS,iopABS,iopBBS
      real(rk)          :: qn, qp

      call self%get_parameter(composition, 'composition', '', 'elemental composition')
      call self%get_parameter(c0, 'c0', 'mg C/m^3', 'background carbon concentration', default=0.0_rk, minimum=0.0_rk)

      if (index(composition,'c')/=0) then
         ! Carbon-based light attenuation [optional, off by default]
         call self%get_parameter(EPS,   'EPS',    'm^2/mg C', 'specific shortwave attenuation', default=0.0_rk, minimum=0.0_rk)
         call self%get_parameter(iopABS,'iopABS', 'm^2/mg C', 'specific shortwave absorption',  default=0.0_rk, minimum=0.0_rk)
         call self%get_parameter(iopBBS,'iopBBS', 'm^2/mg C', 'specific shortwave backscatter', default=0.0_rk, minimum=0.0_rk)

         ! Constant N:C stoichiometry [optional, off by default, and only supported if no explicit nitrogen constituent is present]
         if (index(composition,'n') == 0) then
            call self%get_parameter(qn, 'qn', 'mmol N/mg C', 'nitrogen : carbon ratio', default=0.0_rk, minimum=0.0_rk)
         else
            qn = 0.0_rk
         end if

         ! Constant P:C stoichiometry [optional, off by default, and only supported if no explicit phosphorus constituent is present]
         if (index(composition,'p') == 0) then
            call self%get_parameter(qp, 'qp', 'mmol P/mg C', 'phosphorus : carbon ratio', default=0.0_rk, minimum=0.0_rk)
         else
            qp = 0.0_rk
         end if
      end if
      call self%get_parameter(rRPmX, 'rm', 'm/d', 'sinking velocity', default=0.0_rk)

      call self%initialize_ersem_base(rm=rRPmX, sedimentation=rRPmX>0._rk)

      if (index(composition,'c')/=0) then
         call self%add_constituent('c', 0.0_rk, c0, qn, qp)

         ! Add contributions to light attenuation, absorption, scattering.
         ! Contributions with a scale_factor of 0.0 will automatically be ignored.
         call self%add_to_aggregate_variable(particulate_organic_absorption_coefficient, &
            self%id_c,scale_factor=iopABS,include_background=.true.)
         call self%add_to_aggregate_variable(particulate_organic_backscatter_coefficient, &
            self%id_c,scale_factor=iopBBS,include_background=.true.)
         call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &
            self%id_c,scale_factor=EPS,include_background=.true.)
      end if
      if (index(composition,'n')/=0) call self%add_constituent('n',0.0_rk,qnRPIcX*c0)
      if (index(composition,'p')/=0) call self%add_constituent('p',0.0_rk,qpRPIcX*c0)
      if (index(composition,'s')/=0) then
         call self%get_parameter(s0,'s0','mmol Si/m^3','background silicon concentration',default=qsRPIcX*c0)
         call self%add_constituent('s',0.0_rk,s0)
      end if
      if (index(composition,'f')/=0) call self%add_constituent('f',0.0_rk)

   end subroutine

   subroutine initialize_ersem_base(self,rm,sedimentation)
      class (type_ersem_pelagic_base), intent(inout), target :: self
      real(rk),optional,               intent(in)            :: rm
      logical, optional,               intent(in)            :: sedimentation

      real(rk) :: vel_crit

      ! We are adding a new yaml entry for each ersem_base type related to river dilution behaviour
      call self%get_parameter(self%no_river_dilution,'no_river_dilution','','disable river dilution by setting riverine concentrations equal to those in the receiving grid cell',default=.false.)

      ! Set time unit to d-1
      ! This implies that all rates (sink/source terms, vertical velocities) are given in d-1.
      self%dt = 86400._rk

      ! Sedimentation
      self%ndeposition = 0
      if (present(sedimentation)) then
         ! Explicitly specified whether this particle can sink (and thus be deposited at bed).
         ! If so, allow the user to specify number of pools to deposit in (0 means no deposition)
         ! If not, do not offer option to configure deposition at all.
         if (sedimentation) call self%get_parameter(self%ndeposition,'ndeposition','', 'number of target pools for sedimentation',default=1)
      else
         ! No sinking behaviour specified; get number of pools to deposit in (default: no deposition).
         call self%get_parameter(self%ndeposition,'ndeposition','', 'number of target pools for sedimentation',default=0)
      end if

      if (self%ndeposition>0) then
         call self%get_parameter(vel_crit,'vel_crit','m/s','critical bed shear velocity for deposition',default=0.01_rk)
         self%tdep = vel_crit**2

         call self%register_dependency(self%id_bedstress,standard_variables%bottom_stress)
         call self%register_dependency(self%id_dens,     standard_variables%density)
      end if

      ! Vertical velocity (positive: downwards, negative: upwards)
      if (present(rm)) self%rm = rm

      call self%register_diagnostic_variable(self%id_w_bot,'w_bot','m/s','near-bed vertical velocity',output=output_none,source=source_do_bottom)

   end subroutine initialize_ersem_base

   subroutine add_constituent(self,name,initial_value,background_value,qn,qp)
      class (type_ersem_pelagic_base), intent(inout), target :: self
      character(len=*),                intent(in)            :: name
      real(rk),                        intent(in)            :: initial_value
      real(rk),optional,               intent(in)            :: background_value,qn,qp

      select case (name)
      case ('c')
         call register(self%id_c,'c','mg C','carbon',standard_variables%total_carbon,self%qxc,self%id_cdep,self%id_targetc,1._rk/CMass)
         if (present(qn)) call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_c,scale_factor=qn)
         if (present(qp)) call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_c,scale_factor=qp)
      case ('n')
         call register(self%id_n,'n','mmol N','nitrogen',standard_variables%total_nitrogen,self%qxn,self%id_ndep,self%id_targetn)
      case ('p')
         call register(self%id_p,'p','mmol P','phosphorus',standard_variables%total_phosphorus,self%qxp,self%id_pdep,self%id_targetp)
      case ('s')
         call register(self%id_s,'s','mmol Si','silicate',standard_variables%total_silicate,self%qxs,self%id_sdep,self%id_targets)
      case ('f')
         if (use_iron) call register(self%id_f,'f','umol Fe','iron',standard_variables%total_iron,self%qxf,self%id_fdep,self%id_targetf)
      case ('chl')
         call register(self%id_chl,'Chl','mg','chlorophyll a',total_chlorophyll)
      case default
         call self%fatal_error('add_constituent','Unknown constituent "'//trim(name)//'".')
      end select

   contains

      subroutine register(variable_id,name,base_units,long_name,aggregate_variable,qx,id_xdep,id_dep,scale_factor)
         type (type_state_variable_id),       intent(inout), target     :: variable_id
         character(len=*),                    intent(in)                :: name,base_units,long_name
         type (type_bulk_standard_variable),  intent(in)                :: aggregate_variable
         real(rk),                            intent(inout),allocatable,optional :: qx(:)
         type (type_bottom_state_variable_id),intent(inout),allocatable,optional :: id_dep(:)
         type (type_horizontal_diagnostic_variable_id),intent(inout),allocatable,optional :: id_xdep(:)
         real(rk),                            intent(in), optional      :: scale_factor

         integer :: idep
         character(len=16) :: num

         ! Register state variable
         call self%register_state_variable(variable_id,name,trim(base_units)//'/m^3',long_name, &
            initial_value,minimum=0._rk,vertical_movement=-self%rm/self%dt,background_value=background_value,no_river_dilution=self%no_river_dilution)

         ! Contribute to aggregate variable
         call self%add_to_aggregate_variable(aggregate_variable,variable_id,scale_factor)

         if (self%ndeposition>0.and.present(qx)) then
            ! Create array with fractions per depostion target
            allocate(qx(self%ndeposition))
            do idep = 2,self%ndeposition
               write(num,'(i0)') idep
               call self%get_parameter(qx(idep),'qx'//trim(name)//trim(num),'-','fraction of '//trim(long_name)//' sinking into deposition target '//trim(num))
            end do
            qx(1)=1._rk-sum(qx(2:))

            ! Link to pools in which to deposit matter.
            allocate(id_dep(self%ndeposition))
            allocate(id_xdep(self%ndeposition))

            do idep=1,self%ndeposition
               write(num,'(i0)') idep
               if (qx(idep)/=0) then
                  call self%register_state_dependency(id_dep(idep),'deposition_target'//trim(num)//trim(name),trim(base_units)//'/m^2','target pool '//trim(num)//' for '//trim(long_name)//' deposition')
                  call self%request_coupling_to_model(id_dep(idep),'deposition_target'//trim(num),name)
                  call self%register_diagnostic_variable(id_xdep(idep),'dep'//trim(num)//trim(name),trim(base_units)//'/m^2/d','deposition of '//trim(long_name)//' in target pool '//trim(num),source=source_do_bottom)
               end if
            end do
         end if
      end subroutine register
   end subroutine add_constituent

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

      real(rk) :: tbed,density
      real(rk) :: sdrate
      real(rk) :: w,conc,flux
      integer :: idep

      _HORIZONTAL_LOOP_BEGIN_

         ! Retrieve sinking rate (at centre of cell closest to the bottom)
         ! NB by ERSEM convention, this rate is returned in m/d, positive for sinking, negative for floating!
         w = self%get_sinking_rate(_ARGUMENTS_LOCAL_)

         ! Store near-bed vertical velocity. Use FABM convention (m/s, negative for downward) to
         ! allow this custom fuctionality to be ultimately be replaced by a FABM API.
         ! The near-bed vertical velocity may be used by other modules to compute e.g. near-bed Rouse profiles.
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_w_bot,-w/86400)

         if (self%ndeposition/=0) then

         ! Retrieve bed stress and local density - needed to determine sedimentation rate from sinking rate.
         _GET_HORIZONTAL_(self%id_bedstress,tbed)
         _GET_(self%id_dens,density)

         ! Divide actual bed stress (Pa) by density (kg/m^3) to obtain square of bed shear velocity.
         tbed = tbed/density

         ! Deposition rate based on sinking velocity but mediated by bottom stress (high stress = no deposition)
         ! Original reference: Puls and Suendermann 1990. However, they use the ratio of shear velocities, while we use its square.
         !sdrate = min(fsd*fac,pdepth(I)/timestep) ! Jorn: CFL criterion disabled because FABM does not provide timestep
         sdrate = w * max(0._rk, 1._rk - tbed/self%tdep)

         if (_AVAILABLE_(self%id_c)) then
            _GET_(self%id_c,conc)
            flux = sdrate*conc
            do idep=1,self%ndeposition
              if (self%qxc(idep)/=0) then
                _SET_BOTTOM_ODE_(self%id_targetc(idep), flux*self%qxc(idep))
                _SET_HORIZONTAL_DIAGNOSTIC_(self%id_cdep(idep),flux*self%qxc(idep))
              end if
            end do
            _SET_BOTTOM_EXCHANGE_(self%id_c,-flux)
         end if

         if (_AVAILABLE_(self%id_n)) then
            _GET_(self%id_n,conc)
            flux = sdrate*conc
            do idep=1,self%ndeposition
              if (self%qxn(idep)/=0) then
                _SET_BOTTOM_ODE_(self%id_targetn(idep), flux*self%qxn(idep))
                _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ndep(idep),flux*self%qxn(idep))
              end if
            end do
            _SET_BOTTOM_EXCHANGE_(self%id_n,-flux)
         end if

         if (_AVAILABLE_(self%id_p)) then
            _GET_(self%id_p,conc)
            flux = sdrate*conc
            do idep=1,self%ndeposition
              if (self%qxp(idep)/=0) then
                _SET_BOTTOM_ODE_(self%id_targetp(idep), flux*self%qxp(idep))
                _SET_HORIZONTAL_DIAGNOSTIC_(self%id_pdep(idep),flux*self%qxp(idep))
              end if
            end do
            _SET_BOTTOM_EXCHANGE_(self%id_p,-flux)
         end if

         if (_AVAILABLE_(self%id_s)) then
            _GET_(self%id_s,conc)
            flux = sdrate*conc
            do idep=1,self%ndeposition
               if (self%qxs(idep)/=0) then
                 _SET_BOTTOM_ODE_(self%id_targets(idep), flux*self%qxs(idep))
                 _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sdep(idep),flux*self%qxs(idep))
               end if
            end do
            _SET_BOTTOM_EXCHANGE_(self%id_s,-flux)
         end if

         if (use_iron) then
            if (_AVAILABLE_(self%id_f)) then
               _GET_(self%id_f,conc)
               flux = sdrate*conc
               do idep=1,self%ndeposition
                  if (self%qxf(idep)/=0) then
                    _SET_BOTTOM_ODE_(self%id_targetf(idep), flux*self%qxf(idep))
                    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fdep(idep),flux*self%qxf(idep))
                  end if
               end do
               _SET_BOTTOM_EXCHANGE_(self%id_f,-flux)
            end if
         end if

         if (_AVAILABLE_(self%id_chl)) then
            ! Chlorophyll is simply lost by sinking:
            _GET_(self%id_chl,conc)
            flux = sdrate*conc
            _SET_BOTTOM_EXCHANGE_(self%id_chl,-flux)
         end if

         end if

      _HORIZONTAL_LOOP_END_
   end subroutine

end module
