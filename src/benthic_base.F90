#include "fabm_driver.h"

! Benthic variable that supports resuspension and remineralization.
! Both processes return material to the pelagic.

module ersem_benthic_base

   use fabm_types
   use fabm_particle

   use ersem_shared

   implicit none

!  default: all is private.
   private

   type,extends(type_particle_model),public :: type_ersem_benthic_base
      ! Own state variables
      type (type_bottom_state_variable_id) :: id_c,id_n,id_p,id_f,id_s

      ! Coupled state variables for resuspension and remineralization
      type (type_state_variable_id) :: id_resuspended_c,id_resuspended_n,id_resuspended_p,id_resuspended_s,id_resuspended_f
      type (type_horizontal_diagnostic_variable_id) :: id_resuspension_flux_c,id_resuspension_flux_n,id_resuspension_flux_p,id_resuspension_flux_s,id_resuspension_flux_f,id_bremin
      type (type_state_variable_id) :: id_O3c,id_N1p,id_N3n,id_N4n,id_N5s,id_N7f,id_TA

      ! Dependencies for resuspension
      type (type_horizontal_dependency_id) :: id_bedstress
      type (type_dependency_id)            :: id_dens, id_ETW

      ! Parameters
      real(rk) :: reminQIX,pQIN3X
      logical  :: resuspension
      real(rk) :: er, vel_crit, q10
   contains
      procedure :: initialize
      procedure :: do_bottom

      procedure :: initialize_ersem_benthic_base
      procedure :: add_constituent => benthic_base_add_constituent
   end type

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_ersem_benthic_base), intent(inout), target :: self
   integer,                         intent(in)            :: configunit

   character(len=10) :: composition
   real(rk) :: c0
!
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%initialize_ersem_benthic_base()

      call self%get_parameter(composition, 'composition', '', 'elemental composition')
      call self%get_parameter(self%reminQIX, 'remin', '1/d','remineralisation rate at 10 degrees Celsius',default=0.0_rk)
      if (self%reminQIX /= 0.0_rk) then
         call self%get_parameter(self%q10, 'q10', '-', 'Q_10 temperature coefficient', default=1.0_rk, minimum=1.0_rk)
         if (index(composition,'n') /= 0) &
            call self%get_parameter(self%pQIN3X, 'pN3', '-', 'nitrate fraction of remineralised nitrogen (remainder is ammonium)', default=0.0_rk)
         call self%register_dependency(self%id_ETW, standard_variables%temperature)
      end if
      call self%get_parameter(self%resuspension, 'resuspension', '', 'enable resuspension', default=.false.)

      if (self%resuspension) then
         call self%get_parameter(self%er,'er','1/d','erosion rate',default=0.225_rk)
         call self%get_parameter(self%vel_crit,'vel_crit','m/s','critical shear velocity for resuspension',default=0.02_rk)
         call self%register_dependency(self%id_bedstress,standard_variables%bottom_stress)
         call self%register_dependency(self%id_dens,     standard_variables%density)
      end if

      call self%get_parameter(c0,'c0','mg C/m^2','background carbon concentration',default=0.0_rk)

      if (index(composition,'c')/=0) then
         call self%add_constituent('c',0.0_rk,c0)
         if (self%reminQIX/=0.0_rk) call self%register_state_dependency(self%id_O3c,'O3c','mmol/m^3','dissolved inorganic carbon')
         if (self%reminQIX/=0.0_rk) call self%register_diagnostic_variable(self%id_bremin,'bremin','mg C/m^2/d','carbon remineralization',source=source_do_bottom)
      end if
      if (index(composition,'p')/=0) then
         call self%add_constituent('p',0.0_rk,qpRPIcX*c0)
         if (self%reminQIX/=0.0_rk) then
            call self%register_state_dependency(self%id_N1p,'N1p','mmol P/m^3','phosphate')
            if (.not._VARIABLE_REGISTERED_(self%id_TA)) call self%register_state_dependency(self%id_TA,standard_variables%alkalinity_expressed_as_mole_equivalent)
         end if
      end if
      if (index(composition,'n')/=0) then
         call self%add_constituent('n',0.0_rk,qnRPIcX*c0)
         if (self%reminQIX/=0.0_rk) then
            call self%register_state_dependency(self%id_N3n,'N3n','mmol N/m^3','nitrate')
            call self%register_state_dependency(self%id_N4n,'N4n','mmol N/m^3','ammonium')
            if (.not._VARIABLE_REGISTERED_(self%id_TA)) call self%register_state_dependency(self%id_TA,standard_variables%alkalinity_expressed_as_mole_equivalent)
         end if
      end if
      if (index(composition,'s')/=0) then
         call self%add_constituent('s',0.0_rk,qsRPIcX*c0)
         if (self%reminQIX/=0.0_rk) call self%register_state_dependency(self%id_N5s,'N5s','mmol Si/m^3','silicate')
      end if
      if (index(composition,'f')/=0) then
         call self%add_constituent('f',0.0_rk)
         if (use_iron.and.self%reminQIX/=0.0_rk) call self%register_state_dependency(self%id_N7f,'N7f','umol Fe/m^3','dissolved iron')
      end if
   end subroutine

   subroutine initialize_ersem_benthic_base(self)
      class (type_ersem_benthic_base), intent(inout), target :: self

      self%resuspension = .false.

      ! Set time unit to d-1
      ! This implies that all rates (sink/source terms, vertical velocities) are given in d-1.
      self%dt = 86400._rk
   end subroutine

   subroutine benthic_base_add_constituent(self,name,initial_value,background_value,qn,qp)
      class (type_ersem_benthic_base), intent(inout), target :: self
      character(len=*),                intent(in)            :: name
      real(rk),                        intent(in)            :: initial_value
      real(rk),optional,               intent(in)            :: background_value,qn,qp

      select case (name)
         case ('c')
            call register(self%id_c, self%id_resuspended_c, self%id_resuspension_flux_c, 'c', 'mg C', 'carbon', standard_variables%total_carbon, scale_factor=1._rk/CMass)
            if (present(qn)) call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_c, scale_factor=qn)
            if (present(qp)) call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, scale_factor=qp)
         case ('n')
            call register(self%id_n, self%id_resuspended_n, self%id_resuspension_flux_n, 'n', 'mmol N', 'nitrogen', standard_variables%total_nitrogen)
         case ('p')
            call register(self%id_p, self%id_resuspended_p, self%id_resuspension_flux_p, 'p', 'mmol P', 'phosphorus', standard_variables%total_phosphorus)
         case ('s')
            call register(self%id_s, self%id_resuspended_s, self%id_resuspension_flux_s, 's', 'mmol Si', 'silicate', standard_variables%total_silicate)
         case ('f')
            if (use_iron) call register(self%id_f, self%id_resuspended_f, self%id_resuspension_flux_f, 'f', 'umol Fe', 'iron', standard_variables%total_iron)
         case default
            call self%fatal_error('benthic_base_add_constituent','Unknown constituent "'//trim(name)//'".')
      end select

   contains

      subroutine register(variable_id, resuspended_id, resuspended_flux_id, name, base_units, long_name, aggregate_variable, scale_factor)
         type (type_bottom_state_variable_id),          intent(inout), target :: variable_id
         type (type_state_variable_id),                 intent(inout), target :: resuspended_id
         type (type_horizontal_diagnostic_variable_id), intent(inout), target :: resuspended_flux_id
         character(len=*),                              intent(in)            :: name,base_units,long_name
         type (type_bulk_standard_variable),            intent(in)            :: aggregate_variable
         real(rk),                                      intent(in), optional  :: scale_factor

         call self%register_state_variable(variable_id, name, base_units//'/m^2', long_name, initial_value, minimum=0._rk,background_value=background_value)
         call self%add_to_aggregate_variable(aggregate_variable, variable_id, scale_factor=scale_factor)
         if (self%resuspension) then
            call self%register_state_dependency(resuspended_id, 'resuspended_'//name, base_units//'/m^3', 'pelagic variable taking up resuspended '//long_name)
            call self%request_coupling_to_model(resuspended_id, 'RP', name)
            call self%register_diagnostic_variable(resuspended_flux_id, 'resuspension_flux_'//name, base_units//'/m^2/d', 'resuspension of '//long_name, source=source_do_bottom)
         end if
      end subroutine

   end subroutine benthic_base_add_constituent

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_benthic_base), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: bedstress,density
      real(rk) :: fac,state,resuspension_flux
      real(rk) :: max_rel_res = 4.0_rk            ! max relative loss of matter due to resuspension, in 1/d
      real(rk) :: eT,ETW,reminT

      _HORIZONTAL_LOOP_BEGIN_
         if (self%resuspension) then
            ! The resuspension rate (1/d) is a linear function of shear stress (Pa)
            ! Prefactor "er" (1/d) can be interpreted as c*M/rho_sed*v_crit^2, with:
            ! - c (1/m) = ratio between tracer concentration at the sediment surface and sediemnt-column-integrated tracer
            !   (e.g., c=100 for an exponential profile with penetration depth of 1 cm)
            ! - M = the erosion rate in g*s/m4 (Puls & Suendermann 1990: M=100)
            ! - rho_sed = dry mass of sediment per total volume at the sediment surface. This is grain density (2650 kg/m3 for quartz) multiplied by (1-porosity)
            ! - v_crit being the critical shear velocity (m/s) for resuspension (Puls & Suendermann 1990: v_crit=0.02)
            ! With porosity=0.4, a representative value for er is 100*100/(2650000*0.6)*0.02^2*86400 = 0.225 1/d
            ! Note that the square of bed shear velocity is calculated as the ratio between shear stress (Pa) and water density (kg/m^3).
            _GET_HORIZONTAL_(self%id_bedstress,bedstress)
            _GET_(self%id_dens,density)
            fac = min(max_rel_res,self%er*max(0.0_rk,bedstress/density/self%vel_crit**2 - 1._rk))

            ! Carbon
            if (_VARIABLE_REGISTERED_(self%id_c)) then
               _GET_HORIZONTAL_(self%id_c,state)
               resuspension_flux = fac*state
               _SET_BOTTOM_ODE_(self%id_c,-resuspension_flux)
               _SET_BOTTOM_EXCHANGE_(self%id_resuspended_c,resuspension_flux)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_resuspension_flux_c,resuspension_flux)
            end if

            ! Nitrogen
            if (_VARIABLE_REGISTERED_(self%id_n)) then
               _GET_HORIZONTAL_(self%id_n,state)
               resuspension_flux = fac*state
               _SET_BOTTOM_ODE_(self%id_n,-resuspension_flux)
               _SET_BOTTOM_EXCHANGE_(self%id_resuspended_n,resuspension_flux)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_resuspension_flux_n,resuspension_flux)
            end if

            ! Phosphorus
            if (_VARIABLE_REGISTERED_(self%id_p)) then
               _GET_HORIZONTAL_(self%id_p,state)
               resuspension_flux = fac*state
               _SET_BOTTOM_ODE_(self%id_p,-resuspension_flux)
               _SET_BOTTOM_EXCHANGE_(self%id_resuspended_p,resuspension_flux)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_resuspension_flux_p,resuspension_flux)
            end if

            ! Silicate
            if (_VARIABLE_REGISTERED_(self%id_s)) then
               _GET_HORIZONTAL_(self%id_s,state)
               resuspension_flux = fac*state
               _SET_BOTTOM_ODE_(self%id_s,-resuspension_flux)
               _SET_BOTTOM_EXCHANGE_(self%id_resuspended_s,resuspension_flux)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_resuspension_flux_s,resuspension_flux)
            end if

            ! Iron
            if (_VARIABLE_REGISTERED_(self%id_f)) then
               _GET_HORIZONTAL_(self%id_f,state)
               resuspension_flux = fac*state
               _SET_BOTTOM_ODE_(self%id_f,-resuspension_flux)
               _SET_BOTTOM_EXCHANGE_(self%id_resuspended_f,resuspension_flux)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_resuspension_flux_f,resuspension_flux)
            end if
         end if   ! Resuspension

         ! Remineralization (benthic return)
         if (self%reminQIX/=0.0_rk) then
            _GET_(self%id_ETW,ETW)
            eT = self%q10**((ETW-10._rk)/10._rk)
            reminT=self%reminQIX * eT

            if (_VARIABLE_REGISTERED_(self%id_c)) then
               _GET_HORIZONTAL_(self%id_c,state)
               _SET_BOTTOM_ODE_(self%id_c,-reminT*state)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_bremin,reminT*state)
               _SET_BOTTOM_EXCHANGE_(self%id_O3c,reminT*state/CMass)
            end if
            if (_VARIABLE_REGISTERED_(self%id_p)) then
               _GET_HORIZONTAL_(self%id_p,state)
               _SET_BOTTOM_ODE_(self%id_p,-reminT*state)
               _SET_BOTTOM_EXCHANGE_(self%id_N1p,reminT*state)
               _SET_BOTTOM_EXCHANGE_(self%id_TA, -reminT*state)
            end if
            if (_VARIABLE_REGISTERED_(self%id_n)) then
               _GET_HORIZONTAL_(self%id_n,state)
               _SET_BOTTOM_ODE_(self%id_n,-reminT*state)
               _SET_BOTTOM_EXCHANGE_(self%id_N3n,self%pQIN3X*reminT*state)
               _SET_BOTTOM_EXCHANGE_(self%id_N4n,(1.0_rk-self%pQIN3X)*reminT*state)
               _SET_BOTTOM_EXCHANGE_(self%id_TA,(1.0_rk-self%pQIN3X)*reminT*state - self%pQIN3X*reminT*state)
            end if
            if (_VARIABLE_REGISTERED_(self%id_s)) then
               _GET_HORIZONTAL_(self%id_s,state)
               _SET_BOTTOM_ODE_(self%id_s,-reminT*state)
               _SET_BOTTOM_EXCHANGE_(self%id_N5s,reminT*state)
            end if
            if (_VARIABLE_REGISTERED_(self%id_f)) then
               _GET_HORIZONTAL_(self%id_f,state)
               _SET_BOTTOM_ODE_(self%id_f,-reminT*state)
               _SET_BOTTOM_EXCHANGE_(self%id_N7f,reminT*state)
            end if
         end if   ! Remineralization (benthic return)

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module
