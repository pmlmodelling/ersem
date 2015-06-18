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
      type (type_bottom_state_variable_id) :: id_c,id_n,id_p,id_f,id_s,id_l

      ! Coupled state variables for resuspension and remineralization
      type (type_state_variable_id) :: id_resuspension_c,id_resuspension_n,id_resuspension_p,id_resuspension_s,id_resuspension_f,id_resuspension_l
      type (type_state_variable_id) :: id_O3c,id_N1p,id_N3n,id_N4n,id_N5s,id_N7f
      type (type_model_id)          :: id_resuspension_target

      ! Dependencies for resuspension
      type (type_horizontal_dependency_id) :: id_bedstress,id_wdepth
      type (type_dependency_id)            :: id_dens

      ! Parameters
      character(len=10) :: composition
      real(rk) :: reminQIX,pQIN3X
      logical  :: resuspension
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
!

!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%composition, 'composition', '',   'elemental composition')
      call self%get_parameter(self%reminQIX,    'remin',       '1/d','remineralisation rate',default=0.0_rk)
      if (index(self%composition,'n')/=0 .and. self%reminQIX/=0.0_rk) &
         call self%get_parameter(self%pQIN3X,'pN3','-','nitrate fraction of remineralised nitrogen (remainder is ammonium)',default=0.0_rk)
      call self%get_parameter(self%resuspension,'resuspension','',   'enable resuspension',  default=.false.)

      call self%initialize_ersem_benthic_base()

      if (self%resuspension) then
         call self%register_model_dependency(self%id_resuspension_target,'RP')
         call self%register_dependency(self%id_bedstress,standard_variables%bottom_stress)
         call self%register_dependency(self%id_wdepth,   standard_variables%bottom_depth_below_geoid)
         call self%register_dependency(self%id_dens,     standard_variables%density)
      end if

      if (index(self%composition,'c')/=0) then
         call self%add_constituent('c',0.0_rk)
         if (self%resuspension) then
            call self%register_state_dependency(self%id_resuspension_c,'resuspension_target_c','mg C/m^3','pelagic variable taking up resuspended carbon')
            call self%request_coupling_to_model(self%id_resuspension_c,self%id_resuspension_target,'c')
         end if
         if (self%reminQIX/=0.0_rk) call self%register_state_dependency(self%id_O3c,'O3c','mmol/m^3','dissolved inorganic carbon')
      end if
      if (index(self%composition,'p')/=0) then
         call self%add_constituent('p',0.0_rk)
         if (self%resuspension) then
            call self%register_state_dependency(self%id_resuspension_p,'resuspension_target_p','mmol P/m^3','pelagic variable taking up resuspended phosphorus')
            call self%request_coupling_to_model(self%id_resuspension_p,self%id_resuspension_target,'p')
         end if
         if (self%reminQIX/=0.0_rk) call self%register_state_dependency(self%id_N1p,'N1p','mmol P/m^3','phosphate')
      end if
      if (index(self%composition,'n')/=0) then
         call self%add_constituent('n',0.0_rk)
         if (self%resuspension) then
            call self%register_state_dependency(self%id_resuspension_n,'resuspension_target_n','mmol N/m^3','pelagic variable taking up resuspended nitrogen')
            call self%request_coupling_to_model(self%id_resuspension_n,self%id_resuspension_target,'n')
         end if
         if (self%reminQIX/=0.0_rk) then
            call self%register_state_dependency(self%id_N3n,'N3n','mmol N/m^3','nitrate')
            call self%register_state_dependency(self%id_N4n,'N4n','mmol N/m^3','ammonium')
         end if
      end if
      if (index(self%composition,'s')/=0) then
         call self%add_constituent('s',0.0_rk)
         if (self%resuspension) then
            call self%register_state_dependency(self%id_resuspension_s,'resuspension_target_s','mmol Si/m^3','pelagic variable taking up resuspended silicate')
            call self%request_coupling_to_model(self%id_resuspension_s,self%id_resuspension_target,'s')
         end if
         if (self%reminQIX/=0.0_rk) call self%register_state_dependency(self%id_N5s,'N5s','mmol Si/m^3','silicate')
      end if
      if (index(self%composition,'f')/=0) then
         call self%add_constituent('f',0.0_rk)
         if (use_iron.and.self%resuspension) then
            call self%register_state_dependency(self%id_resuspension_f,'resuspension_target_f','umol Fe/m^3','pelagic variable taking up resuspended iron')
            call self%request_coupling_to_model(self%id_resuspension_f,self%id_resuspension_target,'f')
         end if
         if (use_iron.and.self%reminQIX/=0.0_rk) call self%register_state_dependency(self%id_N7f,'N7f','umol Fe/m^3','dissolved iron')
      end if
      if (index(self%composition,'l')/=0) then
         call self%add_constituent('l',0.0_rk)
         if (self%resuspension) then
            call self%register_state_dependency(self%id_resuspension_l,'resuspension_target_l','mg/m^3','pelagic variable taking up resuspended calcite')
            call self%request_coupling_to_model(self%id_resuspension_l,self%id_resuspension_target,'l')
         end if
         if (index(self%composition,'c')==0) call self%register_state_dependency(self%id_O3c,'O3c','umol/m^3','dissolved inorganic carbon')
      end if
   end subroutine

   subroutine initialize_ersem_benthic_base(self,c_ini,n_ini,p_ini,s_ini,f_ini)
      class (type_ersem_benthic_base), intent(inout), target :: self
      real(rk),optional,               intent(in)            :: c_ini,n_ini,p_ini,s_ini,f_ini

      ! Set time unit to d-1
      ! This implies that all rates (sink/source terms, vertical velocities) are given in d-1.
      self%dt = 86400._rk

      if (present(c_ini)) call self%add_constituent('c',c_ini)
      if (present(n_ini)) call self%add_constituent('n',n_ini)
      if (present(p_ini)) call self%add_constituent('p',p_ini)
      if (present(s_ini)) call self%add_constituent('s',s_ini)
      if (present(f_ini)) call self%add_constituent('f',f_ini)
   end subroutine

   subroutine benthic_base_add_constituent(self,name,initial_value,qn,qp)
      class (type_ersem_benthic_base), intent(inout), target :: self
      character(len=*),                intent(in)            :: name
      real(rk),                        intent(in)            :: initial_value
      real(rk),optional,               intent(in)            :: qn,qp

      select case (name)
         case ('c')
            call self%register_state_variable(self%id_c,'c','mg C/m^2','carbon',initial_value,minimum=0._rk)
            call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_c,scale_factor=1._rk/CMass)
            if (present(qn)) call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_c,scale_factor=qn)
            if (present(qp)) call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_c,scale_factor=qp)
         case ('n')
            call self%register_state_variable(self%id_n, 'n', 'mmol N/m^2', 'nitrogen', initial_value, minimum=0._rk)
            call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_n)
         case ('p')
            call self%register_state_variable(self%id_p, 'p', 'mmol P/m^2', 'phosphorus', initial_value, minimum=0._rk)
            call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_p)
         case ('s')
            call self%register_state_variable(self%id_s,'s','mmol Si/m^2','silicate',initial_value,minimum=0._rk)
            call self%add_to_aggregate_variable(standard_variables%total_silicate,self%id_s)
         case ('f')
            if (use_iron) then
               call self%register_state_variable(self%id_f,'f','umol Fe/m^2','iron',initial_value,minimum=0._rk)
               call self%add_to_aggregate_variable(standard_variables%total_iron,self%id_f)
            end if
         case ('l')
            call self%register_state_variable(self%id_l,'l','mg C/m^2','calcite',initial_value,minimum=0._rk)
            call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_l,scale_factor=1._rk/CMass)
         case default
            call self%fatal_error('benthic_base_add_constituent','Unknown constituent "'//trim(name)//'".')
      end select
   end subroutine benthic_base_add_constituent

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_benthic_base), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: bedstress,density,wdepth
      real(rk) :: Q6cP,Q6nP,Q6pP,Q6sP,Q6fP,Q6lP
      real(rk) :: bedsedXc,bedsedXn,bedsedXp,bedsedXs
      real(rk) :: fac,FerC,FerN,FerP,FerS,FerF

      ! Bed characteristics - from Puls and Sundermann 1990
      ! Critical shear velocity for erosion = 0.02 m/s
      real(rk),parameter :: ter=0.02_rk**2

      ! erosion constant (g/m^2/s)
      real(rk),parameter :: er=100._rk*Ter * 1000._rk*86400._rk  !convert to mg/day

      _HORIZONTAL_LOOP_BEGIN_
         if (self%resuspension) then
            _GET_HORIZONTAL_(self%id_bedstress,bedstress)
            _GET_(self%id_dens,density)

            ! Divide actual stress (Pa) by density (kg/m^3) to obtain square of bed shear velocity.
            ! This allows for comparison with ter.
            bedstress = bedstress/density
         else
            bedstress = 0.0_rk
         end if

         if (bedstress.gt.ter) then
            _GET_HORIZONTAL_(self%id_wdepth,wdepth)
            _GET_HORIZONTAL_(self%id_c,Q6cP)

            ! Inorganic sedment (could replace by transport model)
            ! for now assume 90% is a fixed inorganic component
            bedsedXc=10._rk*(-1069._rk*LOG(wdepth) + 10900._rk)
            bedsedXn=10._rk*(-7.6368_rk*LOG(wdepth) + 78.564_rk)
            bedsedXp=10._rk*(-0.545_rk*LOG(wdepth) + 6.0114_rk)
            bedsedXs=10._rk*(-64.598_rk*LOG(wdepth) + 391.61_rk)

            ! fac = er*(bedstress(k)/ter - 1.0)
            fac = er*(bedstress/ter - 1._rk)/(Q6cP+bedsedXc)

            ! Carbon
            ! FerC=fac*Q6cP(k)/(Q6cP(k)+bedsedXc)
            FerC=fac*Q6cP
            !FerC=max(min(FerC,Q6cP(k)/timestep+ &
            !       min(SQ6c(k)-wsoQ6c(k),0._rk)),0._rk)
            _SET_BOTTOM_ODE_(self%id_c,-FerC)
            _SET_BOTTOM_EXCHANGE_(self%id_resuspension_c,FerC)

            ! Nitrogen
            if (_VARIABLE_REGISTERED_(self%id_n)) then
               _GET_HORIZONTAL_(self%id_n,Q6nP)
               ! FerN=fac*Q6nP(k)/(Q6nP(k)+bedsedXn)
               FerN=fac*Q6nP
               !FerN=max(min(FerN,Q6nP(k)/timestep+ &
               !       min(SQ6n(k)-wsoQ6n(k),0._rk)),0._rk)
               _SET_BOTTOM_ODE_(self%id_n,-FerN)
               _SET_BOTTOM_EXCHANGE_(self%id_resuspension_n,FerN)
            end if

            ! Phosphorus
            if (_VARIABLE_REGISTERED_(self%id_p)) then
               _GET_HORIZONTAL_(self%id_p,Q6pP)
               ! FerP=fac*Q6pP(k)/(Q6pP(k)+bedsedXp)
               FerP=fac*Q6pP
               !FerP=max(min(FerP,Q6pP(k)/timestep+ &
               !       min(SQ6p(k)-wsoQ6p(k),0._rk)),0._rk)
               _SET_BOTTOM_ODE_(self%id_p,-FerP)
               _SET_BOTTOM_EXCHANGE_(self%id_resuspension_p,FerP)
            end if

            ! Silicate
            if (_VARIABLE_REGISTERED_(self%id_s)) then
               _GET_HORIZONTAL_(self%id_s,Q6sP)
               ! FerS=fac*Q6sP(k)/(Q6sP(k)+bedsedXs)
               FerS=fac*Q6sP
               !FerS=max(min(FerS,Q6sP(k)/timestep+ &
               !       min(SQ6s(k)-wsoQ6s(k),0._rk)),0._rk)
               _SET_BOTTOM_ODE_(self%id_s,-FerS)
               _SET_BOTTOM_EXCHANGE_(self%id_resuspension_s,FerS)
            end if

            ! Iron
            if (_VARIABLE_REGISTERED_(self%id_f)) then
               _GET_HORIZONTAL_(self%id_f,Q6fP)
               FerF=fac*Q6fP
               !FerF=max(min(FerF,Q6fP(k)/timestep+ &
               !       min(SQ6f(k)-wsoQ6f(k),0._rk)),0._rk)
               _SET_BOTTOM_ODE_(self%id_f,-FerF)
               _SET_BOTTOM_EXCHANGE_(self%id_resuspension_f,FerF)
            end if
         end if   ! Resuspension

         ! Remineralization (benthic return)
         if (self%reminQIX/=0.0_rk) then
            if (_VARIABLE_REGISTERED_(self%id_c)) then
               _GET_HORIZONTAL_(self%id_c,Q6cP)
               _SET_BOTTOM_ODE_(self%id_c,-self%reminQIX*Q6cP)
               _SET_BOTTOM_EXCHANGE_(self%id_O3c,self%reminQIX*Q6cP/CMass)
            end if
            if (_VARIABLE_REGISTERED_(self%id_p)) then
               _GET_HORIZONTAL_(self%id_p,Q6pP)
               _SET_BOTTOM_ODE_(self%id_p,-self%reminQIX*Q6pP)
               _SET_BOTTOM_EXCHANGE_(self%id_N1p,self%reminQIX*Q6pP)
            end if
            if (_VARIABLE_REGISTERED_(self%id_n)) then
               _GET_HORIZONTAL_(self%id_n,Q6nP)
               _SET_BOTTOM_ODE_(self%id_n,-self%reminQIX*Q6nP)
               _SET_BOTTOM_EXCHANGE_(self%id_N3n,self%pQIN3X*self%reminQIX*Q6nP)
               _SET_BOTTOM_EXCHANGE_(self%id_N4n,(1.0_rk-self%pQIN3X)*self%reminQIX*Q6nP)
            end if
            if (_VARIABLE_REGISTERED_(self%id_s)) then
               _GET_HORIZONTAL_(self%id_s,Q6sP)
               _SET_BOTTOM_ODE_(self%id_s,-self%reminQIX*Q6sP)
               _SET_BOTTOM_EXCHANGE_(self%id_N5s,self%reminQIX*Q6sP)
            end if
            if (_VARIABLE_REGISTERED_(self%id_f)) then
               _GET_HORIZONTAL_(self%id_f,Q6fP)
               _SET_BOTTOM_ODE_(self%id_f,-self%reminQIX*Q6fP)
               _SET_BOTTOM_EXCHANGE_(self%id_N7f,self%reminQIX*Q6fP)
            end if
            if (_VARIABLE_REGISTERED_(self%id_l)) then
               _GET_HORIZONTAL_(self%id_l,Q6lP)
               _SET_BOTTOM_ODE_(self%id_l,-self%reminQIX*Q6lP)
               _SET_BOTTOM_EXCHANGE_(self%id_O3c,self%reminQIX*Q6lP/CMass)
            end if
         end if   ! Remineralization (benthic return)

      _HORIZONTAL_LOOP_END_

   end subroutine

end module