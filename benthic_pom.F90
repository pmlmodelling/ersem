#include "fabm_driver.h"

module pml_ersem_benthic_pom

   use fabm_types
   use pml_ersem_base

   implicit none

!  default: all is private.
   private

   type,extends(type_ersem_benthic_base_model),public :: type_pml_ersem_benthic_pom
      type (type_state_variable_id)        :: id_resuspension_c,id_resuspension_n,id_resuspension_p,id_resuspension_s,id_resuspension_f,id_resuspension_l
      type (type_state_variable_id)        :: id_O3c,id_N1p,id_N3n,id_N4n,id_N5s,id_N7f
      type (type_horizontal_dependency_id) :: id_bedstress,id_wdepth
      type (type_dependency_id)            :: id_dens
      type (type_model_id)                 :: id_resuspension_target

      real(rk) :: reminQIX,pQIN3X
      logical  :: resuspension
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_pml_ersem_benthic_pom), intent(inout), target :: self
   integer,                            intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   character(len=10) :: constituents
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(constituents,     'constituents',default='cnp')
      call self%get_parameter(self%reminQIX,    'reminQIX',    default=0.0_rk)
      call self%get_parameter(self%resuspension,'resuspension',default=.false.)

      call self%initialize_ersem_benthic_base()

      if (self%resuspension) then
         call self%register_model_dependency(self%id_resuspension_target,'RP')
         call self%register_dependency(self%id_bedstress,standard_variables%bottom_stress)
         call self%register_dependency(self%id_wdepth,   standard_variables%bottom_depth_below_geoid)
         call self%register_dependency(self%id_dens,     standard_variables%density)
      end if

      if (index(constituents,'c')/=0) then
         call self%add_constituent('c',0.0_rk)
         if (self%resuspension) then
            call self%register_state_dependency(self%id_resuspension_c,'resuspension_target_c','mmol m-3','pelagic variable taking up resuspended carbon')
            call self%request_coupling(self%id_resuspension_c,'c',source=self%id_resuspension_target)
         end if
         if (self%reminQIX/=0.0_rk) call self%register_state_dependency(self%id_O3c,'O3c','mmol m-3','dissolved inorganic carbon')
      end if
      if (index(constituents,'n')/=0) then
         call self%add_constituent('n',0.0_rk)
         if (self%resuspension) then
            call self%register_state_dependency(self%id_resuspension_n,'resuspension_target_n','mmol m-3','pelagic variable taking up resuspended nitrogen')
            call self%request_coupling(self%id_resuspension_n,'n',source=self%id_resuspension_target)
         end if
         if (self%reminQIX/=0.0_rk) then
            call self%get_parameter(self%pQIN3X,'pQIN3X',default=0.0_rk)
            call self%register_state_dependency(self%id_N3n,'N3n','mmol m-3','nitrate')
            call self%register_state_dependency(self%id_N4n,'N4n','mmol m-3','ammonium')
         end if
      end if
      if (index(constituents,'p')/=0) then
         call self%add_constituent('p',0.0_rk)
         if (self%resuspension) then
            call self%register_state_dependency(self%id_resuspension_p,'resuspension_target_p','mmol m-3','pelagic variable taking up resuspended phosphorus')
            call self%request_coupling(self%id_resuspension_p,'p',source=self%id_resuspension_target)
         end if
         if (self%reminQIX/=0.0_rk) call self%register_state_dependency(self%id_N1p,'N1p','mmol m-3','phosphate')
      end if
      if (index(constituents,'s')/=0) then
         call self%add_constituent('s',0.0_rk)
         if (self%resuspension) then
            call self%register_state_dependency(self%id_resuspension_s,'resuspension_target_s','mmol m-3','pelagic variable taking up resuspended silicate')
            call self%request_coupling(self%id_resuspension_s,'s',source=self%id_resuspension_target)
         end if
         if (self%reminQIX/=0.0_rk) call self%register_state_dependency(self%id_N5s,'N5s','mmol m-3','silicate')
      end if
      if (index(constituents,'f')/=0) then
         call self%add_constituent('f',0.0_rk)
#ifdef IRON
         if (self%resuspension) then
            call self%register_state_dependency(self%id_resuspension_f,'resuspension_target_f','umol m-3','pelagic variable taking up resuspended iron')
            call self%request_coupling(self%id_resuspension_f,'f',source=self%id_resuspension_target)
         end if
         if (self%reminQIX/=0.0_rk) call self%register_state_dependency(self%id_N7f,'N7f','umol m-3','dissolved iron')
#endif
      end if
      if (index(constituents,'l')/=0) then
         call self%add_constituent('l',0.0_rk)
         if (self%resuspension) then
            call self%register_state_dependency(self%id_resuspension_l,'resuspension_target_l','mg m-3','pelagic variable taking up resuspended calcite')
            call self%request_coupling(self%id_resuspension_l,'l',source=self%id_resuspension_target)
         end if
         if (index(constituents,'c')==0) call self%register_state_dependency(self%id_O3c,'O3c','umol m-3','dissolved inorganic carbon')
      end if
   end subroutine

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_pml_ersem_benthic_pom), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
      
      real(rk) :: ter,er,density
      real(rk) :: Q6cP,Q6nP,Q6pP,Q6sP,Q6fP,Q6lP,bedstress,wdepth
      real(rk) :: bedsedXc,bedsedXn,bedsedXp,bedsedXs
      real(rk) :: fac,FerC,FerN,FerP,FerS,FerF

      if (self%resuspension) then
         _HORIZONTAL_LOOP_BEGIN_

   !     Critical stress for erosion (m/s)
   !
         ter=0.02_rk**2
   !
   !     erosion constant (g/m^2/s)
   !
         er=100._rk*Ter * 1000._rk*86400._rk  !convert to mg/day

         _GET_HORIZONTAL_(self%id_bedstress,bedstress)
         _GET_HORIZONTAL_(self%id_wdepth,wdepth)
         _GET_HORIZONTAL_(self%id_c,Q6cP)
         _GET_HORIZONTAL_(self%id_n,Q6nP)
         _GET_HORIZONTAL_(self%id_p,Q6pP)
         _GET_(self%id_dens,density)

         ! From actual stress (Pa) to shear velocity (m/s)
         !bedstress = sqrt(bedstress/density)
         bedstress = 0.0_rk
   !
         if(bedstress.gt.ter) then
   !     Inorganic sedment (could replace by transport model)
   !     for now assume 90% is a fixed inorganic component
         bedsedXc=10._rk*(-1069._rk*LOG(wdepth) + 10900._rk)
         bedsedXn=10._rk*(-7.6368_rk*LOG(wdepth) + 78.564_rk)
         bedsedXp=10._rk*(-0.545_rk*LOG(wdepth) + 6.0114_rk)
         bedsedXs=10._rk*(-64.598_rk*LOG(wdepth) + 391.61_rk)
   !      fac = er*(bedstress(k)/ter - 1.0)
         fac = er*(bedstress/ter - 1._rk)/(Q6cP+bedsedXc)
   ! C
   !      FerC=fac*Q6cP(k)/(Q6cP(k)+bedsedXc)
         FerC=fac*Q6cP
         !FerC=max(min(FerC,Q6cP(k)/timestep+ &
         !       min(SQ6c(k)-wsoQ6c(k),0._rk)),0._rk)
   ! N
   !      FerN=fac*Q6nP(k)/(Q6nP(k)+bedsedXn)
         FerN=fac*Q6nP
         !FerN=max(min(FerN,Q6nP(k)/timestep+ &
         !       min(SQ6n(k)-wsoQ6n(k),0._rk)),0._rk)
   ! P
   !      FerP=fac*Q6pP(k)/(Q6pP(k)+bedsedXp)
         FerP=fac*Q6pP
         !FerP=max(min(FerP,Q6pP(k)/timestep+ &
         !       min(SQ6p(k)-wsoQ6p(k),0._rk)),0._rk)
         if (_VARIABLE_REGISTERED_(self%id_s)) then
   ! S
   !      FerS=fac*Q6sP(k)/(Q6sP(k)+bedsedXs)
         _GET_HORIZONTAL_(self%id_s,Q6sP)
         FerS=fac*Q6sP
         !FerS=max(min(FerS,Q6sP(k)/timestep+ &
         !       min(SQ6s(k)-wsoQ6s(k),0._rk)),0._rk)
         _SET_ODE_BEN_(self%id_s,-FerS)
         _SET_BOTTOM_EXCHANGE_(self%id_resuspension_s,FerS)
         end if
         if (_VARIABLE_REGISTERED_(self%id_f)) then
   ! F
         _GET_HORIZONTAL_(self%id_f,Q6fP)
         FerF=fac*Q6fP
         !FerF=max(min(FerF,Q6fP(k)/timestep+ &
         !       min(SQ6f(k)-wsoQ6f(k),0._rk)),0._rk)
         _SET_ODE_BEN_(self%id_f,-FerF)
         _SET_BOTTOM_EXCHANGE_(self%id_resuspension_f,FerF)
         end if

         _SET_ODE_BEN_(self%id_c,-FerC)
         _SET_ODE_BEN_(self%id_n,-FerN)
         _SET_ODE_BEN_(self%id_p,-FerP)

         _SET_BOTTOM_EXCHANGE_(self%id_resuspension_c,FerC)
         _SET_BOTTOM_EXCHANGE_(self%id_resuspension_n,FerN)
         _SET_BOTTOM_EXCHANGE_(self%id_resuspension_p,FerP)

         endif
   ! end of resuspension bit
      _HORIZONTAL_LOOP_END_
   end if
      
   ! Remineralization (benthic return)
   if (self%reminQIX/=0.0_rk) then
      _HORIZONTAL_LOOP_BEGIN_
         if (_VARIABLE_REGISTERED_(self%id_c)) then
            _GET_HORIZONTAL_(self%id_c,Q6cP)
            _SET_ODE_BEN_(self%id_c,-self%reminQIX*Q6cP)
            _SET_BOTTOM_EXCHANGE_(self%id_O3c,self%reminQIX*Q6cP/CMass)
         end if
         if (_VARIABLE_REGISTERED_(self%id_p)) then
            _GET_HORIZONTAL_(self%id_p,Q6pP)
            _SET_ODE_BEN_(self%id_p,-self%reminQIX*Q6pP)
            _SET_BOTTOM_EXCHANGE_(self%id_N1p,self%reminQIX*Q6pP)
         end if
         if (_VARIABLE_REGISTERED_(self%id_n)) then
            _GET_HORIZONTAL_(self%id_n,Q6nP)
            _SET_ODE_BEN_(self%id_n,-self%reminQIX*Q6nP)
            _SET_BOTTOM_EXCHANGE_(self%id_N3n,self%pQIN3X*self%reminQIX*Q6nP)
            _SET_BOTTOM_EXCHANGE_(self%id_N4n,(1.0_rk-self%pQIN3X)*self%reminQIX*Q6nP)
         end if
         if (_VARIABLE_REGISTERED_(self%id_s)) then
            _GET_HORIZONTAL_(self%id_s,Q6sP)
            _SET_ODE_BEN_(self%id_s,-self%reminQIX*Q6sP)
            _SET_BOTTOM_EXCHANGE_(self%id_N5s,self%reminQIX*Q6sP)
         end if
         if (_VARIABLE_REGISTERED_(self%id_f)) then
            _GET_HORIZONTAL_(self%id_f,Q6fP)
            _SET_ODE_BEN_(self%id_f,-self%reminQIX*Q6fP)
            _SET_BOTTOM_EXCHANGE_(self%id_N7f,self%reminQIX*Q6fP)
         end if
         if (_VARIABLE_REGISTERED_(self%id_l)) then
            _GET_HORIZONTAL_(self%id_l,Q6lP)
            _SET_ODE_BEN_(self%id_l,-self%reminQIX*Q6lP)
            _SET_BOTTOM_EXCHANGE_(self%id_O3c,self%reminQIX*Q6lP/CMass)
         end if
      _HORIZONTAL_LOOP_END_
   end if

   end subroutine

end module