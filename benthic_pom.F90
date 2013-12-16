#include "fabm_driver.h"

module pml_ersem_benthic_pom

   use fabm_types
   use pml_ersem_base

   implicit none

!  default: all is private.
   private

   type,extends(type_ersem_benthic_base_model),public :: type_pml_ersem_benthic_pom
      type (type_state_variable_id)        :: id_resuspenion_c,id_resuspenion_n,id_resuspenion_p,id_resuspenion_s,id_resuspenion_f
      type (type_state_variable_id)        :: id_O3c,id_N1p,id_N3n,id_N4n,id_N5s,id_N7f
      type (type_horizontal_dependency_id) :: id_bedstress,id_wdepth
      type (type_dependency_id)            :: id_dens

      real(rk) :: reminQIX,pQIN3X
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

   real(rk),parameter :: CMass = 12._rk

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
   logical :: has_s,has_f
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(has_s,'has_s',default=.true.)
      call self%get_parameter(has_f,'has_f',default=.true.)
      call self%get_parameter(self%reminQIX,'reminQIX',default=0.0_rk)

      if (has_s.and.has_f) then
         call self%initialize_ersem_benthic_base(c_ini=0._rk,p_ini=0._rk,n_ini=0._rk,s_ini=0._rk,f_ini=0._rk)
      elseif (has_s) then
         call self%initialize_ersem_benthic_base(c_ini=0._rk,p_ini=0._rk,n_ini=0._rk,s_ini=0._rk)
      elseif (has_f) then
         call self%initialize_ersem_benthic_base(c_ini=0._rk,p_ini=0._rk,n_ini=0._rk,f_ini=0._rk)
      else
         call self%initialize_ersem_benthic_base(c_ini=0._rk,p_ini=0._rk,n_ini=0._rk)
      end if

      call self%register_state_dependency(self%id_resuspenion_c,'resuspension_target_c','mmol m-3','pelagic variable taking up resuspended carbon',    required=.false.)
      call self%register_state_dependency(self%id_resuspenion_n,'resuspension_target_n','mmol m-3','pelagic variable taking up resuspended nitrogen',  required=.false.)
      call self%register_state_dependency(self%id_resuspenion_p,'resuspension_target_p','mmol m-3','pelagic variable taking up resuspended phosphorus',required=.false.)
      if (has_s) call self%register_state_dependency(self%id_resuspenion_s,'resuspension_target_s','mmol m-3','pelagic variable taking up resuspended silicate',required=.false.)
      if (has_f) call self%register_state_dependency(self%id_resuspenion_f,'resuspension_target_f','umol m-3','pelagic variable taking up resuspended iron',required=.false.)

      if (self%reminQIX/=0.0_rk) then
         call self%get_parameter(self%pQIN3X,'pQIN3X',default=0.0_rk)
         call self%register_state_dependency(self%id_O3c,'O3c','mmol m-3','dissolved inorganic carbon')
         call self%register_state_dependency(self%id_N1p,'N1p','mmol m-3','phosphate')
         call self%register_state_dependency(self%id_N3n,'N3n','mmol m-3','nitrate')
         call self%register_state_dependency(self%id_N4n,'N4n','mmol m-3','ammonium')
         if (has_s) call self%register_state_dependency(self%id_N5s,'N5s','mmol m-3','silicate')
         if (has_f) call self%register_state_dependency(self%id_N7f,'N7f','umol m-3','dissolved iron')
      end if

      call self%register_dependency(self%id_bedstress,standard_variables%bottom_stress)
      call self%register_dependency(self%id_wdepth,   standard_variables%bottom_depth_below_geoid)
      call self%register_dependency(self%id_dens,     standard_variables%density)
   end subroutine

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_pml_ersem_benthic_pom), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
      
      real(rk) :: ter,er,density
      real(rk) :: Q6cP,Q6nP,Q6pP,Q6sP,Q6fP,bedstress,wdepth
      real(rk) :: bedsedXc,bedsedXn,bedsedXp,bedsedXs
      real(rk) :: fac,FerC,FerN,FerP,FerS,FerF

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
      _SET_BOTTOM_EXCHANGE_(self%id_resuspenion_s,FerS)
      end if
      if (_VARIABLE_REGISTERED_(self%id_f)) then
! F
      _GET_HORIZONTAL_(self%id_f,Q6fP)
      FerF=fac*Q6fP
      !FerF=max(min(FerF,Q6fP(k)/timestep+ &
      !       min(SQ6f(k)-wsoQ6f(k),0._rk)),0._rk)
      _SET_ODE_BEN_(self%id_f,-FerF)
      _SET_BOTTOM_EXCHANGE_(self%id_resuspenion_f,FerF)
      end if

      _SET_ODE_BEN_(self%id_c,-FerC)
      _SET_ODE_BEN_(self%id_n,-FerN)
      _SET_ODE_BEN_(self%id_p,-FerP)

      _SET_BOTTOM_EXCHANGE_(self%id_resuspenion_c,FerC)
      _SET_BOTTOM_EXCHANGE_(self%id_resuspenion_n,FerN)
      _SET_BOTTOM_EXCHANGE_(self%id_resuspenion_p,FerP)

      endif
! end of resuspension bit

      if (self%reminQIX/=0.0_rk) then
         ! Remineralization (benthic return)
         _SET_ODE_BEN_(self%id_c,-self%reminQIX*Q6cP)
         _SET_ODE_BEN_(self%id_p,-self%reminQIX*Q6pP)
         _SET_ODE_BEN_(self%id_n,-self%reminQIX*Q6nP)
         _SET_BOTTOM_EXCHANGE_(self%id_O3c,self%reminQIX*Q6cP/CMass)
         _SET_BOTTOM_EXCHANGE_(self%id_N1p,self%reminQIX*Q6pP)
         _SET_BOTTOM_EXCHANGE_(self%id_N3n,self%pQIN3X*self%reminQIX*Q6nP)
         _SET_BOTTOM_EXCHANGE_(self%id_N4n,(1.0_rk-self%pQIN3X)*self%reminQIX*Q6nP)
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
      end if
   _HORIZONTAL_LOOP_END_

   end subroutine

end module