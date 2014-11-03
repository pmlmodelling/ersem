#include "fabm_driver.h"

module ersem_benthic_column

   use fabm_types

   use ersem_shared

   implicit none

   private

   ! Model that specifies the structure of the three-layer sediment column
   type,extends(type_base_model),public :: type_ersem_benthic_column
      type (type_bottom_state_variable_id)          :: id_D1m,id_D2m
      type (type_horizontal_diagnostic_variable_id) :: id_poro,id_Dtot,id_diff(3),id_EDZ_mixX

      real(rk) :: EDZ_1X,EDZ_2X,EDZ_3X,d_totX
      real(rk) :: qPWX,EDZ_mixX
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

   ! Model for nitrification
   type,extends(type_base_model),public :: type_ersem_benthic_nitrification
      type (type_bottom_state_variable_id) :: id_K3n,id_K4n,id_G2o
      type (type_state_variable_id)        :: id_N4n
      type (type_dependency_id)            :: id_ETW,id_phx
      type (type_horizontal_dependency_id) :: id_D1m

      real(rk) :: q10nitX,hM4M3X,sM4M3X,xno3X
      integer :: ISWphx
   contains
      procedure :: initialize => benthic_nitrification_initialize
      procedure :: do_bottom  => benthic_nitrification_do_bottom
   end type

contains

   subroutine initialize(self,configunit)
      class (type_ersem_benthic_column),intent(inout),target :: self
      integer,                                     intent(in)           :: configunit
      
      ! Set time unit to d-1. This implies that all rates (sink/source terms) are given in d-1.
      self%dt = 86400._rk

      call self%get_parameter(self%EDZ_1X,'EDZ_1X','m^2/d','diffusivity in 1st (oxygenated) layer')
      call self%get_parameter(self%EDZ_2X,'EDZ_2X','m^2/d','diffusivity in 2nd (oxidized) layer')
      call self%get_parameter(self%EDZ_3X,'EDZ_3X','m^2/d','diffusivity in 3rd (anoxic) layer')
      call self%get_parameter(self%qPWX,'qPWX','-','fraction of pore water in the sediment')
      call self%get_parameter(self%EDZ_mixX,'EDZ_mixX','d/m','equilibrium diffusive speed between sediment surface water')
      call self%get_parameter(self%d_totX,'d_totX','m','depth of sediment column')

      call self%register_state_variable(self%id_D1m,'D1m','m','depth of bottom interface of 1st layer',standard_variable=depth_of_bottom_interface_of_layer_1)
      call self%register_state_variable(self%id_D2m,'D2m','m','depth of bottom interface of 2nd layer',standard_variable=depth_of_bottom_interface_of_layer_2)
      call self%register_diagnostic_variable(self%id_poro,'poro','-','porosity',standard_variable=sediment_porosity,missing_value=self%qPWX)
      call self%register_diagnostic_variable(self%id_Dtot,'Dtot','m','depth of sediment column',missing_value=self%d_totX,standard_variable=depth_of_sediment_column)
      call self%register_diagnostic_variable(self%id_diff(1),'diff1','m^2/s','diffusivity in layer 1',standard_variable=diffusivity_in_sediment_layer_1,missing_value=self%EDZ_1X)
      call self%register_diagnostic_variable(self%id_diff(2),'diff2','m^2/s','diffusivity in layer 2',standard_variable=diffusivity_in_sediment_layer_2,missing_value=self%EDZ_2X)
      call self%register_diagnostic_variable(self%id_diff(3),'diff3','m^2/s','diffusivity in layer 3',standard_variable=diffusivity_in_sediment_layer_3,missing_value=self%EDZ_3X)
      call self%register_diagnostic_variable(self%id_EDZ_mixX,'cmix','s/m','equilibrium diffusive speed between sediment surface water',standard_variable=pelagic_benthic_transfer_constant,missing_value=self%EDZ_mixX)
   end subroutine

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_benthic_column),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
   
      _HORIZONTAL_LOOP_BEGIN_

      _HORIZONTAL_LOOP_END_
   end subroutine do_bottom

   subroutine benthic_nitrification_initialize(self,configunit)
      class (type_ersem_benthic_nitrification),intent(inout),target :: self
      integer,                                     intent(in)           :: configunit

      ! Set time unit to d-1. This implies that all rates (sink/source terms) are given in d-1.
      self%dt = 86400._rk

      call self%get_parameter(self%q10nitX,'q10nitX','-',            'Q_10 temperature coefficient mfor nitrification')
      call self%get_parameter(self%hM4M3X, 'hM4M3X', 'mmol/m^3',     'Michaelis-Menten constant for nitrate limitation of nitrification')
      call self%get_parameter(self%ISWphx, 'ISWphx', '',             'pH influence on nitrification',default=0)
      call self%get_parameter(self%sM4M3X, 'sM4M3X', '1/d',          'maximum nitrification rate')
      call self%get_parameter(self%xno3X,  'xno3X',  'mol O_2/mol N','oxygen consumed per nitrate produced in nitrification')

      call self%register_state_dependency(self%id_K3n,'K3n','mmol/m^2','nitrate')
      call self%register_state_dependency(self%id_K4n,'K4n','mmol/m^2','ammonium')
      call self%register_state_dependency(self%id_G2o,'G2o','mmol/m^2','oxygen')
      call self%register_state_dependency(self%id_N4n,'N4n','mmol N/m^3','pelagic ammonium')
      call self%register_dependency(self%id_D1m,'D1m','m','depth of bottom interface of 1st layer',standard_variable=depth_of_bottom_interface_of_layer_1)

      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      if (self%ISWphx==1) call self%register_dependency(self%id_phx,standard_variables%ph_reported_on_total_scale)
   end subroutine benthic_nitrification_initialize

   subroutine benthic_nitrification_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_benthic_nitrification),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: K3n,K3nP,K4nP,N4n
      real(rk) :: ETW,phx
      real(rk) :: Mu_m,eT,eN,Fph,jM4M3n,D1m

      _HORIZONTAL_LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_K3n,K3n) !TODO background
         _GET_HORIZONTAL_(self%id_K3n,K3nP)
         _GET_HORIZONTAL_(self%id_K4n,K4nP)
         _GET_(self%id_N4n,N4n)

         _GET_HORIZONTAL_(self%id_D1m,D1m)
      
         _GET_(self%id_ETW,ETW)

         !* Nitrification (layer 1 only):
         ! Reaction: (NH4+,OH-) + 2*O2 -> (NO3-,H+) + 2*H2O
         ! treated as first order reaction for NH4 in the oxygenated layer
         ! Average nitrate density (mmol/m3): MU_m

         ! Average nitrate density in first/oxygenated layer
         MU_m = max(0.0_rk,K3n)/D1m

         ! Limitation factor for temperature (Q_10)
         eT = self%q10nitX**((ETW-10._rk)/10._rk)

         ! Limitation factor for nitrate (hyperbolic, 1 at zero nitrate, dropping to 0 at high nitrate).
         eN = self%hM4M3X/(self%hM4M3X+MU_m)

         ! Ph influence on nitrification - empirical equation
         if(self%ISWphx==1) then
            _GET_(self%id_phx,phx)
           Fph = min(2._rk,max(0._rk,0.6111_rk*phx-3.8889_rk))
           jM4M3n = Fph * self%sM4M3X * K4nP * eT * eN
         else
           jM4M3n = self%sM4M3X*K4nP*eT*eN
         end if

         ! Impose nitrification source-sink terms for benthic ammonium, nitrate, oxygen.
         _SET_BOTTOM_ODE_(self%id_K3n,jM4M3n)
         _SET_BOTTOM_ODE_(self%id_K4n,-jM4M3n)
         _SET_BOTTOM_ODE_(self%id_G2o,-self%xno3X*jM4M3n)
      _HORIZONTAL_LOOP_END_

   end subroutine benthic_nitrification_do_bottom

end module
