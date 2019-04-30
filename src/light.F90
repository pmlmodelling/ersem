#include "fabm_driver.h"

module ersem_light

   use fabm_types

   implicit none

   private

   type,extends(type_base_model),public :: type_ersem_light
      ! Identifiers for diagnostic variables
      type (type_diagnostic_variable_id)   :: id_EIR, id_parEIR, id_xEPS
      type (type_dependency_id)            :: id_dz, id_xEPSp, id_ESS
      type (type_horizontal_dependency_id) :: id_I_0

      ! Parameters
      real(rk) :: EPSESSX,EPS0X,pEIR_eowX
   contains
!     Model procedures
      procedure :: initialize
      procedure :: get_light
   end type type_ersem_light

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_ersem_light),intent(inout),target :: self
      integer,                     intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%EPS0X,    'EPS0r',   '1/m',   'background shortwave attenuation', default=4.E-2_rk)
      call self%get_parameter(self%EPSESSX,  'EPSESS',  'm^2/mg','specific shortwave attenuation of silt', default=4.E-5_rk)
      call self%get_parameter(self%pEIR_eowX,'pEIR_eow','-',     'photosynthetically active fraction of shortwave radiation', default=.5_rk)

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_EIR,'EIR','W/m^2','shortwave radiation', &
              standard_variable=standard_variables%downwelling_shortwave_flux,source=source_do_column)
      call self%register_diagnostic_variable(self%id_parEIR,'parEIR','W/m^2','photosynthetically active radiation', &
              standard_variable=standard_variables%downwelling_photosynthetic_radiative_flux,source=source_do_column)
      call self%register_diagnostic_variable(self%id_xEPS,'xEPS','1/m','attenuation coefficient of shortwave flux', &
              source=source_do_column)

      ! Register environmental dependencies (temperature, shortwave radiation)
      call self%register_dependency(self%id_I_0,standard_variables%surface_downwelling_shortwave_flux)
      call self%register_dependency(self%id_dz, standard_variables%cell_thickness)
      call self%register_dependency(self%id_xEPSp,standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_ESS, type_bulk_standard_variable(name='mass_concentration_of_silt'))
   end subroutine

   subroutine get_light(self,_ARGUMENTS_VERTICAL_)
      class (type_ersem_light),intent(in) :: self
      _DECLARE_ARGUMENTS_VERTICAL_

      real(rk) :: buffer,dz,xEPS,xtnc,EIR,ESS

      _GET_HORIZONTAL_(self%id_I_0,buffer)

      if (buffer.lt.0._rk) buffer=0._rk

      _VERTICAL_LOOP_BEGIN_
         _GET_(self%id_dz,dz)     ! Layer height (m)
         _GET_(self%id_xEPSp,xEPS) ! Extinction coefficient of shortwave radiation, due to particulate organic material (m-1)
         _GET_(self%id_ESS,ESS)   ! Suspended silt
         xEPS = xEPS + self%EPS0X + self%EPSESSX*ESS
         xtnc = xEPS*dz
         EIR = buffer/xtnc*(1.0_rk-exp(-xtnc))  ! Note: this computes the vertical average, not the value at the layer centre.
         buffer = buffer*exp(-xtnc)
         _SET_DIAGNOSTIC_(self%id_EIR,EIR)                     ! Local shortwave radiation
         _SET_DIAGNOSTIC_(self%id_parEIR,EIR*self%pEIR_eowX)   ! Local photosynthetically active radiation
         _SET_DIAGNOSTIC_(self%id_xEPS,xEPS)                   ! Vertical attenuation of shortwave radiation
      _VERTICAL_LOOP_END_

   end subroutine get_light

end module
