#include "fabm_driver.h"

! This module calculates the solar zenith angle as a function of longitude, latitude, day of year, and time of day.
! The code that performs this calculation has been adapted from the General Ocean Turbulence Model (GOTM, http://www.gotm.net).
! Original GOTM source file: src/airsea/solar_zenith_angle.F90.

module ersem_zenith_angle

   use fabm_types
   use ersem_shared

   implicit none

   private

   type,extends(type_base_model),public :: type_ersem_zenith_angle
      ! Variables
      type (type_horizontal_dependency_id)          :: id_lon,id_lat
      type (type_horizontal_diagnostic_variable_id) :: id_zenith_angle
      type (type_global_dependency_id)              :: id_yday
   contains
!     Model procedures
      procedure :: initialize
      procedure :: do_surface
   end type

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_ersem_zenith_angle),intent(inout),target :: self
      integer,                        intent(in)           :: configunit
!
!EOP
!-----------------------------------------------------------------------
!BOC
      ! Register environmental dependencies (lon,lat,zenith,day of year)
      call self%register_horizontal_dependency(self%id_lon,standard_variables%longitude)
      call self%register_horizontal_dependency(self%id_lat,standard_variables%latitude)
      call self%register_horizontal_diagnostic_variable(self%id_zenith_angle,'zenithA','degrees','zenith_angle',standard_variable=zenith_angle,source=source_do_surface,missing_value=0._rk)
      call self%register_global_dependency(self%id_yday,standard_variables%number_of_days_since_start_of_the_year)

   end subroutine initialize

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_ersem_zenith_angle), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: lon,lat,h,th0,th02,th03,sundec,yday,zenithA
      integer :: iday

      ! Retrieve time since beginning of the year
      _GET_GLOBAL_(self%id_yday,yday)

      !Get day integer and day fraction in hours:
      iday=int(yday)+1
      h=(yday+1-iday)*24._rk

      ! Leap year not considered:
      th0 = pi*iday/182.5_rk
      th02 = 2._rk*th0
      th03 = 3._rk*th0

      ! Sun declination :
      sundec = 0.006918_rk - 0.399912_rk*cos(th0) + 0.070257_rk*sin(th0) &
            - 0.006758_rk*cos(th02) + 0.000907_rk*sin(th02)              &
            - 0.002697_rk*cos(th03) + 0.001480_rk*sin(th03)

      ! Leave spatial loops (if any)
      _HORIZONTAL_LOOP_BEGIN_

         ! Retrieve longitude, latitude
         _GET_HORIZONTAL_(self%id_lon,lon)
         _GET_HORIZONTAL_(self%id_lat,lat)

         ! Sun hour angle :
         zenithA = ((h-12._rk)*15._rk + lon)*deg2rad

         ! Cosine of the solar zenith angle :
         zenithA =sin(lat*deg2rad)*sin(sundec)+cos(lat*deg2rad)*cos(sundec)*cos(zenithA)
         zenithA =acos(max(0._rk,zenithA))/deg2rad
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_zenith_angle,zenithA)

      ! Leave spatial loops (if any)
      _HORIZONTAL_LOOP_END_

   end subroutine do_surface

end module
