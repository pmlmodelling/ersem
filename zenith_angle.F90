#include "fabm_driver.h"

module ersem_zenith_angle

   use fabm_types
   use ersem_shared
   use ersem_pelagic_base

   implicit none

   private

   type,extends(type_ersem_pelagic_base),public :: type_ersem_zenith_angle
      ! Variables
      type (type_horizontal_dependency_id)     :: id_lon,id_lat
      type (type_horizontal_diagnostic_variable_id)     :: id_zenith_angle
      type (type_global_dependency_id)     :: id_yday
   contains
!     Model procedures
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_ersem_zenith_angle),intent(inout),target :: self
      integer,                         intent(in)           :: configunit
!
!EOP
!-----------------------------------------------------------------------
!BOC
      ! Initialize pelagic base model (this also sets the time unit to per day, instead of the default per second)
      call self%initialize_ersem_base(sedimentation=.false.)

      ! Register environmental dependencies (lon,lat,zenith,day of year)
      call self%register_horizontal_dependency(self%id_lon,standard_variables%longitude)
      call self%register_horizontal_dependency(self%id_lat,standard_variables%latitude)
      call self%register_horizontal_diagnostic_variable(self%id_zenith_angle,'zenithA','degrees','zenith_angle',standard_variable=zenith_angle)
      call self%register_global_dependency(self%id_yday,standard_variables%number_of_days_since_start_of_the_year)

   end subroutine initialize
   
   subroutine do(self,_ARGUMENTS_DO_)

      class (type_ersem_zenith_angle),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: lon,lat,h,th0,th02,th03,coszen,sundec,thsun,yday,zenithA
      integer :: iday

      ! Leave spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve longitue, latitude and time since beginning of the year
         _GET_HORIZONTAL_(self%id_lon,lon)
         _GET_HORIZONTAL_(self%id_lat,lat)
         _GET_GLOBAL_(self%id_yday,yday)

         !Get day integer and day fraction in hours:
         iday=int(yday)+1
         h=(yday+1-iday)*24._rk

         ! Leap year not considered:
         th0 = pi*iday/182.5_rk
         th02 = 2._rk*th0
         th03 = 3._rk*th0

         ! Sun declination :
         sundec = 0.006918_rk - 0.399912_rk*cos(th0) + 0.070257_rk*sin(th0)         &
               - 0.006758_rk*cos(th02) + 0.000907_rk*sin(th02)                 &
               - 0.002697_rk*cos(th03) + 0.001480_rk*sin(th03)

         ! Sun hour angle :
         zenithA = (h-12._rk)*15._rk*deg2rad + lon

         ! Cosine of the solar zenith angle :
         zenithA =sin(lat)*sin(sundec)+cos(lat)*cos(sundec)*cos(zenithA)
         zenithA =acos(max(0._rk,zenithA))/deg2rad
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_zenith_angle,zenithA)
      ! Leave spatial loops (if any)
      _LOOP_END_

   end subroutine do
   
end module
