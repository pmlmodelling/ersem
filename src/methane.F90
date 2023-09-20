#include "fabm_driver.h"

!------------------------------------------------------------------------
! This module calculates saturation concentrations and air-sea
! exchange of methane.
!------------------------------------------------------------------------

module ersem_methane

   use fabm_types
   use fabm_particle

   use ersem_shared
   use ersem_pelagic_base

   implicit none

   private

   type,extends(type_ersem_pelagic_base),public :: type_ersem_methane
      ! Variables
      type (type_dependency_id)                     :: id_ETW,id_X1X,id_chltot
      type (type_horizontal_dependency_id)          :: id_wnd,id_pCH4a
      type (type_horizontal_diagnostic_variable_id) :: id_airsea
      type (type_diagnostic_variable_id)            :: id_satp

      integer :: iswCH4

   contains
!     Model procedures
      procedure :: initialize
      procedure :: do_surface
      procedure :: do
   end type

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_ersem_methane),      intent(inout),target :: self
      integer,                         intent(in)           :: configunit
!
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%initialize_ersem_base(sedimentation=.false.)

      call self%add_constituent('c',0.0_rk)

      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      call self%register_dependency(self%id_X1X,standard_variables%practical_salinity)
      call self%register_dependency(self%id_wnd,standard_variables%wind_speed)
      call self%register_dependency(self%id_pCH4a,partial_pressure_of_ch4)
      call self%register_diagnostic_variable(self%id_satp,'satp','%','methane % saturation')
      call self%register_diagnostic_variable(self%id_airsea,'airsea','mmol CH4/m^2/d','airsea flux of CH4',source=source_do_surface)
      call self%get_parameter(self%iswCH4,'iswCH4','','air-sea flux switch (0: off, 1: on)',default=1)
   end subroutine initialize

   subroutine do(self,_ARGUMENTS_DO_)
      class (type_ersem_methane), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: CH4c,ETW,X1X,chltot
      real(rk) :: koCH4,pCH4a

       _LOOP_BEGIN_
          _GET_(self%id_c,CH4c)
          _GET_(self%id_ETW,ETW)
          _GET_(self%id_X1X,X1X)
          _GET_HORIZONTAL_(self%id_pCH4a,pCH4a)
   
         koCH4 = ch4_transfer_coefficient(self,ETW,X1X)

       _SET_DIAGNOSTIC_(self%id_satp,(100._rk*CH4c/(KoCH4*pCH4a*1.e-9/1000.)))

      _LOOP_END_
   end subroutine

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_ersem_methane), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: CH4c,ETW,X1X,wnd, fwind
      real(rk) :: koCH4,sc,ch4flux,pCH4a

      _HORIZONTAL_LOOP_BEGIN_
         _GET_(self%id_c,CH4c)
         _GET_(self%id_ETW,ETW)
         _GET_(self%id_X1X,X1X)
         _GET_HORIZONTAL_(self%id_wnd,wnd)
         _GET_HORIZONTAL_(self%id_pCH4a,pCH4a)

! Schmidt number for CH4 (Wanninkhof, 2014)
 
        sc=2101.2_rk-131.54_rk*ETW+4.4931_rk*ETW**2._rk-0.08676_rk*ETW**3._rk+0.00070663_rk*ETW**4._rk
        fwind =  0.251_rk * wnd**2 *(sc/660._rk)**(-0.5_rk)
        fwind=fwind*24._rk/100._rk  ! convert to m/day

        koCH4 = ch4_transfer_coefficient(self,ETW,X1X)

!        pCH4a is given in natm, so it is here converted to atm  by means of a 1.-9 factor.
!        To convert nmol/L to mmol/m3 we divide it by 1000.

        if (self%iswCH4 .eq. 1) then
         ch4flux = fwind*(KoCH4*pCH4a*1.e-9_rk /  1000._rk - CH4c)
        else
         ch4flux=0._rk
        endif

         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_airsea,ch4flux)

         _SET_SURFACE_EXCHANGE_(self%id_c,ch4flux)

      _HORIZONTAL_LOOP_END_
   end subroutine do_surface

   function ch4_transfer_coefficient(self,ETW,X1X) result(koCH4)
       class (type_ersem_methane), intent(in) :: self
       real(rk),                         intent(in) :: ETW,X1X
       real(rk)                                     :: koCH4
       real(rk)           :: tk,tk100

! Coefficients for temperature and salinity dependence of methane solubility 
! according to Wiesenburg and Guinasso (1979)

       real(rk),parameter :: A1 = -415.2807_rk  !-68.8862_rk
       real(rk),parameter :: A2 = 596.8104_rk   !101.4956_rk
       real(rk),parameter :: A3 = 379.2599_rk   !28.7314_rk
       real(rk),parameter :: A4 = -62.0757_rk
       real(rk),parameter :: B1 = -0.059160_rk  !-0.076146_rk
       real(rk),parameter :: B2 = 0.032174_rk   !0.043970_rk
       real(rk),parameter :: B3 = -0.0048198_rk !-0.006872_rk

       TK=ETW+273.15_rk
       TK100=TK/100._rk

        koCH4 = exp(A1+A2/tk100 + A3 * log(tk100) + A4 * tk100 + &
        &       X1X * (B1 + B2 * tk100 + B3 * tk100 ** 2._rk))
    end function
   
end module
