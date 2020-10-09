#include "fabm_driver.h"

! Benthic variable that supports resuspension and remineralization.
! Both processes return material to the pelagic.

module ersem_benthic_calcite

   use fabm_types

   use fabm_particle
   use ersem_shared
   use ersem_benthic_base

   implicit none

!  default: all is private.
   private

   type,extends(type_ersem_benthic_base),public :: type_ersem_benthic_calcite
      type (type_horizontal_diagnostic_variable_id) :: id_dissolution
      type (type_dependency_id)                     :: id_Om_Cal

      ! Parameters
      real(rk) :: fdissmax, fdissmin, ndiss, KcalomX
      integer  :: iswcal
      
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
   class (type_ersem_benthic_calcite), intent(inout), target :: self
   integer,                         intent(in)            :: configunit

   real(rk) :: c0
!
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%initialize_ersem_benthic_base()

      call self%get_parameter(self%iswcal,'iswcal','','dissolution dependence on calcite saturation (0: none, 1: power law, 2: hyperbolic)', minimum=0, maximum=2)
      select case (self%iswcal)
         case (1)
            call self%get_parameter(self%ndiss,'ndiss','-','power of the dissolution law (Keir 1980)')
         case (2)
            call self%get_parameter(self%KcalomX,'KcalomX','-','half-saturation constant for calcification limitation from saturation state')
      end select
      if (self%iswcal == 0) then
         call self%get_parameter(self%fdissmin, 'fdiss', '1/d','specific dissolution rate', default=0.0_rk)
         self%fdissmax = 0.0_rk
      else
         call self%get_parameter(self%fdissmax, 'fdissmax', '1/d','maximum specific dissolution rate', minimum=0._rk, default=0.0_rk)
         call self%get_parameter(self%fdissmin, 'fdissmin', '1/d','minimum specific dissolution rate', minimum=0._rk, default=0.001_rk * self%fdissmax)
         call self%register_dependency(self%id_Om_Cal,'Om_Cal','-','calcite saturation')
      end if
      call self%get_parameter(c0,'c0','mg C/m^2','background calcite concentration',default=0.0_rk)

      call self%add_constituent('c',0.0_rk,c0)
      call self%register_state_dependency(self%id_O3c,'O3c','mmol/m^3','dissolved inorganic carbon')
      call self%register_state_dependency(self%id_TA,standard_variables%alkalinity_expressed_as_mole_equivalent)
      call self%register_diagnostic_variable(self%id_dissolution,'dissolution','mg C/m^2/d','dissolution',source=source_do_bottom)

   end subroutine

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)

      class (type_ersem_benthic_calcite),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: bL2c
      real(rk) :: Om_Cal
      real(rk) :: fdiss

      _HORIZONTAL_LOOP_BEGIN_

         _GET_HORIZONTAL_(self%id_c, bL2c)
         _GET_(self%id_om_cal, om_cal)

         if (self%iswcal==0) then  ! NB select case would be cleaner but makes vectorization impossible for ifort 14
            fdiss = 0._rk
         elseif (self%iswcal==1) then
            fdiss = (max(1._rk-om_cal,0._rk))**self%ndiss
         else
            fdiss = max(0._rk,(1._rk-om_cal)/(1._rk-om_cal+self%KcalomX))
         end if

         fdiss = max(fdiss * self%fdissmax, self%fdissmin)

         _SET_BOTTOM_ODE_(self%id_c, -fdiss*bL2c)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dissolution, -fdiss*bL2c)
         _SET_BOTTOM_EXCHANGE_(self%id_O3c, fdiss*bL2c/CMass)
         _SET_BOTTOM_EXCHANGE_(self%id_TA, 2*fdiss*bL2c/CMass)  ! Dissolution of CaCO3 increases alkalinity by 2 units

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module
