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

      ! Coupled state variables for resuspension and remineralization

      type (type_horizontal_diagnostic_variable_id) :: id_ben_dissolution
      type (type_dependency_id)     :: id_Om_Cal

      ! Parameters
      real(rk) :: sL2O3X,ndiss,KcalomX
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

   character(len=10) :: composition
   real(rk) :: c0
!
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%initialize_ersem_benthic_base()

      call self%get_parameter(self%sL2O3X, 'sL2O3X', '1/d','dissolution rate at 10 degrees Celsius',default=0.0_rk)
      call self%get_parameter(c0,'c0','mg C/m^2','background carbon concentration',default=0.0_rk)
      call self%get_parameter(self%iswcal,'iswcal','','calcification/dissolution dependence on calcite saturation (1: power law, 2: hyperbolic)')
      
      if (self%iswcal<0.or.self%iswcal>2) then
         call self%log_message('ISWcal set to 1')
         self%iswcal = 1
      end if
      select case (self%iswcal)
         case (1)
            call self%get_parameter(self%ndiss,'ndiss','-','power of the dissolution law (Keir 1980)')
         case (2)
            call self%get_parameter(self%KcalomX,'KcalomX','-','half-saturation constant for calcification limitation from saturation state')
      end select
      
      call self%add_constituent('c',0.0_rk,c0)
      call self%register_dependency(self%id_Om_Cal,'Om_Cal','-','calcite saturation')
      call self%register_state_dependency(self%id_O3c,'O3c','mmol/m^3','dissolved inorganic carbon')
      call self%register_state_dependency(self%id_TA,standard_variables%alkalinity_expressed_as_mole_equivalent)
      call self%register_diagnostic_variable(self%id_ben_dissolution,'ben_dissolution','mg C/m^2/d','calcite benthic dissolution',source=source_do_bottom)
      
   end subroutine

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)

      class (type_ersem_benthic_calcite),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: bL2c
      real(rk) :: Om_Cal, O3c, TA
      real(rk) :: fdiss

      _HORIZONTAL_LOOP_BEGIN_

         if (legacy_ersem_compatibility) then
            ! Legacy ERSEM includes background value, but this is inappropriate as it is used in a sink term.
            _GET_HORIZONTAL_WITH_BACKGROUND_(self%id_c,bL2c)
         else
            _GET_HORIZONTAL_(self%id_c,bL2c)
         end if
         _GET_(self%id_om_cal,om_cal)

         if (self%iswcal==0) then  ! NB select case would be cleaner but makes vectorization impossible for ifort 14
            fdiss = 0._rk
         elseif (self%iswcal==1) then
            fdiss = (max(1._rk-om_cal,0._rk))**self%ndiss
         else
            fdiss = max(0._rk,(1._rk-om_cal)/(1._rk-om_cal+self%KcalomX))
         end if

         fdiss = max(fdiss,0.001_rk )* self%sL2O3X

         _SET_BOTTOM_ODE_(self%id_c,  -fdiss*bL2c)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ben_dissolution,-fdiss*bL2c)
         _SET_BOTTOM_EXCHANGE_(self%id_O3c, fdiss*bL2c/CMass)
         _SET_BOTTOM_EXCHANGE_(self%id_TA,2*fdiss*bL2c/CMass)  ! Dissolution of CaCO3 increases alkalinity by 2 units
      
      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module
