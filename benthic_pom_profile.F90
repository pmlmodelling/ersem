#include "fabm_driver.h"

! Benthic particulate organic matter variable with idealized profile.
!
! Its profile is defined by its density (quantity/m2) and its penetration depth.
! By assuming an exponential distribution of matter (constant positive concentration
! at sediment surface, tending to zero at infinite depth, these two variables suffice
! to specify the concentration profile.
!
! From the idealized concentration profile, this module computes average concentrations
! per layer. Other modules can retirev these values, but also provide their rate of change.
! This rate of change is then converted into a change in total depth-integrated matter
! and a chjnage in pnentration depth.

module pml_ersem_benthic_pom_profile

   use fabm_types
   use fabm_particle

   use pml_ersem_shared
   use pml_ersem_benthic_base

   implicit none

!  default: all is private.
   private

   type,extends(type_ersem_benthic_base_model),public :: type_ersem_benthic_pom_profile_model
      type (type_bottom_state_variable_id) :: id_pen_depth
      type (type_horizontal_dependency_id) :: id_D1m, id_D2m
      type (type_horizontal_dependency_id) :: id_c_l1_sms
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

   type,extends(type_base_model) :: type_layer_content_calculator
      ! Dependencies for layer-specific sink and source terms
      type (type_horizontal_dependency_id) :: id_c_tot,id_pen_depth
      type (type_horizontal_diagnostic_variable_id) :: id_c
   contains
      procedure :: do_bottom  => layer_content_calculator_do_bottom
   end type

contains

   subroutine initialize(self,configunit)
   class (type_ersem_benthic_pom_profile_model), intent(inout), target :: self
   integer,                                      intent(in)            :: configunit

   character(len=10) :: composition
   class (type_layer_content_calculator),pointer :: translator
!EOP
!-----------------------------------------------------------------------
!BOC
   call self%get_parameter(composition,'composition','','elemental composition',default='cnp')
   if (index(composition,'c')/=0) then
      call self%add_constituent('c',0.0_rk)
      call self%register_state_variable(self%id_pen_depth,'pen_depth','cm','penetration depth')
      call self%register_dependency(self%id_c_l1_sms,'c_l1_sms','mg C/m^2/s','carbon in layer 1 sinks-sources')
      call self%request_coupling(self%id_c_l1_sms,'l1/c_sms')
      call self%register_dependency(self%id_D1m,'D1m','m','depth of bottom interface of 1st layer')
      call self%register_dependency(self%id_D2m,'D2m','m','depth of bottom interface of 2nd layer')

      ! Create a module that will compute the contents of specific layers
      allocate(translator)
      call self%add_child(translator,'l1',configunit=configunit)
      call translator%register_dependency(translator%id_c_tot,    'c_tot',    'mg C/m^2/s','depth-integrated carbon')
      call translator%register_dependency(translator%id_pen_depth,'pen_depth','cm',        'penetration depth')
      call translator%register_diagnostic_variable(translator%id_c,'c','mg C/m^3','carbon')
      call translator%act_as_state_variable(translator%id_c)
      call self%request_coupling(translator%id_pen_depth,'pen_depth')
   end if
   end subroutine initialize

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_benthic_pom_profile_model), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: z_pen,D1m,D2m
      real(rk) :: density
      real(rk) :: sms1
      
      _HORIZONTAL_LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_c,        density)
         _GET_HORIZONTAL_(self%id_pen_depth,z_pen)
         _GET_HORIZONTAL_(self%id_D1m,      D1m)
         _GET_HORIZONTAL_(self%id_D2m,      D2m)
         
         ! Change in penetration depth can be derived by considering that the current penetration depth (z_pen)
         ! and mass density (c_int) are perturbed by addition of mass (delta_c) at some known depth (z_sms)
         ! New penetration depth = (z_pen*c_int + z_sms*delta_c)/(c_int+delta_c) = z_pen + delta_z
         ! delta_z = (z_pen*c_int + z_sms*delta_c)/(c_int+delta_c) - z_pen
         !         = (z_sms-z_pen)*delta_c/(c_int+delta_c)
         ! As we are considering change over an infinitessimal time, delta_c<<c_int, and we have obtain
         ! delta_z = (z_sms-z_pen)*delta_c/c_int

         ! Process change in layer 1
         _GET_HORIZONTAL_(self%id_c_l1_sms, sms1)
         _SET_BOTTOM_ODE_(self%id_c, sms1)
         _SET_BOTTOM_ODE_(self%id_pen_depth, (0.5_rk*D1m-z_pen)*sms1/density)
      _HORIZONTAL_LOOP_END_
   end subroutine do_bottom

   subroutine layer_content_calculator_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_layer_content_calculator), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: pen_depth
      real(rk) :: density
      real(rk) :: d_top, d_bot, d_totX

      _HORIZONTAL_LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_c_tot,    density)
         _GET_HORIZONTAL_(self%id_pen_depth,pen_depth)

         ! TODO: determine d_top, d_bot, d_totX

         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_c, density*partQ(pen_depth, d_top, d_bot, d_totX))

      _HORIZONTAL_LOOP_END_
   end subroutine layer_content_calculator_do_bottom

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: partQ \label{sec:partQ}
!
! !DESCRIPTION:
!  TODO - CHECK THIS
!
!  Calculates the fraction of detritus between d\_top and d\_bot
!
!  IN:  Penetration depth of detrital component...........d\_pen
!       Top of detrital layer.............................d\_top
!       Bottom of detrital layer..........................d\_bot
!       Maximum depth of detrital layer...................d\_max
!
!  OUT: Fraction of detritus between d\_top and t\_bot......partQ
!
!  NOTE: Does not treat silicate
!
!  NOTE: Factor 13.8
!        This mystic factor takes care that exp(-..) won't become too
!        small so that the simulation crashes.  -13.8 is the smallest
!        number so that exp(-13.8) is evaluated correctly. It depends
!        on the accuracy of the implemented FORTRAN. (Cora)
!
!        We had with our FORTRAN some problems,
!        it was unable to calculate $exp(-x)$ with $x>13.8$.
!        This should be 0 ($<10^(-30)$ or something). But on runtime
!        we got an algebraic underflow. So we checked x before
!        going to the exponential function. If you leave out
!        this comparison, check first, how your FORTRAN treats
!        exp(-15),exp(-50). (Wolfgang)
!\\
!\\
! !INTERFACE:
   real(rk) function partQ( d_pen, d_top, d_bot, d_max )
!
! !INPUT PARAMETERS:
!     ! TODO - Document these
      real(rk), intent(in) :: d_pen, d_top, d_bot, d_max
!
! !LOCAL VARIABLES:
!     ! TODO - document these
      real(rk) :: norm, d_top1, d_bot1, d_max1
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!

      d_max1 = d_max  ! removed outdated protection against underflow: MIN( d_pen*13.8_rk, d_max )
      d_bot1 = min(d_bot, d_max1)
      d_top1 = min(d_top, d_bot1)

      if ( d_max1>0._rk ) then
         ! Penetration depth > 0: integrate idealized [exponential] distribution over desired depth interval.
         
         ! Compute normalization factor: integral of exponential distribution from surface to bottom of column.
         norm = 1._rk - exp(-d_max1/d_pen)

         ! Compute integral of exponential over desired depth interval and normalize to obtain fraction between 0 and 1.
         partQ = (exp(-d_top1/d_pen) - exp(-d_bot1/d_pen)) / norm
      else
         ! Penetration depth = 0 (or < 0, but that's an artefact): all mass in surface layer of zero thickness.
         if (d_top==0._rk) then
            ! Desired depth range includes surface, so include all.
            partQ = 1._rk
         else
            ! Desired depth range excludes surface, so exclude all.
            partQ = 0._rk
         end if
      end if

   end function partQ
!
!EOC
!-----------------------------------------------------------------------

end module