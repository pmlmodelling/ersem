#include "fabm_driver.h"

module ersem_benthic_column_particulate_matter

! Particulate organic matter with idealized [exponential] profile.
!
! This file contains two modules that can be instantiated by the user (from fabm.yaml):
! type_ersem_benthic_column_particulate_matter: describes particulate organic matter in terms of column-integrated mass and penetration depth
! type_ersem_benthic_pom_layer:        describes particulate organic matter within a user-specified depth interval
!
! The idealized profile is defined only in terms of its density (quantity/m^2) and its penetration depth.
! By assuming an exponential distribution of matter (constant positive concentration
! at sediment surface, tending to zero at infinite depth), these two variables suffice
! to specify the concentration profile.
!
! From the idealized concentration profile, we first compute densities (quantity/m^2) per layer.
! These layer-specific densities act as state variable: other modules can retrieve their value,
! but also provide their rate of change. We finally retrieve these rates of change, and convert them
! into a change in total depth-integrated matter *and* a change in penetration depth.
!
! The "within a single depth interval" functionality is partitioned over two modules:
!
! type_ersem_benthic_pom_layer is the master module that describes everything there is to know about
! the POM contents in a single depth interval. This is the module that can be instantiated from fabm.yaml.
! Its function is to collect interval-specific sink-source terms and translate those into the change in
! column-integrated POM and penetration depth.
!
! type_layer_content_calculator simply computes the mass with the specified depth interval.
! It is created automatically as a submodel of type_ersem_benthic_pom_layer, and does not
! interact with the user.
!
! The split of functionality across two modules is needed because the mass within the interval must be known
! early, so all other modules can use it, while the change in column-integrated mass and penetration depth
! must be computed late, after all other idividual modules have computed interval-specific rates of change.
! Thus we have the dependency chain:
!      1 mass within the desired depth interval (type_layer_content_calculator)
!   -> 2 interval-specific sink-source terms (all other modules)
!   -> 3 change in column-integrated mass and penetration depth (type_ersem_benthic_pom_layer)
! Putting 1 and 3 in the same module creates a circular dependency [that's BAD].

   use fabm_types
   use fabm_particle

   use ersem_shared
   use ersem_benthic_base

   implicit none

!  default: all is private.
   private

   ! Module for particulate organic matter class with idealized profile (e.g., Q6 Q7)
   type,extends(type_ersem_benthic_base),public :: type_ersem_benthic_column_particulate_matter
      type (type_bottom_state_variable_id) :: id_penetration_c,id_penetration_n,id_penetration_p,id_penetration_s
   contains
      procedure :: initialize
   end type

   ! Module for particulate organic matter within a single, user-specified depth interval.
   type,extends(type_particle_model),public :: type_ersem_benthic_pom_layer
      type (type_bottom_state_variable_id) :: id_c_int,id_n_int,id_p_int,id_s_int
      type (type_bottom_state_variable_id) :: id_pen_depth_c,id_pen_depth_n,id_pen_depth_p,id_pen_depth_s
      type (type_horizontal_dependency_id) :: id_c_sms,id_n_sms,id_p_sms,id_s_sms

      type (type_model_id) :: id_Q

      ! Layer extents
      integer :: surface_boundary_type   ! 0: constant depth, 1: dynamic depth from variable (e.g., depth of oxygenated layer)
      integer :: bottom_boundary_type    ! 0: constant depth, 1: dynamic depth from variable (e.g., depth of oxygenated layer)
      real(rk) :: surface_boundary_depth                                ! only used if surface_boundary_type==0
      real(rk) :: bottom_boundary_depth                                 ! only used if bottom_boundary_type==0
      type (type_horizontal_dependency_id) :: id_surface_boundary_depth ! only used if surface_boundary_type==1
      type (type_horizontal_dependency_id) :: id_bottom_boundary_depth  ! only used if bottom_boundary_type==1
   contains
      procedure :: initialize => layer_initialize
      procedure :: do_bottom  => layer_do_bottom
   end type

   type,extends(type_base_model) :: type_layer_content_calculator
      type (type_horizontal_diagnostic_variable_id) :: id_c
      type (type_horizontal_dependency_id) :: id_c_int,id_pen_depth, id_d_tot

      ! Layer extents
      integer :: surface_boundary_type   ! 0: constant depth, 1: dynamic depth from variable (e.g., depth of oxygenated layer)
      integer :: bottom_boundary_type    ! 0: constant depth, 1: dynamic depth from variable (e.g., depth of oxygenated layer)
      real(rk) :: surface_boundary_depth ! only used if surface_boundary_type==0
      real(rk) :: bottom_boundary_depth  ! only used if bottom_boundary_type==0
      type (type_horizontal_dependency_id) :: id_surface_boundary_depth ! only used if surface_boundary_type==1
      type (type_horizontal_dependency_id) :: id_bottom_boundary_depth  ! only used if bottom_boundary_type==1
   contains
      procedure :: do_bottom  => layer_content_calculator_do_bottom
   end type

contains
   
   subroutine initialize(self,configunit)
      class (type_ersem_benthic_column_particulate_matter), intent(inout), target :: self
      integer,                                     intent(in)            :: configunit

      class (type_ersem_benthic_pom_layer),pointer :: single_layer

      ! Perform normal benthic initialization (i.e., for POM without profile or penentration depth)
      call self%type_ersem_benthic_base%initialize(configunit)

      ! Add penetration depths for all active constituents.
      if (_VARIABLE_REGISTERED_(self%id_c)) call self%register_state_variable(self%id_penetration_c,'pen_depth_c','m','penetration depth of carbon')
      if (_VARIABLE_REGISTERED_(self%id_n)) call self%register_state_variable(self%id_penetration_n,'pen_depth_n','m','penetration depth of nitrogen')
      if (_VARIABLE_REGISTERED_(self%id_p)) call self%register_state_variable(self%id_penetration_p,'pen_depth_p','m','penetration depth of phosphorus')
      if (_VARIABLE_REGISTERED_(self%id_s)) call self%register_state_variable(self%id_penetration_s,'pen_depth_s','m','penetration depth of silicate')

      ! Create a submodel for particulate organic matter at the sediment surface
      ! This will receive sinks and sources associated with benthic-pelagic exchange - sedimentation, resuspension.
      allocate(single_layer)
      call single_layer%parameters%set('surface_boundary_type',0)
      call single_layer%parameters%set('bottom_boundary_type',0)
      call single_layer%parameters%set('surface_boundary_depth',0.0_rk)
      call single_layer%parameters%set('bottom_boundary_depth',0.0_rk)
      call single_layer%parameters%set('composition',self%composition)
      call single_layer%request_coupling('Q',self%name)
      call self%add_child(single_layer,'surface',configunit=configunit)
   end subroutine initialize

   subroutine layer_initialize(self,configunit)
   class (type_ersem_benthic_pom_layer), intent(inout), target :: self
   integer,                              intent(in)            :: configunit

   character(len=10) :: composition
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(composition,'composition', '','elemental composition',default='cnp')
      call self%get_parameter(self%surface_boundary_type,'surface_boundary_type','','surface boundary type (0: constant depth, 1: variable depth)',default=0)
      call self%get_parameter(self%bottom_boundary_type, 'bottom_boundary_type', '','bottom boundary type (0: constant depth, 1: variable depth)', default=0)
      if (self%surface_boundary_type==0) then
         call self%get_parameter(self%surface_boundary_depth,'surface_boundary_depth','m','surface boundary depth',default=0.0_rk)
      else
         call self%register_dependency(self%id_surface_boundary_depth,'surface_boundary_depth','m','surface boundary depth')
      end if
      if (self%bottom_boundary_type==0) then
         call self%get_parameter(self%bottom_boundary_depth,'bottom_boundary_depth','m','bottom boundary depth',default=0.0_rk)
      else
         call self%register_dependency(self%id_bottom_boundary_depth,'bottom_boundary_depth','m','bottom boundary depth')
      end if

      call self%register_model_dependency(self%id_Q,'Q')
      if (index(composition,'c')/=0) call layer_add_constituent(self,'c','mg C','carbon',    self%id_c_int,self%id_pen_depth_c,self%id_c_sms)
      if (index(composition,'n')/=0) call layer_add_constituent(self,'n','mmol','nitrogen',  self%id_n_int,self%id_pen_depth_n,self%id_n_sms)
      if (index(composition,'p')/=0) call layer_add_constituent(self,'p','mmol','phosphorus',self%id_p_int,self%id_pen_depth_p,self%id_p_sms)
      if (index(composition,'s')/=0) call layer_add_constituent(self,'s','mmol','silicate',  self%id_s_int,self%id_pen_depth_s,self%id_s_sms)

   end subroutine layer_initialize

   subroutine layer_add_constituent(self,name,units,long_name,id_c_int,id_pen_depth,id_sms)
      class (type_ersem_benthic_pom_layer), intent(inout),target :: self
      character(len=*),                     intent(in)           :: name,units,long_name
      type (type_bottom_state_variable_id), intent(inout),target :: id_c_int,id_pen_depth
      type (type_horizontal_dependency_id), intent(inout),target :: id_sms

      class (type_layer_content_calculator),pointer :: layer_content_calculator

      ! Register key inputs: column-integrated mass, penetration depth, and combined source-sink terms.
      call self%register_state_dependency(id_c_int,trim(name)//'_int',trim(units)//'/m^2', 'column-integrated '//trim(long_name))
      call self%register_state_dependency(id_pen_depth,'pen_depth_'//trim(name),'m','penetration depth for '//trim(long_name))
      call self%register_dependency(id_sms,trim(name)//'_sms',trim(units)//'/m^2/s', 'sinks-sources for '//trim(long_name))

      ! Create a module that will compute the mass in the desired depth interval.
      ! This module will also aggregate source-sink terms for thgis interval.
      allocate(layer_content_calculator)
      call self%add_child(layer_content_calculator,'content_calculator_'//trim(name),configunit=-1)

      ! Copy information on interval boundaries from master module.
      layer_content_calculator%surface_boundary_type  = self%surface_boundary_type
      layer_content_calculator%bottom_boundary_type   = self%bottom_boundary_type
      layer_content_calculator%surface_boundary_depth = self%surface_boundary_depth
      layer_content_calculator%bottom_boundary_depth  = self%surface_boundary_depth
      if (layer_content_calculator%surface_boundary_type==1) then
         ! Dynamic top boundary
         call layer_content_calculator%register_dependency(layer_content_calculator%id_surface_boundary_depth,'surface_boundary_depth','m','surface boundary depth')
         call layer_content_calculator%request_coupling(layer_content_calculator%id_surface_boundary_depth,'surface_boundary_depth')
      end if
      if (layer_content_calculator%bottom_boundary_type==1) then
         ! Dynamic bottom boundary
         call layer_content_calculator%register_dependency(layer_content_calculator%id_bottom_boundary_depth,'bottom_boundary_depth','m','bottom boundary depth')
         call layer_content_calculator%request_coupling(layer_content_calculator%id_bottom_boundary_depth,'bottom_boundary_depth')
      end if

      ! Get column-integrated mass and penetration depth from master module.
      call layer_content_calculator%register_dependency(layer_content_calculator%id_c_int,'c_int',trim(units)//'/m^2','column-integrated mass')
      call layer_content_calculator%register_dependency(layer_content_calculator%id_pen_depth,'pen_depth','m', 'penetration depth')
      call layer_content_calculator%request_coupling(layer_content_calculator%id_c_int,trim(name)//'_int')
      call layer_content_calculator%request_coupling(layer_content_calculator%id_pen_depth,'pen_depth_'//trim(name))
      call layer_content_calculator%register_dependency(layer_content_calculator%id_d_tot,'d_tot','m','depth of sediment column',standard_variable=depth_of_sediment_column)

      ! Register the only diagnostic exported by layer_content_calculator: mass integrated over desired depth interval.
      call layer_content_calculator%register_diagnostic_variable(layer_content_calculator%id_c,'c',trim(units)//'/m^2','mass')
      call layer_content_calculator%act_as_state_variable(layer_content_calculator%id_c)

      ! Link source-sink terms associated with layer-specific mass to the master module.
      call self%request_coupling(id_sms,'content_calculator_'//trim(name)//'/c_sms')

      ! Allow bulk coupling to a particulate organic matter module (rather than requiring the user to provide links per constituent and penetration depth)
      call self%request_coupling_to_model(id_c_int,self%id_Q,name)
      call self%request_coupling_to_model(id_pen_depth,self%id_Q,'pen_depth_'//trim(name))

      ! Create an alias in the master model for the layer-integrated density computed by layer_content_calculator.
      call self%add_horizontal_variable(name,trim(units)//'m^2','layer-integrated '//trim(long_name),domain=domain_bottom,act_as_state_variable=.true.)
      call self%request_coupling(name,'content_calculator_'//trim(name)//'/c')
   end subroutine layer_add_constituent

   subroutine layer_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_benthic_pom_layer), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: d_pen,d_top,d_bot,d_sms
      real(rk) :: c_int,sms

      _HORIZONTAL_LOOP_BEGIN_
         ! Determine top and bottom of desired depth interval.

         if (self%surface_boundary_type==0) then
            ! Constant top depth
            d_top = self%surface_boundary_depth
         else
            ! Variable top depth
            _GET_HORIZONTAL_(self%id_surface_boundary_depth,d_top)
         end if

         if (self%bottom_boundary_type==0) then
            ! Constant bottom depth
            d_bot = self%bottom_boundary_depth
         else
            ! Variable bottom depth
            _GET_HORIZONTAL_(self%id_bottom_boundary_depth,d_bot)
         end if

         ! Assume sinks-sources are homogenously distributed over desired depth interval.
         ! Thus, average depth of mass insertion/removal is the average of surface and bottom depths.
         d_sms = (d_top+d_bot)/2

         ! For each constituent: contribute to depth-integrated sink-source terms, contribute to change in penetration depth.
         if (_VARIABLE_REGISTERED_(self%id_c_int)) call layer_process_constituent_changes(self,_ARGUMENTS_LOCAL_,d_sms,self%id_c_int,self%id_pen_depth_c,self%id_c_sms)
         if (_VARIABLE_REGISTERED_(self%id_n_int)) call layer_process_constituent_changes(self,_ARGUMENTS_LOCAL_,d_sms,self%id_n_int,self%id_pen_depth_n,self%id_n_sms)
         if (_VARIABLE_REGISTERED_(self%id_p_int)) call layer_process_constituent_changes(self,_ARGUMENTS_LOCAL_,d_sms,self%id_p_int,self%id_pen_depth_p,self%id_p_sms)
         if (_VARIABLE_REGISTERED_(self%id_s_int)) call layer_process_constituent_changes(self,_ARGUMENTS_LOCAL_,d_sms,self%id_s_int,self%id_pen_depth_s,self%id_s_sms)
      _HORIZONTAL_LOOP_END_
   end subroutine layer_do_bottom

   subroutine layer_process_constituent_changes(self,_ARGUMENTS_LOCAL_,d_sms,id_c_int,id_pen_depth,id_sms)
      class (type_ersem_benthic_pom_layer), intent(in) :: self
      _DECLARE_ARGUMENTS_LOCAL_
      type (type_bottom_state_variable_id),intent(in) :: id_c_int,id_pen_depth
      type (type_horizontal_dependency_id),intent(in) :: id_sms
      real(rk),                            intent(in) :: d_sms

      real(rk) :: c_int,d_pen,sms

      ! Change in penetration depth can be derived by considering that the current penetration depth (d_pen)
      ! and mass density (c_int) are perturbed by addition of mass (delta_c) at some known depth (d_sms)
      ! New penetration depth = (d_pen*c_int + d_sms*delta_c)/(c_int+delta_c) = d_pen + delta_z
      ! delta_z = (d_pen*c_int + d_sms*delta_c)/(c_int+delta_c) - z_pen
      !         = (d_sms-z_pen)*delta_c/(c_int+delta_c)
      ! As we are considering change over an infinitesimal time, delta_c<<c_int, and we obtain
      ! delta_z = (d_sms-d_pen)*delta_c/c_int

      ! Retrieve depth-integrated mass, penetration depth, sinks-sources.
      _GET_HORIZONTAL_(id_c_int,c_int)
      _GET_HORIZONTAL_(id_pen_depth,d_pen)
      _GET_HORIZONTAL_(id_sms,sms)

      ! Apply sinks-sources to depth-integrated mass, compute change in penetration depth.
      _SET_BOTTOM_ODE_(id_c_int,sms)
      _SET_BOTTOM_ODE_(id_pen_depth, (d_sms-d_pen)*sms/c_int)
   end subroutine layer_process_constituent_changes

   subroutine layer_content_calculator_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_layer_content_calculator), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: d_pen
      real(rk) :: c_int
      real(rk) :: d_top, d_bot, d_totX

      _HORIZONTAL_LOOP_BEGIN_
         ! Retrieve column-integrated mass, penetration depth, and depth of modelled sediment column.
         _GET_HORIZONTAL_(self%id_c_int,c_int)
         _GET_HORIZONTAL_(self%id_pen_depth,d_pen)
         _GET_HORIZONTAL_(self%id_d_tot,d_totX)

         ! Determine top and bottom of desired depth interval.
         if (self%surface_boundary_type==0) then
            d_top = self%surface_boundary_depth
         else
            _GET_HORIZONTAL_(self%id_surface_boundary_depth,d_top)
         end if
         if (self%bottom_boundary_type==0) then
            d_bot = self%bottom_boundary_depth
         else
            _GET_HORIZONTAL_(self%id_bottom_boundary_depth,d_bot)
         end if

         ! Compute depth-integrated mass in desired depth interval.
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_c, c_int*partQ(d_pen, d_top, d_bot, d_totX))

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
      real(rk) :: norm, d_top1, d_bot1
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
      d_bot1 = min(d_bot, d_max)
      d_top1 = min(d_top, d_bot1)

      if ( d_pen>0._rk ) then
         ! Penetration depth > 0: integrate idealized [exponential] distribution over desired depth interval.

         ! Compute normalization factor: integral of exponential distribution from surface to bottom of column.
         ! This interval must by definition contain 100 % of the modelled mass.
         norm = 1._rk - exp(-d_max/d_pen)

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