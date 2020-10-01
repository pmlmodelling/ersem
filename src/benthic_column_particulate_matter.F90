#include "fabm_driver.h"

module ersem_benthic_column_particulate_matter

! Particulate organic matter with idealized [exponential] profile.
!
! -------------------------------------------------
! Shape of the profile
! -------------------------------------------------
!
! The concentration of matter, in mass per unit sediment volume, is assumed to be an
! exponential function of depth $z$:
!
!   C(z) = C0*exp(-b*z)
!
! What value do constants $C0$ and $b$ take?
!
! Let's define the mean depth of matter within entire sediment column (depth range $0$ to $\infty$):
!
!   z_mean = \int_0^\infty z C(z) dz / \int_0^\infty C(z) dz
!
! When we insert the exponential distribution of $C(z)$, $C0$ drops out and we obtain
!
!   z_mean = \int_0^\infty z exp(-b*z) dz / \int_0^\infty exp(-b*z) dz
!
! For the denominator of z_mean (concentration integrated from surface to infinite depth)
! we find through standard integration:
!
!   \int_0^\infty exp(-b*z) dz = [-1/b exp(-b*z)]_0^\infty = 1/b
!
! Numerator $\int_0^\infty z exp(-b*z) dz$ can be found through integration by parts:
! General rule: $\int u v dz = u \int v dz - \int (u' \int v dz) dz$
! We can apply this to find the antiderivative of $z exp(-b*z)$:
!
!   \int z exp(-b*z) dz = -z/b exp(-b*z) - 1/b^2 exp(-b*z) = -(z+1/b)/b exp(-b*z)
!
! Verify by differentiation (apply chain rule):
!
!   -1/b exp(-b*z) + z exp(-b*z) + 1/b exp(-b*z) = z exp(-b*z) [OK]
!
! Integrating from $0$ to $\intfy$ while assuming $b>0$ we obtain
!
!   \int_0^\infty z exp(-b*z) dz = [-(z+1/b)/b exp(-b*z)]_0^\infty = 1/b^2
!
! Combining the expressions for the numerator and denominator of $z_mean$:
!
!   z_mean = (1/b^2)/(1/b) = 1/b
!
! Thus, the exponential decay constant $b$ is equal to $1/z_mean$, and
!
!   C(z) = C0*exp(-z/z_mean)
!
! NB this is a standard result (mean of the exponential distribution); the above derivation
! is given simply for completeness.
!
! The integral of C(z) from 0 to the bottom of the modelled column, z_bot, should equal the modelled density of mass:
!
!   \int_0^z_bot C(z) dz = [-z_mean*C0*exp(-z/z_mean)]_0^z_bot = z_mean*C0*(1-exp(-z_bot/z_mean)) = C_int
!
! Thus, surface concentration $C0$ can be rewritten in terms of the depth-integrated concentration $C_int$,
! (integrated up to $z_bot$, not $\infty$!):
!
!   C0 = C_int/z_mean/(1-exp(-z_bot/z_mean))
!
! Note: Ebenhoeh et al. 1995 define the penetration depth as mean depth within the model domain,
! i.e., with lower boundary $z_bot$ rather than $\infty$:
!
!   z_mean' = \int_0^z_bot z C(z) dz / \int_0^z_bot C(z) dz
!           = [-(z+1/b)/b exp(-b*z)]_0^z_bot / [-1/b exp(-b*z)]_0^z_bot
!           = [1/b^2-(z_bot+1/b)/b exp(-b*z_bot)] / [1/b-1/b exp(-b*z_bot)]
!           = 1/b [1-(z_bot/b+1) exp(-b*z_bot)] / [1-exp(-b*z_bot)]
!           = 1/b [1-z_bot/b exp(-b*z_bot)/[1-exp(-b*z_bot)]
!
! This will NOT produce an explicit expression for b as a function of z_mean'.
! Instead, Ebenhoeh et al. implicitly assume z_bot>>z_mean', which leads to
! $b \approx 1/z_mean'$. However, this approximation becomes invalid when
! z_mean becomes large, and that can happen in the model. To avoid this problem,
! we define penetration depth as the mean depth of mass between 0 and \infty.
! Since our exponential function remains the same, expressions based on that
! (e.g., impact of bioturbation) are the same as in the Oldenburg model.
! Expressions based on the interpretation of penetration depth (the mean depth of mass
! between 0 and \infty in our case, and between 0 and z_bot for the Oldenburg model)
! will differ slightly, as e.g. in the next section.
!
! -------------------------------------------------
! Impact of sources and sinks at different depths
! -------------------------------------------------
!
! Sources and sinks change $C(z)$, and therefore also penetration depth
!
!   z_mean = \int_0^\infty z C(z) dz / \int_0^\infty C(z) dz
!
! The time derivative of this expression is found by applying the chain rule
!
!   d/dt z_mean = [\int_0^\infty z d/dt C(z) dz - \int_0^\infty d/dt C(z) dz \int_0^\infty z C(z) dz / \int_0^\infty C(z) dz] / \int_0^\infty C(z) dz
!
! Introducing depth-integrated concentration
!
!   C_int_\infty = \int_0^\infty C(z) dz,
!
! depth-integrated sources minus sinks
!
!   sms = \int_0^\infty d/dt C(z) dz,
!
! and the mean depth of the sources minus sinks
!
!   z_sms = \int_0^\infty z d/dt C(z) dz / \int_0^\infty d/dt C(z) dz,
!
! we can simplify this to
!
!   d/dt z_mean = (z_sms - z_mean) sms/C_int_\infty
!
! It is worth noting that the final division is by the concentration integrated from
! surface to *infinite depth*, $C_int_\infty$:
!
!   $C_int_\infty = z_mean C0
!
! That is not the same as $C_int$, i.e., the concentration integrated over the modelled
! depth interval ($0$ to $z_bot$), which equals
!
!    C_int = z_mean C0 (1-exp(-z_bot/z_mean))
!
! In the original Oldenburg implementation, penetration depth is defined as the mean depth of
! mass between 0 and z_bot (not 0 to \infty). As a result, $C_int$ is substituted for $C_int_\infty$.
! One consequence of this difference is that penetration depth in the Oldenburg model can "run away"
! (tend to infinity), while that is prevented in our expression due to the additional multiplication of the
! rate of change with 1-exp(-z_bot/z_mean).
!
! It is worth noting that the resulting expression for d/dt z_mean does NOT depend on distribution
! $C(z)$. That is, it is valid for ANY vertical distribution, exponential or otherwise.
!
! -------------------------------------------------
! Impact of bioturbation
! -------------------------------------------------
!
! Modelled as a diffusion process, bioturbation changes $C(z)$, and therefore also penetration depth
!
!   z_mean = \int_0^\infty z C(z) dz / \int_0^\infty C(z) dz
!
! Let us try this first for the denominator:
!
!   d/dt \int_0^\infty C(z) dz = \int_0^\infty d/dt C(z) dz
!
! For d/dt C(z) we have the normal diffusion equation:
!
!   d/dt C(z) = d/dz(D d/dz C(z)) = d/dz(-D/z_mean C0 exp(-z/z_mean)) = D C0/z_mean^2 exp(-z/z_mean)
!
! Inserting this expression we obtain for the change in depth-integrated mass:
!
!   d/dt \int_0^\infty C(z) dz = \int_0^\infty D C0/z_mean^2 exp(-z/z_mean) dz
!                              = [-D C0/z_mean exp(-z/z_mean)]_0^\infty
!                              = D C0/z_mean
!
! However, we KNOW that diffusion witin the column should not affect the mass integral. Why is this then non-zero?
! The reason for this is that we have not accounted for the no-flux boundary conditions. As a result, we are implicity
! using a non-zero inward flux at the surface that is determined by the gradient:
!
!   -D d/dz C(0) = D C0/z_mean exp(-z/z_mean) = D C0/z_mean
!
! In other words, to close the column for mass, we need to subtract this surface flux.
! (NB the gradient is zero at infinite depth, i.e., the bottom flux is zero)
!
! With this knowledge, we can revisit the numerator in the change in penetration depth $z_mean$.
! Its time derivative equals
!
!   d/dt \int_0^\infty z C(z) dz = \int_0^\infty z d/dt C(z) dz
!
! From the diffusion equation for $C(z)$ we derived $d/dt C(z) =  D C0/z_mean^2 exp(-z/z_mean)$.
! Inserting this in $d/dt z_mean$, while limiting bioturbation to maximum depth $z_tur$, we obtain
!
!   d/dt \int_0^z_tur z C(z) dz = D C0/z_mean^2 \int_0^z_tur z exp(-z/z_mean) dz
!
! For the expression within the integral we previously found antiderivative $-(z+z_mean)z_mean exp(-z/z_mean)$.
! Thus,
!
!   \int_0^z_tur z exp(-z/z_mean) dz = [-(z+z_mean)z_mean exp(-z/z_mean)]_0^z_tur
!                                    = -(z_tur+z_mean)z_mean exp(-z_tur/z_mean) + z_mean^2
!
! Using this antiderivative to solve the integral from $0$ to $z_tur$:
!
!    d/dt \int_0^\infty z C(z) dz = D C0 [1-(z_tur/z_mean+1)] exp(-z_tur/z_mean)
!
! This includes surface flux $-D C0 z/z_mean d/dz exp(-z/z_mean) = 0$, as well as non-zero bottom flux
! $-D C0 z_tur/z_mean exp(-z_tur/z_mean)$. Subtracting the latter we obtain:
!
!    d/dt \int_0^\infty z C(z) dz = D C0 [1 - exp(-z_tur/z_mean)]
!
! Finally, $d/dt z_mean$ is given by the ratio of this expression to $\int_0^\infty C(z) dz = C0 z_mean$:
!
!   d/dt z_mean = D/z_mean [1 - exp(-z_tur/z_mean)]
!
! This describes how the penetration depth $z_mean$ changes due to bioturbation up to depth $z_tur$, with
! the intensity of bioturbation described by diffusivity $D$. It is worth noting that $d/dt z_mean \to \infty$
! when $z_mean \to 0$. Thus, in the analytical solution of the model penetration depth cannot become negative
! as long as $D>0$.
!
! -------------------------------------------------
! Impact of burial
! -------------------------------------------------
!
! A change in penetration depth without an associated change in total mass $C_int_\infty$
! will cause a change in the mass in the depth interval described by the model, i.e.,
! \int_0^z_bot C0 exp(-z/z_mean). Inserting $C0 = C_int_\infty/z_mean$:
!
!   \int_0^z_bot C0 exp(-z/z_mean) = C_int_\infty/z_mean \int_0^z_bot exp(-z/z_mean)
!
! Taking the time derivative, with only z_mean depending on time (C_int_\infty is constant!):
!
!   d/dt \int_0^z_bot C0 exp(-z/z_mean) = d/dt z_mean -C_int_\infty/z_mean^2 \int_0^z_bot exp(-z/z_mean)
!                                         + C_int_\infty/z_mean \int_0^z_bot z/z_mean^2 d/dt z_mean exp(-z/z_mean) dz
!
! For the first term (first line), we have:
!
!    d/dt z_mean -C_int_\infty/z_mean^2 \int_0^z_bot exp(-z/z_mean) = - 1/z_mean d/dt z_mean C_int
!
! with C_int being the concentration integrated from surface to model bottom [not infinite depth!]
!
! We can rewrite the second term by moving constants outside the integral:
!
!    C_int_\infty/z_mean^3 d/dt z_mean \int_0^z_bot z exp(-z/z_mean) dz
!
! Now we can solve the integral with our previously found antiderivative:
!
!    \int_0^\z_bot z exp(-z/z_mean) dz = [-(z+z_mean)z_mean exp(-z/z_mean)]_0^z_bot
!                                      = -(z_bot+z_mean)z_mean exp(-z_bot/z_mean) + z_mean^2
!
! Inserting this
!
!    C_int_\infty/z_mean d/dt z_mean [-(z_bot/z_mean+1) exp(-z_bot/z_mean) + 1] = C_int_\infty/z_mean d/dt z_mean [1 - exp(-z_bot/z_mean) - z_bot/z_mean exp(-z_bot/z_mean)]
!
! Since $C_int = C_int_\infty [1 - exp(-z_bot/z_mean)]$, we can rewrite this as
!
!    1/z_mean d/dt z_mean C_int - C_int_\infty/z_mean d/dt z_mean z_bot/z_mean exp(-z_bot/z_mean)
!
! Combining the two contributions to d/dt \int_0^z_bot C0 exp(-z/z_mean), we are left with
!
!    d/dt \int_0^z_bot C0 exp(-z/z_mean) = - C_int_\infty z_bot/z_mean^2 exp(-z_bot/z_mean) d/dt z_mean
!
! JB 28/11/2014: this solution was validated in Python with sympy, using the following code
!
!    C_int_infty,z_bot,z,t = sympy.symbols('C_int_infty,z_bot,z,t',positive=True,real=True)
!    z_mean = sympy.Function('z_mean',positive=True,real=True)(t)
!    C = C_int_infty/z_mean*sympy.exp(-z/z_mean)
!    dC_dt = sympy.diff(C,t)
!    print sympy.integrate(dC_dt,(z,0,z_bot),conds='none')
!
! Inserting $C_int_\infty = C_int/[1 - exp(-z_bot/z_mean)]$
!
!    d/dt \int_0^z_bot C0 exp(-z/z_mean) = - C_int/[1 - exp(-z_bot/z_mean)] z_bot/z_mean^2 exp(-z_bot/z_mean) d/dt z_mean
!
! Since $C0 = C_int/[1 - exp(-z_bot/z_mean)]/z_mean$, this is equivalent to
!
!    d/dt \int_0^z_bot C0 exp(-z/z_mean) = C(z_bot) z_bot/z_mean d/dt z_mean
!
! In the Oldenburg code, the final expression for the change in mass within model domain due to burial is
!
!    - C_int/[1-exp(-z_bot/z_mean)] [exp(-(z_bot-d/dt z_mean)/z_mean-exp(-z_bot/z_mean)]
!
! Both our expression and the Oldenburg one contain -C_int/[1-exp(-z_bot/z_mean)] exp(-z_bot/z_mean).
! A discrepancy remains in the scale factor applied to it:
!
!    our derivation: z_bot/z_mean^2 d/dt z_mean
!    Oldenburg:      exp(d/dt z_mean/z_mean) - 1
!
! It's worth noting that the Oldenburg model here has a unit problem: the factor in the exponential
! is not dimensionless, but has unit time^-1.
!
! -------------------------------------------------
!
! This file contains two modules that can be instantiated by the user (from fabm.yaml):
! type_ersem_benthic_column_particulate_matter: describes particulate organic matter in terms of column-integrated mass and penetration depth
! type_ersem_benthic_pom_layer:                 describes particulate organic matter within a user-specified depth interval
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
! type_ersem_benthic_pom_layer computes the density of all constituents within the specified depth interval.
! This is the module that can be instantiated from fabm.yaml.
!
! type_constituent_for_single_layer_change collects interval-specific sink-source terms for a single constituent
! and translates those into the change in column-integrated matter and penetration depth. It is created
! automatically as a submodel of type_ersem_benthic_pom_layer and does not interact with the user.
!
! The split of functionality across two modules is needed because the mass within the interval must be known
! early, so all other modules can use it, while the change in column-integrated mass and penetration depth
! must be computed late, after all other individual modules have computed interval-specific rates of change.
! Thus we have the dependency chain:
!      1 mass within the desired depth interval (type_ersem_benthic_pom_layer)
!   -> 2 interval-specific sink-source terms (all other modules)
!   -> 3 change in column-integrated mass and penetration depth (type_constituent_for_single_layer_change)
! Putting 1 and 3 in the same module creates a circular dependency [that's BAD].

   use fabm_types
   use fabm_particle

   use ersem_shared
   use ersem_benthic_base

   implicit none

!  default: all is private.
   private

   ! Module for particulate organic matter class with idealized profile (e.g., Q6, Q7)
   ! Also handles the change in average penetration depth due to bioturbation, as well as burial in
   ! the form of organic matter moving beyond the bottom of the modelled column due to an increase
   ! in penetration depth. This module will always create a "surface" submodel for organic matter at
   ! the very surface of the sediment; this is typically used as target for sedimention/resuspension.
   type,extends(type_ersem_benthic_base),public :: type_ersem_benthic_column_particulate_matter
      type (type_bottom_state_variable_id) :: id_penetration_c,id_penetration_n,id_penetration_p,id_penetration_s
      type (type_bottom_state_variable_id) :: id_buried_c,id_buried_n,id_buried_p,id_buried_s
      type (type_bottom_state_variable_id) :: id_resuspendable_c,id_resuspendable_p,id_resuspendable_n,id_resuspendable_s
      type (type_horizontal_dependency_id) :: id_D, id_z_tur, id_d_tot, id_er

      logical :: burial
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

   ! Information for computing depth-integrated mass of a single constituent within a user-specified layer.
   type type_constituent_in_single_layer
      type (type_horizontal_diagnostic_variable_id) :: id_local     ! layer-integrated mass
      type (type_horizontal_dependency_id)          :: id_int       ! column-integrated mass
      type (type_horizontal_dependency_id)          :: id_pen_depth ! penetration depth
   end type

   ! Submodel that computes depth-integrated particulate organic matter within a single, user-specified depth interval.
   ! Top and bottom boundaries of this interval can be absolute depths or they can be derived from external variables.
   ! Supports remineralization, in which case the model must be coupled to sinks for each remineralized constituent.
   type,extends(type_particle_model),public :: type_ersem_benthic_pom_layer
      type (type_constituent_in_single_layer), allocatable :: constituents(:)
      type (type_horizontal_dependency_id) :: id_d_tot ! only used if variable_maximum_depth is set
      type (type_bottom_state_variable_id) :: id_remin_ox

      ! Layer extents
      logical :: variable_minimum_depth                        ! if not set, use minimum depth; if set, take dynamic minimum depth from variable (e.g., depth of oxygenated layer)
      logical :: variable_maximum_depth                        ! if not set, use maximum depth; if set, take dynamic maximum depth from variable (e.g., depth of oxygenated layer)
      real(rk) :: minimum_depth                                ! only used if variable_minimum_depth is not set
      real(rk) :: maximum_depth                                ! only used if variable_maximum_depth is not set
      type (type_horizontal_dependency_id) :: id_minimum_depth ! only used if variable_minimum_depth is set
      type (type_horizontal_dependency_id) :: id_maximum_depth ! only used if variable_maximum_depth is set
   contains
      procedure :: initialize => layer_initialize
      procedure :: do_bottom  => layer_do_bottom
   end type

   ! Submodel that converts layer-integrated sources of a single constituent into changes in depth-integrated mass and penetration depth.
   type,extends(type_base_model) ::  type_constituent_for_single_layer_change
      type (type_horizontal_dependency_id) :: id_sms           ! layer-integrated rate of change
      type (type_horizontal_dependency_id) :: id_local         ! layer-integrated mass
      type (type_bottom_state_variable_id) :: id_int           ! column-integrated mass
      type (type_bottom_state_variable_id) :: id_pen_depth     ! penetration depth
      type (type_bottom_state_variable_id) :: id_remin_target  ! target variable for remineralized matter
      type (type_bottom_state_variable_id) :: id_remin_ox      ! oxygen source for mineralisation
      type (type_horizontal_diagnostic_variable_id) :: id_remin_flux
      type(type_dependency_id)             :: id_ETW           ! temperature (controls remineralisation rate)
      type (type_horizontal_dependency_id) :: id_d_tot         ! depth of entire sediment column
      type (type_horizontal_dependency_id) :: id_pen_depth_c   ! penetration depth of carbon (only for source_depth_distribution==3!)

      ! Layer extents
      logical :: variable_minimum_depth                        ! if not set, use minimum depth; if set, take dynamic minimum depth from variable (e.g., depth of oxygenated layer)
      logical :: variable_maximum_depth                        ! if not set, use maximum depth; if set, take dynamic maximum depth from variable (e.g., depth of oxygenated layer)
      real(rk) :: minimum_depth                                ! only used if variable_minimum_depth is not set
      real(rk) :: maximum_depth                                ! only used if variable_maximum_depth is not set
      type (type_horizontal_dependency_id) :: id_minimum_depth ! only used if variable_minimum_depth is set
      type (type_horizontal_dependency_id) :: id_maximum_depth ! only used if variable_maximum_depth is set

      real(rk) :: remin
      real(rk) :: q10
      integer :: source_depth_distribution
   contains
      procedure :: do_bottom  => constituent_for_single_layer_change_do_bottom
   end type
   
contains

   subroutine initialize(self,configunit)
      class (type_ersem_benthic_column_particulate_matter), intent(inout), target :: self
      integer,                                              intent(in)            :: configunit

      character(len=10)                             :: composition
      real(rk)                                      :: c0
      class (type_ersem_benthic_pom_layer), pointer :: single_layer

      call self%initialize_ersem_benthic_base()

      call self%get_parameter(composition, 'composition', '', 'elemental composition')
      call self%get_parameter(c0,'c0','mg C/m^2','background carbon concentration',default=0.0_rk)

      call self%get_parameter(self%resuspension, 'resuspension', '', 'enable resuspension', default=.false.)

      if (index(composition,'c')/=0) call self%add_constituent('c',0.0_rk,c0)
      if (index(composition,'p')/=0) call self%add_constituent('p',0.0_rk,qpRPIcX*c0)
      if (index(composition,'n')/=0) call self%add_constituent('n',0.0_rk,qnRPIcX*c0)
      if (index(composition,'s')/=0) call self%add_constituent('s',0.0_rk,qsRPIcX*c0)
      if (index(composition,'f')/=0) call self%add_constituent('f',0.0_rk)

      ! Add penetration depths for all active constituents.
      if (_VARIABLE_REGISTERED_(self%id_c)) call self%register_state_variable(self%id_penetration_c,'pen_depth_c','m','penetration depth of carbon',minimum=0.0_rk)
      if (_VARIABLE_REGISTERED_(self%id_n)) call self%register_state_variable(self%id_penetration_n,'pen_depth_n','m','penetration depth of nitrogen',minimum=0.0_rk)
      if (_VARIABLE_REGISTERED_(self%id_p)) call self%register_state_variable(self%id_penetration_p,'pen_depth_p','m','penetration depth of phosphorus',minimum=0.0_rk)
      if (_VARIABLE_REGISTERED_(self%id_s)) call self%register_state_variable(self%id_penetration_s,'pen_depth_s','m','penetration depth of silicate',minimum=0.0_rk)

      ! Bioturbation
      ! Link to related standard variables (defined by ERSEM in shared.F90; not by FABM!)
      call self%register_dependency(self%id_D,particulate_diffusivity_due_to_bioturbation)
      call self%register_dependency(self%id_z_tur,bioturbation_depth)

      ! Total column depth is needed for burial formulation
      call self%register_dependency(self%id_d_tot,depth_of_sediment_column)

      ! Burial
      call self%get_parameter(self%burial,'burial','','enable burial',default=.false.)
      if (self%burial) then
         ! Add burial targets for all constituents
         if (_VARIABLE_REGISTERED_(self%id_c)) then
            call self%register_state_dependency(self%id_buried_c,'buried_c','mg C/m^2','buried carbon')
            call self%request_coupling_to_model(self%id_buried_c,'burial_target','c')
         end if
         if (_VARIABLE_REGISTERED_(self%id_p)) then
            call self%register_state_dependency(self%id_buried_p,'buried_p','mmol P/m^2','buried phosphorus')
            call self%request_coupling_to_model(self%id_buried_p,'burial_target','p')
         end if
         if (_VARIABLE_REGISTERED_(self%id_n)) then
            call self%register_state_dependency(self%id_buried_n,'buried_n','mmol N/m^2','buried nitrogen')
            call self%request_coupling_to_model(self%id_buried_n,'burial_target','n')
         end if
         if (_VARIABLE_REGISTERED_(self%id_s)) then
            call self%register_state_dependency(self%id_buried_s,'buried_s','mmol Si/m^2','buried silicate')
            call self%request_coupling_to_model(self%id_buried_s,'burial_target','s')
         end if
      end if

      if (self%resuspension) then
         ! Add dependencies on state variables that will be resuspended, and couple these to the surface layer of our own state variabes.
         ! This ensures that resuspension does not only update the biomass, but also the penetration depth.
         call self%register_dependency(self%id_er,sediment_erosion)
         if (_VARIABLE_REGISTERED_(self%id_c)) then
            call self%register_state_dependency(self%id_resuspendable_c,'resuspendable_c','mg C/m^2','source of resuspended carbon')
            call self%request_coupling(self%id_resuspendable_c,'surface/c')
         end if
         if (_VARIABLE_REGISTERED_(self%id_p)) then
            call self%register_state_dependency(self%id_resuspendable_p,'resuspendable_p','mmol P/m^2','source of resuspended phosphorus')
            call self%request_coupling(self%id_resuspendable_p,'surface/p')
         end if
         if (_VARIABLE_REGISTERED_(self%id_n)) then
            call self%register_state_dependency(self%id_resuspendable_n,'resuspendable_n','mmol N/m^2','source of resuspended nitrogen')
            call self%request_coupling(self%id_resuspendable_n,'surface/n')
         end if
         if (_VARIABLE_REGISTERED_(self%id_s)) then
            call self%register_state_dependency(self%id_resuspendable_s,'resuspendable_s','mmol Si/m^2','source of resuspended silicate')
            call self%request_coupling(self%id_resuspendable_s,'surface/s')
         end if
      end if

      ! Create a submodel for particulate organic matter at the sediment surface
      ! This will receive sinks and sources associated with benthic-pelagic exchange - sedimentation, resuspension.
      allocate(single_layer)
      call single_layer%parameters%set('variable_minimum_depth',.false.)
      call single_layer%parameters%set('variable_maximum_depth',.false.)
      call single_layer%parameters%set('minimum_depth',0.0_rk)
      call single_layer%parameters%set('maximum_depth',0.0_rk)
      call single_layer%parameters%set('composition',composition)
      call single_layer%couplings%set_string('Q',self%name)
      call self%add_child(single_layer,'surface',configunit=configunit)
   end subroutine initialize

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_benthic_column_particulate_matter), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: D,z_tur,z_bot
      real(rk) :: z_mean,z_mean_sms,burial_flux,C_int,v_er,resuspension_flux
      real(rk) :: max_pen_depth_change = 0.05_rk  ! max change in penetration depth due to bioturbation, in m/d
      real(rk) :: max_rel_res = 4.0_rk            ! max relative loss of matter due to resuspension, in 1/d

      _HORIZONTAL_LOOP_BEGIN_
         ! Get diffusivity and maximum depth of bioturbation
         _GET_HORIZONTAL_(self%id_D,D)
         _GET_HORIZONTAL_(self%id_z_tur,z_tur)
         _GET_HORIZONTAL_(self%id_d_tot,z_bot)

         ! Get sediment erosion rate in m/d
         if (self%resuspension) then
            _GET_HORIZONTAL_(self%id_er,v_er)
         end if

         ! Apply change in penetration depth due to bioturbation.
         ! See its derivation in the comments at the top of the file, section "Impact of bioturbation".
         ! If burial is active, also translate the change in penetration depth due to bioturbation to a burial
         ! flux at the bottom of the sediment column. See its derivation at the top of the file, section "Impact of burial".
         if (_VARIABLE_REGISTERED_(self%id_c)) then
            _GET_HORIZONTAL_(self%id_penetration_c,z_mean)
            z_mean_sms = min(D/z_mean,max_pen_depth_change)*(1.0_rk - exp(-z_tur/z_mean))
            _SET_BOTTOM_ODE_(self%id_penetration_c,z_mean_sms)
            if (self%burial.and.z_mean>0.0_rk) then
               _GET_HORIZONTAL_(self%id_c,C_int)
               burial_flux = C_int/(1.0_rk - exp(-z_bot/z_mean))*exp(-z_bot/z_mean)*z_bot/z_mean*z_mean_sms/z_mean
               _SET_BOTTOM_ODE_(self%id_c,-burial_flux)
               _SET_BOTTOM_ODE_(self%id_buried_c,burial_flux)
            end if
            if (self%resuspension) then
               _GET_HORIZONTAL_(self%id_c,C_int)
               resuspension_flux = C_int*min(max_rel_res, v_er/(z_mean+1.e-8_rk)/(1.0_rk - exp (-z_bot/z_mean)))
               _SET_BOTTOM_ODE_(self%id_resuspendable_c,-resuspension_flux)
               _SET_BOTTOM_EXCHANGE_(self%id_resuspended_c,resuspension_flux)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_resuspension_flux_c,resuspension_flux)
            end if
         end if
         if (_VARIABLE_REGISTERED_(self%id_p)) then
            _GET_HORIZONTAL_(self%id_penetration_p,z_mean)
            z_mean_sms = min(D/z_mean,max_pen_depth_change)*(1.0_rk - exp(-z_tur/z_mean))
            _SET_BOTTOM_ODE_(self%id_penetration_p,z_mean_sms)
            if (self%burial.and.z_mean>0.0_rk) then
               _GET_HORIZONTAL_(self%id_p,C_int)
               burial_flux = C_int/(1.0_rk - exp(-z_bot/z_mean))*exp(-z_bot/z_mean)*z_bot/z_mean*z_mean_sms/z_mean
               _SET_BOTTOM_ODE_(self%id_p,-burial_flux)
               _SET_BOTTOM_ODE_(self%id_buried_p,burial_flux)
            end if
            if (self%resuspension) then
               _GET_HORIZONTAL_(self%id_p,C_int)
               resuspension_flux = C_int*min(max_rel_res, v_er/(z_mean+1.e-8_rk)/(1.0_rk - exp (-z_bot/z_mean)))
               _SET_BOTTOM_ODE_(self%id_resuspendable_p,-resuspension_flux)
               _SET_BOTTOM_EXCHANGE_(self%id_resuspended_p,resuspension_flux)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_resuspension_flux_p,resuspension_flux)
            end if
         end if
         if (_VARIABLE_REGISTERED_(self%id_n)) then
            _GET_HORIZONTAL_(self%id_penetration_n,z_mean)
            z_mean_sms = min(D/z_mean,max_pen_depth_change)*(1.0_rk - exp(-z_tur/z_mean))
            _SET_BOTTOM_ODE_(self%id_penetration_n,z_mean_sms)
            if (self%burial.and.z_mean>0.0_rk) then
               _GET_HORIZONTAL_(self%id_n,C_int)
               burial_flux = C_int/(1.0_rk - exp(-z_bot/z_mean))*exp(-z_bot/z_mean)*z_bot/z_mean*z_mean_sms/z_mean
               _SET_BOTTOM_ODE_(self%id_n,-burial_flux)
               _SET_BOTTOM_ODE_(self%id_buried_n,burial_flux)
            end if
            if (self%resuspension) then
               _GET_HORIZONTAL_(self%id_n,C_int)
               resuspension_flux = C_int*min(max_rel_res, v_er/(z_mean+1.e-8_rk)/(1.0_rk - exp (-z_bot/z_mean)))
               _SET_BOTTOM_ODE_(self%id_resuspendable_n,-resuspension_flux)
               _SET_BOTTOM_EXCHANGE_(self%id_resuspended_n,resuspension_flux)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_resuspension_flux_n,resuspension_flux)
            end if
         end if
         if (_VARIABLE_REGISTERED_(self%id_s)) then
            _GET_HORIZONTAL_(self%id_penetration_s,z_mean)
            z_mean_sms = min(D/z_mean,max_pen_depth_change)*(1.0_rk - exp(-z_tur/z_mean))
            _SET_BOTTOM_ODE_(self%id_penetration_s,z_mean_sms)
            if (self%burial.and.z_mean>0.0_rk) then
               _GET_HORIZONTAL_(self%id_s,C_int)
               burial_flux = C_int/(1.0_rk - exp(-z_bot/z_mean))*exp(-z_bot/z_mean)*z_bot/z_mean*z_mean_sms/z_mean
               _SET_BOTTOM_ODE_(self%id_s,-burial_flux)
               _SET_BOTTOM_ODE_(self%id_buried_s,burial_flux)
            end if
            if (self%resuspension) then
               _GET_HORIZONTAL_(self%id_s,C_int)
               resuspension_flux = C_int*min(max_rel_res, v_er/(z_mean+1.e-8_rk)/(1.0_rk - exp (-z_bot/z_mean)))
               _SET_BOTTOM_ODE_(self%id_resuspendable_s,-resuspension_flux)
               _SET_BOTTOM_EXCHANGE_(self%id_resuspended_s,resuspension_flux)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_resuspension_flux_s,resuspension_flux)
            end if
         end if
      _HORIZONTAL_LOOP_END_
   end subroutine do_bottom

   subroutine layer_initialize(self,configunit)
      class (type_ersem_benthic_pom_layer), intent(inout), target :: self
      integer,                              intent(in)            :: configunit

      character(len=10) :: composition
      real(rk)          :: remin, q10
      integer           :: source_depth_distribution
      integer           :: iconstituent
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(composition,'composition', '','elemental composition')
      call self%get_parameter(self%variable_minimum_depth,'variable_minimum_depth','','minimum depth is variable', default=.false.)
      if (self%variable_minimum_depth) then
         call self%register_dependency(self%id_minimum_depth,'minimum_depth','m','minimum depth')
      else
         call self%get_parameter(self%minimum_depth,'minimum_depth','m','minimum depth',default=0.0_rk)
      end if
      call self%get_parameter(self%variable_maximum_depth,'variable_maximum_depth','','maximum depth is variable', default=.false.)
      if (self%variable_maximum_depth) then
         call self%register_dependency(self%id_maximum_depth,'maximum_depth','m','maximum depth')
      else
         call self%get_parameter(self%maximum_depth,'maximum_depth','m','maximum depth',default=0.0_rk)
      end if
      call self%get_parameter(remin,'remin','1/d','remineralization rate at 10 degrees Celsius',default=0.0_rk)
      if  (remin /= 0._rk) then
         call self%get_parameter(q10, 'q10', '-', 'Q_10 temperature coefficient', default=1.0_rk, minimum=1.0_rk)
      else
         q10=1.0_rk
      endif

      call self%get_parameter(source_depth_distribution,'source_depth_distribution', '','vertical distribution of changes (1: constant absolute rate, 2: constant relative rate, 3: constant carbon-based relative rate)',default=1)

      call self%register_dependency(self%id_d_tot,depth_of_sediment_column)

      allocate(self%constituents(len_trim(composition)))
      do iconstituent=1,len_trim(composition)
         select case (composition(iconstituent:iconstituent))
         case ('c')
            call layer_initialize_constituent(self,self%constituents(iconstituent),'c','mg C',   'carbon',    remin,q10,source_depth_distribution,standard_variables%total_carbon,1.0_rk/CMass)
         case ('n')
            call layer_initialize_constituent(self,self%constituents(iconstituent),'n','mmol N', 'nitrogen',  remin,q10,source_depth_distribution,standard_variables%total_nitrogen)
         case ('p')
            call layer_initialize_constituent(self,self%constituents(iconstituent),'p','mmol P', 'phosphorus',remin,q10,source_depth_distribution,standard_variables%total_phosphorus)
         case ('s')
            call layer_initialize_constituent(self,self%constituents(iconstituent),'s','mmol Si','silicate',  remin,q10,source_depth_distribution,standard_variables%total_silicate)
         case ('f')
         case default
            call self%fatal_error('layer_initialize','unknown constituent '//composition(iconstituent:iconstituent)//' specified.')
         end select
      end do

   end subroutine layer_initialize

   subroutine layer_initialize_constituent(self,info,name,units,long_name,remin,q10,source_depth_distribution,aggregate_target,aggregate_scale_factor)
      class (type_ersem_benthic_pom_layer),     intent(inout), target :: self
      class (type_constituent_in_single_layer), intent(inout)         :: info
      character(len=*),                         intent(in)            :: name,units,long_name
      real(rk),                                 intent(in)            :: remin,q10
      integer,                                  intent(in)            :: source_depth_distribution
      type (type_bulk_standard_variable),       intent(in)            :: aggregate_target
      real(rk),optional,                        intent(in)            :: aggregate_scale_factor

      class (type_constituent_for_single_layer_change), pointer :: change_processor

      ! Register the only diagnostic for this constituent: mass integrated over desired depth interval.
      ! This diagnostic acts like a state variable, so that other models can provides sinks and sources.
      ! These are converted by the "change_processor" submodel into appropriate changes in colmn-integrated mass and penetration depth.
      call self%register_diagnostic_variable(info%id_local,name,units//'/m^2',long_name, &
         act_as_state_variable=.true.,domain=domain_bottom,output=output_none,source=source_do_bottom)
      call self%add_to_aggregate_variable(aggregate_target,info%id_local,aggregate_scale_factor)

      ! Register dependencies on column-integrated mass and penetration depth.
      ! Make sure these can only be linked to state variables, as both mass and penetration depth will change dynamically.
      call self%register_dependency(info%id_int,name//'_int',units//'/m^2','column-integrated '//long_name)
      call self%register_dependency(info%id_pen_depth,'pen_depth_'//name,'m', 'penetration depth of '//long_name)
      info%id_int%link%target%fake_state_variable = .true.
      info%id_pen_depth%link%target%fake_state_variable = .true.

      ! Allow bulk coupling to a particulate organic matter module
      ! (rather than requiring the user to provide links per constituent and penetration depth)
      call self%request_coupling_to_model(info%id_int,'Q',name)
      call self%request_coupling_to_model(info%id_pen_depth,'Q','pen_depth_'//name)

      ! Create a submodel that converts layer-integrated sources into a change in column-integrated mass and penetration depth.
      allocate(change_processor)
      change_processor%dt = 86400._rk
      change_processor%remin = remin
      change_processor%q10 = q10
      change_processor%source_depth_distribution = source_depth_distribution
      call self%add_child(change_processor,'change_in_'//name//'_processor',configunit=-1)
      call change_processor%register_state_dependency(change_processor%id_int,name//'_int',units//'/m^2', 'column-integrated '//long_name)
      call change_processor%register_state_dependency(change_processor%id_pen_depth,'pen_depth_'//name,'m','penetration depth for '//long_name)
      call change_processor%register_dependency(change_processor%id_sms,name//'_sms',units//'/m^2/s', 'sinks-sources for '//long_name)

      ! Copy information about layer boundaries to submodel.
      change_processor%variable_minimum_depth = self%variable_minimum_depth
      change_processor%variable_maximum_depth = self%variable_maximum_depth
      change_processor%minimum_depth = self%minimum_depth
      change_processor%maximum_depth = self%maximum_depth
      if (change_processor%variable_minimum_depth) then
         ! Dynamic top boundary
         call change_processor%register_dependency(change_processor%id_minimum_depth,'minimum_depth','m','minimum depth')
         call change_processor%request_coupling(change_processor%id_minimum_depth,'../minimum_depth')
      end if
      if (change_processor%variable_maximum_depth) then
         ! Dynamic bottom boundary
         call change_processor%register_dependency(change_processor%id_maximum_depth,'maximum_depth','m','maximum depth')
         call change_processor%request_coupling(change_processor%id_maximum_depth,'../maximum_depth')
      end if
      call change_processor%register_dependency(change_processor%id_d_tot,depth_of_sediment_column)

      ! Link source-sink terms associated with layer-specific mass to the master module.
      call change_processor%request_coupling(change_processor%id_int,      '../'//name//'_int')
      call change_processor%request_coupling(change_processor%id_pen_depth,'../pen_depth_'//name)
      call change_processor%request_coupling(change_processor%id_sms,      '../'//name//'_sms_tot')

      if (legacy_ersem_compatibility.and.source_depth_distribution==3) then
         ! Source location depends on penetration depth of carbon - register a dependency on it.
         call change_processor%register_dependency(change_processor%id_pen_depth_c,'pen_depth_c','m','penetration depth of carbon')
         call change_processor%request_coupling(change_processor%id_pen_depth_c,'../pen_depth_c')
      end if

      ! Register a dependency on the diagnostic that holds the layer-integrated density.
      ! This will be used to compute remineralization
      if (remin/=0.0_rk) then
         ! Register submodel dependencies: layer-integrated mass (subject to remineralization) and sink for remineralized matter.
         call change_processor%register_dependency(change_processor%id_local,name,units//'/m^2','layer-integrated '//long_name)
         call change_processor%register_state_dependency(change_processor%id_remin_target,name//'_remin_target',units//'/m^2','sink for remineralized '//long_name)
         if (name=='c') then 
             call self%register_state_dependency(self%id_remin_ox,'o_remin_source','mmol O_2/m^2','oxygen')
             call change_processor%register_state_dependency(change_processor%id_remin_ox,'o_remin_source','mmol O_2/m^2','oxygen')
             call change_processor%request_coupling(change_processor%id_remin_ox, '../o_remin_source')
         end if
         call change_processor%register_diagnostic_variable(change_processor%id_remin_flux,'remin_flux',units//'/m^2/d','mineralisation flux',source=source_do_bottom)
         ! Couple submodel dependencies to top-level ("self") equivalents.
         ! For the remineralization target, register an alias at the top level so that it can be coupled directly from fabm.yaml.
         call change_processor%request_coupling(change_processor%id_local,'../'//name)
         call self%add_horizontal_variable(name//'_remin_target',units//'/m^2','sink for remineralized '//long_name,act_as_state_variable=.true.)
         call change_processor%request_coupling(change_processor%id_remin_target,'../'//name//'_remin_target')
         call change_processor%register_dependency(change_processor%id_ETW, standard_variables%temperature)
      end if

   end subroutine layer_initialize_constituent

   subroutine layer_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_benthic_pom_layer), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      integer :: iconstituent

      do iconstituent=1,size(self%constituents)
         if (_VARIABLE_REGISTERED_(self%constituents(iconstituent)%id_local)) &
            call layer_process_constituent(self,_ARGUMENTS_DO_BOTTOM_,self%constituents(iconstituent))
      end do
   end subroutine layer_do_bottom

   subroutine layer_process_constituent(self,_ARGUMENTS_DO_BOTTOM_,info)
      class (type_ersem_benthic_pom_layer),   intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
      type (type_constituent_in_single_layer),intent(in) :: info

      real(rk) :: c_int
      real(rk) :: d_pen
      real(rk) :: d_min, d_max, d_totX

      _HORIZONTAL_LOOP_BEGIN_

         ! Retrieve column-integrated mass, penetration depth, and depth of modelled sediment column.
         _GET_HORIZONTAL_(info%id_int,c_int)
         _GET_HORIZONTAL_(info%id_pen_depth,d_pen)
         _GET_HORIZONTAL_(self%id_d_tot,d_totX)

         ! Determine top and bottom of desired depth interval.
         if (self%variable_minimum_depth) then
            _GET_HORIZONTAL_(self%id_minimum_depth,d_min)
         else
            d_min = self%minimum_depth
         end if
         if (self%variable_maximum_depth) then
            _GET_HORIZONTAL_(self%id_maximum_depth,d_max)
         else
            d_max = self%maximum_depth
         end if

         ! Compute depth-integrated mass within desired depth interval.
         _SET_HORIZONTAL_DIAGNOSTIC_(info%id_local, c_int*partQ(d_pen, d_min, d_max, d_totX))

      _HORIZONTAL_LOOP_END_
   end subroutine layer_process_constituent

   subroutine constituent_for_single_layer_change_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_constituent_for_single_layer_change), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: d_min,d_max
      real(rk) :: c_int,d_pen,d_pen_c,sms,d_sms
      real(rk) :: c_int_local,sms_remin, d_tot
      real(rk),parameter :: max_relax = 4.0_rk   ! max relaxation rate for penetration depth (1/d)
      real(rk) :: relax,eT,ETW

      _HORIZONTAL_LOOP_BEGIN_
         ! Determine top and bottom of desired depth interval.

         if (self%variable_minimum_depth) then
            ! Minimum depth from coupled variable.
            _GET_HORIZONTAL_(self%id_minimum_depth,d_min)
         else
            ! Constant minimum depth
            d_min = self%minimum_depth
         end if

         if (self%variable_maximum_depth) then
            ! Maximum depth from coupled variable.
            _GET_HORIZONTAL_(self%id_maximum_depth,d_max)
         else
            ! Constant maximum depth
            d_max = self%maximum_depth
         end if

         ! Ensure that d_max >= d_min:
         d_max = max(d_max,d_min)

         ! Retrieve depth-integrated mass, penetration depth, sinks-sources.
         _GET_HORIZONTAL_(self%id_int,c_int)
         _GET_HORIZONTAL_(self%id_pen_depth,d_pen)
         _GET_HORIZONTAL_(self%id_sms,sms)

         ! Convert sources-sinks in per second, as returned by FABM, to our own time unit.
         sms = sms*self%dt

         ! Add local remineralization
         if (self%remin/=0.0_rk) then
            _GET_(self%id_ETW, ETW)
            eT  = self%q10**((ETW-10._rk)/10._rk)
            _GET_HORIZONTAL_(self%id_local,c_int_local)
            sms_remin = self%remin*c_int_local*eT
            _SET_BOTTOM_ODE_(self%id_remin_target,sms_remin)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_remin_flux,sms_remin)
            if (_VARIABLE_REGISTERED_(self%id_remin_ox)) _SET_BOTTOM_ODE_(self%id_remin_ox,-sms_remin/CMass)
            sms = sms - sms_remin
         end if

         if (self%source_depth_distribution==1) then
            ! Assume sinks-sources are homogenously distributed over desired depth interval.
            ! That is, d/dt C(z) within the layer does not depend on depth.
            ! Thus, average depth of mass insertion/removal is the average of surface and bottom depths.
            d_sms = (d_min+d_max)/2
         else
            ! Assume the relative rate of change in mass is depth-independent:
            ! 1/C(z) d/dt C(z) within the active layer does not depend on depth.
            ! Thus we can define source terms as d/dt C(z) = r C(z)
            ! The mean depth of these source terms then equals
            !    \int_d_min^d_max z d/dt C(z) dz / \int_d_min^d_max d/dt C(z) dz
            ! As d/dt C(z) = r C(z), this is equal to
            !    \int_d_min^d_max z C(z) dz / \int_d_min^d_max C(z) dz
            ! In words, the mean depth of source terms within the specified depth interval is equal
            ! to the mean depth of matter within the same depth interval.
            !
            ! Now we need to find the mean depth of matter. With C(z) = C0*exp(-z/z_mean), this can be written as
            !    \int_d_min^d_max exp(-z/z_mean) z dz / \int_d_min^d_max exp(-z/z_mean) dz
            ! Anti-derivatives of these expressions are derived near the top of this file.
            ! For the denominator we found
            !    \int_d_min^d_max exp(-z/z_mean) dz = [-z_mean exp(-z/z_mean)]_d_min^\d_max
            !                                       = z_mean [exp(-d_min/z_mean)-exp(-d_max/z_mean)]
            ! For the numerator we found
            !    \int_d_min^d_max z exp(-z/z_mean) dz = [-(z+z_mean)z_mean exp(-z/z_mean)]_d_min^\d_max
            !                                         = z_mean [(d_min+z_mean) exp(-d_min/z_mean)-(d_max+z_mean) exp(-d_max/z_mean)]
            !                                         = z_mean^2 [exp(-d_min/z_mean)-exp(-d_max/z_mean)] + z_mean [d_min exp(-d_min/z_mean)-d_max exp(-d_max/z_mean)]
            ! Combining numerator and denominator, we obtain
            !    z_mean + [d_min exp(-d_min/z_mean)-d_max exp(-d_max/z_mean)] / [exp(-d_min/z_mean)-exp(-d_max/z_mean)]
            ! This can be rearranged to
            !    z_mean + d_min - (d_max-d_min)/(exp((d_max-d_min)/z_mean)-1)
            if (legacy_ersem_compatibility) then
               ! If $d_max-d_min>>z_mean$, the above expression simplifies to z_mean + d_min.
               ! This is the expression used in legacy ERSEM, which thus assumes that
               ! bottom interface of the active layer is much deeper than the penetration depth.
               if (self%source_depth_distribution==2) then
                  d_sms = d_min + d_pen
               elseif (self%source_depth_distribution==3) then
                  ! As for self%source_depth_distribution==2, but derive depth of change from
                  ! the penetration depth of carbon, even when changing nitrogen, phosphorous, silicate.
                  ! This was done in legacy ERSEM, but likely only because Q?Distribution routines
                  ! take a single penetration depth argument that is then applied to all constituents.
                  _GET_HORIZONTAL_(self%id_pen_depth_c,d_pen_c)
                  d_sms = d_min + d_pen_c
               end if
            else
               if ( (d_max-d_min)/d_pen > ZeroX ) then
                  d_sms = d_min + d_pen - (d_max-d_min)/(exp((d_max-d_min)/d_pen)-1)
               else
                  ! Safety valve for pathological cases where d_max == d_min
                  ! (dmax >= d_min is ensured above)
                  d_sms = d_min
               end if
            end if
         end if

         ! Compute change in penetration depth. See its derivation in the comments at the top of the file,
         ! section "Impact of sources and sinks at different depths".
         _GET_HORIZONTAL_(self%id_d_tot,d_tot)
         if (legacy_ersem_compatibility) then
            ! Legacy ERSEM clips d_sms to avoid positive feedback that leads to $d_pen->\infty$
            ! This feedback very like stems from the assumption that d_pen<<d_tot.
            if (d_pen>d_tot .and. sms<0._rk) d_sms = d_tot
            _SET_BOTTOM_ODE_(self%id_pen_depth, (d_sms-d_pen)*sms/c_int)
         else
            ! New formulation defines d_pen explicitly as mean depth between 0 and \infty, rather than
            ! between 0 and d_tot [see derivation near top of file]. Thus, it does not assume d_pen<<d_tot.
            ! The result of this is that C_int_infty (not C_int) appears in the denominator.
            ! This should already prevent a positive feedback, since the last term in the expression below
            ! ensures that the change goes to zero when d_pen>>d_tot. When d_pen<<d_tot, the old result
            ! is recovered. JornB 11/3/2015

            ! In the expression below, we protect against extreme changes in penetration
            ! depth at small c_int. Note per the analytical derivation, relax = sms/c_int.
            if (sms==0.0_rk) then
               ! Separately handle sms==0, in order to prevent 0/0 when sms and c_int are both 0.
               relax = 0.0_rk
            else
               ! The expression below handles c_int==0 gracefully provided the user does not
               ! activate runtime floating point exceptions that trap division by 0.
               relax = sign(min(abs(sms)/c_int,max_relax),sms)
            end if
            _SET_BOTTOM_ODE_(self%id_pen_depth, (d_sms-d_pen)*relax*(1-exp(-d_tot/d_pen)))
         end if

         ! Apply sinks-sources to depth-integrated mass
         _SET_BOTTOM_ODE_(self%id_int,sms)
      _HORIZONTAL_LOOP_END_
   end subroutine constituent_for_single_layer_change_do_bottom

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
