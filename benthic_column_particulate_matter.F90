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
! The mean depth of matter equals:
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
! In the original Oldenburg implementation, the difference between these quantities ignored.
! That is, in $d/dt z_mean$, $C_int$ is substituted for $C_int_\infty$.
!
! Worth noting that the resulting expression for d/dt z_mean does NOT depend on distribution
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
      type (type_bottom_state_variable_id) :: id_buried_c,id_buried_n,id_buried_p,id_buried_s
      type (type_horizontal_dependency_id) :: id_D, id_z_tur
      type (type_model_id)                 :: id_buried_Q
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

   type,extends(type_base_model) :: type_burial
      type (type_bottom_state_variable_id) :: id_mass,id_buried_mass
      type (type_horizontal_dependency_id) :: id_penetration_depth,id_penetration_depth_sms,id_d_tot
   contains
      procedure :: do_bottom => burial_do_bottom
   end type

   ! Module for particulate organic matter within a single, user-specified depth interval.
   ! Supports remineralization, in which case the model must be coupled to sinks for each remineralized constituent.
   type,extends(type_particle_model),public :: type_ersem_benthic_pom_layer
      type (type_bottom_state_variable_id) :: id_c_int,id_n_int,id_p_int,id_s_int
      type (type_bottom_state_variable_id) :: id_pen_depth_c,id_pen_depth_n,id_pen_depth_p,id_pen_depth_s
      type (type_horizontal_dependency_id) :: id_c_sms,id_n_sms,id_p_sms,id_s_sms
      type (type_bottom_state_variable_id) :: id_c_remin_target,id_n_remin_target,id_p_remin_target,id_s_remin_target
      type (type_horizontal_dependency_id) :: id_c_local,id_n_local,id_p_local,id_s_local

      type (type_model_id) :: id_Q

      real(rk) :: remin

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
      integer,                                              intent(in)            :: configunit

      class (type_ersem_benthic_pom_layer),pointer :: single_layer
      logical                                      :: bury

      ! Perform normal benthic initialization (i.e., for POM without profile or penentration depth)
      call self%type_ersem_benthic_base%initialize(configunit)

      ! Add penetration depths for all active constituents.
      if (_VARIABLE_REGISTERED_(self%id_c)) call self%register_state_variable(self%id_penetration_c,'pen_depth_c','m','penetration depth of carbon')
      if (_VARIABLE_REGISTERED_(self%id_n)) call self%register_state_variable(self%id_penetration_n,'pen_depth_n','m','penetration depth of nitrogen')
      if (_VARIABLE_REGISTERED_(self%id_p)) call self%register_state_variable(self%id_penetration_p,'pen_depth_p','m','penetration depth of phosphorus')
      if (_VARIABLE_REGISTERED_(self%id_s)) call self%register_state_variable(self%id_penetration_s,'pen_depth_s','m','penetration depth of silicate')

      ! Bioturbation
      ! Link to related standard variables (defined by ERSEM in shared.F90; not by FABM!)
      call self%register_dependency(self%id_D,particulate_diffusivity_due_to_bioturbation)
      call self%register_dependency(self%id_z_tur,bioturbation_depth)

      ! Burial
      call self%get_parameter(bury,'burial','','activate burial',default=.false.)
      if (bury) then
         ! Enable wholesale linking to model "burial_target", which would hook up all buried mass constituents.
         call self%register_model_dependency(self%id_buried_Q,'burial_target')

         ! Add burial modules for all constituents
         if (_VARIABLE_REGISTERED_(self%id_c)) call add_burial(self,self%id_c,self%id_penetration_c,self%id_buried_c)
         if (_VARIABLE_REGISTERED_(self%id_n)) call add_burial(self,self%id_n,self%id_penetration_n,self%id_buried_n)
         if (_VARIABLE_REGISTERED_(self%id_p)) call add_burial(self,self%id_p,self%id_penetration_p,self%id_buried_p)
         if (_VARIABLE_REGISTERED_(self%id_s)) call add_burial(self,self%id_s,self%id_penetration_s,self%id_buried_s)
      end if

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

   subroutine add_burial(self,id_mass,id_penetration,id_buried_mass)
      class (type_ersem_benthic_column_particulate_matter), intent(inout), target :: self
      type (type_bottom_state_variable_id),                 intent(in)            :: id_mass,id_penetration
      type (type_bottom_state_variable_id),                 intent(out)           :: id_buried_mass

      class (type_burial),pointer :: burial

      ! Create a dependency on the POM level for buried mass.
      ! That will enable coupling to be specified on this level.
      call self%register_state_dependency(id_buried_mass,'buried_'//trim(id_mass%link%name),id_mass%link%target%units,'buried '//trim(id_mass%link%target%long_name))
      call self%request_coupling_to_model(id_buried_mass,self%id_buried_Q,id_mass%link%name)

      allocate(burial)
      call self%add_child(burial,trim(id_mass%link%name)//'_burial',configunit=-1)

      call burial%register_state_dependency(burial%id_mass,'mass',id_mass%link%target%units,'mass')
      call burial%register_state_dependency(burial%id_buried_mass,'buried_mass',id_mass%link%target%units,'buried mass')
      call burial%register_dependency(burial%id_penetration_depth,'penetration_depth','m','penetration depth')
      call burial%register_dependency(burial%id_penetration_depth_sms,'penetration_depth_sms','m/s','change in penetration depth')
      call burial%register_dependency(burial%id_d_tot,depth_of_sediment_column)

      ! Link to variables in parent model
      call burial%request_coupling(burial%id_buried_mass,'buried_'//trim(id_mass%link%name))
      call burial%request_coupling(burial%id_mass,id_mass%link%name)
      call burial%request_coupling(burial%id_penetration_depth,id_penetration%link%name)
      call burial%request_coupling(burial%id_penetration_depth_sms,trim(id_penetration%link%name)//'_sms')
   end subroutine

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_benthic_column_particulate_matter), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: D,z_mean,z_tur

      _HORIZONTAL_LOOP_BEGIN_
         ! Get diffusivity and maximum depth of bioturbation
         _GET_HORIZONTAL_(self%id_D,D)
         _GET_HORIZONTAL_(self%id_z_tur,z_tur)

         ! Apply change in penentration depth due to bioturbation.
         ! See its derivation in the comments at the top of the file,
         ! section "Impact of sources and sinks at different depths".
         if (_VARIABLE_REGISTERED_(self%id_c)) then
            _GET_HORIZONTAL_(self%id_penetration_c,z_mean)
            _SET_BOTTOM_ODE_(self%id_penetration_c,D/z_mean*(1.0_rk - exp(-z_tur/z_mean)))
         end if
         if (_VARIABLE_REGISTERED_(self%id_n)) then
            _GET_HORIZONTAL_(self%id_penetration_n,z_mean)
            _SET_BOTTOM_ODE_(self%id_penetration_n,D/z_mean*(1.0_rk - exp(-z_tur/z_mean)))
         end if
         if (_VARIABLE_REGISTERED_(self%id_p)) then
            _GET_HORIZONTAL_(self%id_penetration_p,z_mean)
            _SET_BOTTOM_ODE_(self%id_penetration_p,D/z_mean*(1.0_rk - exp(-z_tur/z_mean)))
         end if
         if (_VARIABLE_REGISTERED_(self%id_s)) then
            _GET_HORIZONTAL_(self%id_penetration_s,z_mean)
            _SET_BOTTOM_ODE_(self%id_penetration_s,D/z_mean*(1.0_rk - exp(-z_tur/z_mean)))
         end if
      _HORIZONTAL_LOOP_END_
   end subroutine do_bottom

   subroutine burial_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_burial), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: z_mean,z_bot,C_int,z_mean_sms,burial_flux

      _HORIZONTAL_LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_d_tot,z_bot)
         _GET_HORIZONTAL_(self%id_mass,C_int)
         _GET_HORIZONTAL_(self%id_penetration_depth,z_mean)
         _GET_HORIZONTAL_(self%id_penetration_depth_sms,z_mean_sms)
         burial_flux = C_int/(1.0_rk - exp(-z_bot/z_mean))*exp(-z_bot/z_mean)*z_bot/z_mean*z_mean_sms/z_mean
         _SET_BOTTOM_ODE_(self%id_mass,-burial_flux)
         _SET_BOTTOM_ODE_(self%id_buried_mass,burial_flux)
      _HORIZONTAL_LOOP_END_
   end subroutine burial_do_bottom

   subroutine layer_initialize(self,configunit)
      class (type_ersem_benthic_pom_layer), intent(inout), target :: self
      integer,                              intent(in)            :: configunit

      character(len=10) :: composition
!EOP
!-----------------------------------------------------------------------
!BOC
      self%dt = 86400._rk

      call self%get_parameter(composition,'composition', '','elemental composition',default='cnp')
      call self%get_parameter(self%remin,'remin','d-1','remineralization rate',default=0.0_rk)
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
      if (index(composition,'c')/=0) call layer_add_constituent(self,'c','mg C','carbon',    self%id_c_int,self%id_pen_depth_c,self%id_c_sms,self%id_c_local,self%id_c_remin_target)
      if (index(composition,'n')/=0) call layer_add_constituent(self,'n','mmol','nitrogen',  self%id_n_int,self%id_pen_depth_n,self%id_n_sms,self%id_n_local,self%id_n_remin_target)
      if (index(composition,'p')/=0) call layer_add_constituent(self,'p','mmol','phosphorus',self%id_p_int,self%id_pen_depth_p,self%id_p_sms,self%id_p_local,self%id_p_remin_target)
      if (index(composition,'s')/=0) call layer_add_constituent(self,'s','mmol','silicate',  self%id_s_int,self%id_pen_depth_s,self%id_s_sms,self%id_s_local,self%id_s_remin_target)

   end subroutine layer_initialize

   subroutine layer_add_constituent(self,name,units,long_name,id_c_int,id_pen_depth,id_sms,id_local,id_remin_target)
      class (type_ersem_benthic_pom_layer), intent(inout),target :: self
      character(len=*),                     intent(in)           :: name,units,long_name
      type (type_bottom_state_variable_id), intent(inout),target :: id_c_int,id_pen_depth,id_remin_target
      type (type_horizontal_dependency_id), intent(inout),target :: id_sms,id_local

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
      layer_content_calculator%bottom_boundary_depth  = self%bottom_boundary_depth
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
      call layer_content_calculator%register_diagnostic_variable(layer_content_calculator%id_c,'c',trim(units)//'/m^2','mass',act_as_state_variable=.true.,domain=domain_bottom)

      ! Link source-sink terms associated with layer-specific mass to the master module.
      call self%request_coupling(id_sms,'content_calculator_'//trim(name)//'/c_sms')

      ! Allow bulk coupling to a particulate organic matter module (rather than requiring the user to provide links per constituent and penetration depth)
      call self%request_coupling_to_model(id_c_int,self%id_Q,name)
      call self%request_coupling_to_model(id_pen_depth,self%id_Q,'pen_depth_'//trim(name))

      ! Create an alias in the master model for the layer-integrated density computed by layer_content_calculator.
      call self%add_horizontal_variable(name,trim(units)//'m^2','layer-integrated '//trim(long_name),domain=domain_bottom,act_as_state_variable=.true.)
      call self%request_coupling(name,'content_calculator_'//trim(name)//'/c')

      ! Register a dependency on the diagnostic that holds the layer-integrated density.
      ! This will be used to compute remineralization
      call self%register_dependency(id_local,trim(name)//'_local',trim(units)//'m^2','layer-integrated '//trim(long_name))
      call self%request_coupling(id_local,'content_calculator_'//trim(name)//'/c')
      if (self%remin/=0.0_rk) call self%register_state_dependency(id_remin_target,trim(name)//'_remin_target',trim(units)//'m^2','sink for remineralized '//trim(long_name))

      select case (name)
         case ('c'); call layer_content_calculator%add_to_aggregate_variable(standard_variables%total_carbon,layer_content_calculator%id_c,1.0_rk/CMass)
         case ('n'); call layer_content_calculator%add_to_aggregate_variable(standard_variables%total_nitrogen,layer_content_calculator%id_c)
         case ('p'); call layer_content_calculator%add_to_aggregate_variable(standard_variables%total_phosphorus,layer_content_calculator%id_c)
         case ('s'); call layer_content_calculator%add_to_aggregate_variable(standard_variables%total_silicate,layer_content_calculator%id_c)
      end select
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
         if (_VARIABLE_REGISTERED_(self%id_c_int)) call layer_process_constituent_changes(self,_ARGUMENTS_LOCAL_,d_sms,self%id_c_int,self%id_pen_depth_c,self%id_c_sms,self%id_c_local,self%id_c_remin_target)
         if (_VARIABLE_REGISTERED_(self%id_n_int)) call layer_process_constituent_changes(self,_ARGUMENTS_LOCAL_,d_sms,self%id_n_int,self%id_pen_depth_n,self%id_n_sms,self%id_n_local,self%id_n_remin_target)
         if (_VARIABLE_REGISTERED_(self%id_p_int)) call layer_process_constituent_changes(self,_ARGUMENTS_LOCAL_,d_sms,self%id_p_int,self%id_pen_depth_p,self%id_p_sms,self%id_p_local,self%id_p_remin_target)
         if (_VARIABLE_REGISTERED_(self%id_s_int)) call layer_process_constituent_changes(self,_ARGUMENTS_LOCAL_,d_sms,self%id_s_int,self%id_pen_depth_s,self%id_s_sms,self%id_s_local,self%id_s_remin_target)
      _HORIZONTAL_LOOP_END_
   end subroutine layer_do_bottom

   subroutine layer_process_constituent_changes(self,_ARGUMENTS_LOCAL_,d_sms,id_c_int,id_pen_depth,id_sms,id_local,id_remin_target)
      class (type_ersem_benthic_pom_layer), intent(in) :: self
      _DECLARE_ARGUMENTS_LOCAL_
      type (type_bottom_state_variable_id),intent(in) :: id_c_int,id_pen_depth,id_remin_target
      type (type_horizontal_dependency_id),intent(in) :: id_sms,id_local
      real(rk),                            intent(in) :: d_sms

      real(rk) :: c_int,d_pen,sms
      real(rk) :: c_int_local,sms_remin

      ! Retrieve depth-integrated mass, penetration depth, sinks-sources.
      _GET_HORIZONTAL_(id_c_int,c_int)
      _GET_HORIZONTAL_(id_pen_depth,d_pen)
      _GET_HORIZONTAL_(id_sms,sms)

      ! Convert sources-sinks in per second, as returned by FABM, to our own time unit.
      sms = sms*self%dt

      ! Add local remineralization
      if (self%remin/=0.0_rk) then
         _GET_HORIZONTAL_(id_c_int,c_int_local)
         sms_remin = self%remin*c_int_local
         _SET_BOTTOM_ODE_(id_remin_target,sms_remin)
         sms = sms - sms_remin
      end if

      ! Apply sinks-sources to depth-integrated mass, compute change in penetration depth.
      ! For the change in penentration depth, see its derivation in the comments at the top of the file,
      ! section "Impact of bioturbation".
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