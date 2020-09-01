#include "fabm_driver.h"

module ersem_benthic_column_dissolved_matter

   use fabm_types
   use fabm_builtin_models

   use ersem_shared

   implicit none

   private

   integer, parameter :: nlayers = 3

   type type_single_constituent
      type (type_bottom_state_variable_id) :: id_int             ! depth-integrated mass in the benthic column
      type (type_bottom_state_variable_id) :: id_int_deep        ! depth-integrated [negative] mass below the zero-concentration isocline
      type (type_state_variable_id)        :: id_pel             ! concentration in bottom-most pelagic cell
      type (type_horizontal_dependency_id) :: id_sms(nlayers)    ! sources-sinks specified for layer-integrated variables (pore water and adsorbed)
      type (type_horizontal_dependency_id) :: id_pw_sms(nlayers) ! sources-sinks specified for layer-integrated pore water concentration
      type (type_horizontal_diagnostic_variable_id) :: id_pbf    ! pelagic-benthic flux
      logical :: nonnegative = .true.                            ! Flag specifying whether constituent is non-negative
   end type

   ! Model for dissolved matter in sediment, using idealized equilibrium profiles
   ! to determine pelagic-benthic diffusive flux, and to determine equilibrium depths of
   ! all but the last layer (NB layer depth will be relaxed to equilibrium value).
   type,extends(type_base_model),public :: type_ersem_benthic_column_dissolved_matter
      type (type_horizontal_dependency_id) :: id_Dm(nlayers)    ! depth of bottom interface of individual layers (last is total column height)
      type (type_horizontal_dependency_id) :: id_poro           ! porosity (currently one single value is used across all layers)
      type (type_horizontal_dependency_id) :: id_cmix           ! pelagic-benthic transfer in bottom boundary layer (height of BBL divided by diffusivity within BBL)
      type (type_horizontal_dependency_id) :: id_diff(nlayers)  ! effective diffusivity within individual layers [includes bioirrigation contribution, if any]
      type (type_bottom_state_variable_id) :: id_layer          ! depth of bottom interface of own layer (where own concentration drops to zero)
      type (type_horizontal_diagnostic_variable_id) :: id_conc_eq(nlayers)  ! mean equilibrium pore water concentration in individual layers
      type (type_horizontal_diagnostic_variable_id) :: id_conc_tot(nlayers) ! mean pore water concentration in individual layers
      real(rk) :: ads(nlayers)
      real(rk) :: relax, minD
      integer :: last_layer
      logical :: correction
      type (type_single_constituent),allocatable :: constituents(:)
   contains
      procedure :: initialize => benthic_dissolved_matter_initialize
      procedure :: do_bottom  => benthic_dissolved_matter_do_bottom
   end type

   type type_single_constituent_estimates
      type (type_horizontal_diagnostic_variable_id) :: id_conc                        ! depth-averaged pore water concentration
      type (type_horizontal_dependency_id)          :: id_int                         ! depth-integrated total mass (mass/m2) in entire column [pore water + adsorbed]
      type (type_horizontal_diagnostic_variable_id) :: id_per_layer_total(nlayers)    ! depth-integrated total mass (mass/m2) per layer [pore water + adsorbed]
      type (type_horizontal_diagnostic_variable_id) :: id_per_layer_pw_total(nlayers) ! depth-integrated pore water mass (mass/m2) per layer
   end type

   ! Model that provides estimates of the depth-integrated mass within individual layers.
   ! This model is not meant to be instantiated by users, but is created automatically
   ! as child "per_layer" of the main model, type_ersem_benthic_column_dissolved_matter.
   ! Other models can access the layer-specific mass variables and provide them with
   ! source and sink terms.
   type,extends(type_base_model) :: type_dissolved_matter_per_layer
      type (type_horizontal_dependency_id) :: id_Dm(nlayers)
      type (type_horizontal_dependency_id) :: id_poro
      real(rk) :: ads(nlayers)
      integer :: last_layer
      type (type_single_constituent_estimates),allocatable :: constituents(:)
   contains
      procedure :: do_bottom => dissolved_matter_per_layer_do_bottom
   end type

contains

   subroutine benthic_dissolved_matter_initialize(self,configunit)
      class (type_ersem_benthic_column_dissolved_matter),intent(inout),target :: self
      integer,                                           intent(in)           :: configunit

      character(len=10) :: composition
      integer           :: ilayer, iconstituent
      character(len=16) :: index
      real(rk)          :: c0
      type (type_horizontal_standard_variable) :: standard_variable

      class (type_dissolved_matter_per_layer), pointer :: profile

      ! Set time unit to d-1. This implies that all rates (sink/source terms) are given in d-1.
      self%dt = 86400._rk

      ! Obtain parameter values
      call self%get_parameter(composition,'composition','','composition (any combination of c,n,p,s,o,a)')
      call self%get_parameter(c0,'c0','mg C/m^2','background carbon concentration',default=0.0_rk)
      call self%get_parameter(self%last_layer,'last_layer','','sediment layer where concentration drops to zero',default=nlayers)
      if (composition=='') call self%fatal_error('benthic_dissolved_matter_initialize','composition must include at least one chemical constituent')
      if (self%last_layer/=nlayers .and. len_trim(composition)>1) call self%fatal_error('benthic_dissolved_matter_initialize','last_layer cannot be set for solutes with more than one chemical constituent')
      if (self%last_layer/=nlayers) then
         call self%get_parameter(self%relax,'relax','1/d','rate of relaxation towards equilibrium concentration profile')
         call self%get_parameter(self%minD, 'minD','m',  'minimum depth of zero-concentration isocline')

         write (index,'(i0)') self%last_layer
         standard_variable%name = 'depth_of_bottom_interface_of_layer_'//trim(index)
         call self%register_state_dependency(self%id_layer, standard_variable)
      end if
      self%ads = 1.0_rk
      do ilayer=1,self%last_layer
         write (index,'(i0)') ilayer
         call self%get_parameter(self%ads(ilayer),'ads'//trim(index),'-','adsorption in layer '//trim(index)//' (total:dissolved)',default=1.0_rk)
      end do
      call self%get_parameter(self%correction,'correction','','move losses in oxygenic layer to deeper layers if pelagic concentration is limiting',default=.false.)

      ! Create model that computes concentrations per benthic layer.
      allocate(profile)
      call self%add_child(profile,'per_layer',configunit=-1)
      profile%ads = self%ads
      profile%last_layer = self%last_layer

      ! Create constituent-specific variables.
      allocate(self%constituents(len_trim(composition)))
      allocate(profile%constituents(len_trim(composition)))
      do iconstituent=1,len_trim(composition)
         select case (composition(iconstituent:iconstituent))
         case ('c')
            call initialize_constituent(self,self%constituents(iconstituent),profile,profile%constituents(iconstituent),'c','mmol/m^2','carbon',standard_variables%total_carbon,c0)
         case ('n')
            call initialize_constituent(self,self%constituents(iconstituent),profile,profile%constituents(iconstituent),'n','mmol/m^2','nitrogen',standard_variables%total_nitrogen,c0)
         case ('p')
            call initialize_constituent(self,self%constituents(iconstituent),profile,profile%constituents(iconstituent),'p','mmol/m^2','phosphorus',standard_variables%total_phosphorus,c0)
         case ('s')
            call initialize_constituent(self,self%constituents(iconstituent),profile,profile%constituents(iconstituent),'s','mmol/m^2','silicate',standard_variables%total_silicate,c0)
         case ('o')
            call initialize_constituent(self,self%constituents(iconstituent),profile,profile%constituents(iconstituent),'o','mmol/m^2','oxygen', nonnegative=legacy_ersem_compatibility)
         case ('a')
            call initialize_constituent(self,self%constituents(iconstituent),profile,profile%constituents(iconstituent),'a','mEq/m^2','alkalinity')
         case default
            call self%fatal_error('benthic_dissolved_matter_initialize','Invalid value for parameter "composition". Permitted: c,n,p,s,o,a.')
         end select
      end do

      ! Dependencies
      call self%register_dependency(self%id_cmix,pelagic_benthic_transfer_constant)
      call self%register_dependency(self%id_poro,sediment_porosity)
      call profile%register_dependency(profile%id_poro,sediment_porosity)
      do ilayer=1,nlayers
         write (index,'(i0)') ilayer
         standard_variable%name = 'diffusivity_in_sediment_layer_'//trim(index)
         call self%register_dependency(self%id_diff(ilayer),    standard_variable)
         standard_variable%name = 'depth_of_bottom_interface_of_layer_'//trim(index)
         call self%register_dependency(self%id_Dm(ilayer),      standard_variable)
         call profile%register_dependency(profile%id_Dm(ilayer),standard_variable)
      end do
      call self%request_coupling   (self%id_Dm(nlayers),   depth_of_sediment_column)
      call profile%request_coupling(profile%id_Dm(nlayers),depth_of_sediment_column)

   end subroutine benthic_dissolved_matter_initialize

   subroutine initialize_constituent(self,info,profile,profile_info,name,units,long_name,aggregate_target,background_value,nonnegative)
      class (type_ersem_benthic_column_dissolved_matter),intent(inout),target :: self
      type (type_single_constituent),                    intent(inout),target :: info
      class (type_dissolved_matter_per_layer),           intent(inout),target :: profile
      type (type_single_constituent_estimates),          intent(inout),target :: profile_info
      character(len=*),                                  intent(in)           :: name,units,long_name
      type (type_bulk_standard_variable),optional,       intent(in)           :: aggregate_target
      real(rk),optional,                                 intent(in)           :: background_value
      logical,optional,                                  intent(in)           :: nonnegative

      integer           :: ilayer
      character(len=16) :: index

      if (present(nonnegative)) info%nonnegative = nonnegative

      ! State variable for depth-integrated matter [adsorbed + in pore water]
      call self%register_state_variable(info%id_int,name,units,long_name,background_value=background_value)
      if (.not.info%nonnegative) call self%register_state_variable(info%id_int_deep,name//'_deep',units,long_name//' below zero isocline')

      ! Contribution of state variable to aggregate quantity (if any).
      if (present(aggregate_target)) call self%add_to_aggregate_variable(aggregate_target,info%id_int)

      ! Dependency on lowermost pelagic concentration.
      call self%register_state_dependency(info%id_pel,trim(name)//'_pel','mmol/m^3','pelagic '//trim(long_name))

      ! Diagnostic for pelagic-benthic flux.
      call self%register_diagnostic_variable(info%id_pbf,trim(name)//'_pb_flux','mmol/m^2/d','flux of '//trim(long_name)//' from benthos to pelagic',source=source_do_bottom)

      ! Register new constituent with child model that computes mass per benthic layer.
      call profile%register_dependency(profile_info%id_int,trim(name)//'_int',units,'depth-integrated '//trim(long_name))
      call profile%request_coupling(profile_info%id_int,name)
      do ilayer=1,nlayers
         write (index,'(i0)') ilayer

         ! Register diagnostics for total [adsorbed + pore water] layer integral, and pore water integral.
         ! These act as state variables, in order to allow other modules to provide source terms for them.
         call profile%register_diagnostic_variable(profile_info%id_per_layer_total(ilayer),trim(name)//trim(index),units,'total '//trim(long_name)//' in layer '//trim(index)//' (absorbed + dissolved)',act_as_state_variable=.true.,domain=domain_bottom,output=output_none,source=source_do_bottom)
         call profile%register_diagnostic_variable(profile_info%id_per_layer_pw_total(ilayer),trim(name)//trim(index)//'_pw',units,'total '//trim(long_name)//' in pore water of layer '//trim(index),act_as_state_variable=.true.,domain=domain_bottom,output=output_none,source=source_do_bottom)
         if (ilayer <= self%last_layer .or. info%nonnegative) then
            ! Redirect sources-sinks for this layer to the depth-integrated value (or to the integral up to the zero isocline, if present)
            call copy_horizontal_fluxes(profile, profile_info%id_per_layer_total(ilayer), name)
            call copy_horizontal_fluxes(profile, profile_info%id_per_layer_pw_total(ilayer), name)
         else
            ! Redirect sources-sinks for this layer to the depth-integrated value below the zero isocline
            call copy_horizontal_fluxes(profile, profile_info%id_per_layer_total(ilayer), name//'_deep')
            call copy_horizontal_fluxes(profile, profile_info%id_per_layer_pw_total(ilayer), name//'_deep')
         end if

         ! Make sure that sources-sinks of layer-integrated mass are counted in conservation checks on source/sink basis (e.g., with check_conservation).
         if (present(aggregate_target)) then
            call profile%add_to_aggregate_variable(aggregate_target,profile_info%id_per_layer_total(ilayer))
            call profile%add_to_aggregate_variable(aggregate_target,profile_info%id_per_layer_pw_total(ilayer))
         end if

         ! Collect sources for depth-integrated total [pore water + adsorbed] matter
         call self%register_dependency(info%id_sms(ilayer),trim(name)//'_sms_l'//trim(index),trim(units)//'/s','sources-sinks of total '//trim(long_name)//' in layer '//trim(index))
         call self%request_coupling(info%id_sms(ilayer),'per_layer/'//trim(name)//trim(index)//'_sms_tot')

         ! Collect sources for depth-integrated pore water concentration
         call self%register_dependency(info%id_pw_sms(ilayer),trim(name)//'_pw_sms_l'//trim(index),trim(units)//'/s','sources-sinks of dissolved '//trim(long_name)//' in layer '//trim(index))
         call self%request_coupling(info%id_pw_sms(ilayer),'per_layer/'//trim(name)//trim(index)//'_pw_sms_tot')

         ! Register diagnostics for layer-specific pore water concentrations (equilibrium and actual)
         if (self%last_layer==nlayers) then
            call self%register_diagnostic_variable(self%id_conc_eq(ilayer),trim(name)//trim(index)//'_eq_conc','mmol/m^3','equilibrium pore water concentration of '//trim(long_name)// ' in layer '//trim(index),domain=domain_bottom,source=source_do_bottom)
            call self%register_diagnostic_variable(self%id_conc_tot(ilayer),trim(name)//trim(index)//'_conc','mmol/m^3','pore water concentration of '//trim(long_name)// ' in layer '//trim(index),domain=domain_bottom,source=source_do_bottom)
         end if
      end do

      ! Create diagnostic for depth-averaged pore water concentration (currently only for solutes that span entire column)
      if (profile%last_layer==nlayers) call profile%register_diagnostic_variable(profile_info%id_conc,trim(name)//'_conc',trim(units)//'/m','depth-averaged pore water concentration of '//trim(long_name),domain=domain_bottom,source=source_do_bottom)

   end subroutine initialize_constituent

   subroutine benthic_dissolved_matter_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_benthic_column_dissolved_matter),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      integer :: iconstituent

      do iconstituent=1,size(self%constituents)
         call process_constituent(self,_ARGUMENTS_DO_BOTTOM_,self%constituents(iconstituent))
      end do
   end subroutine benthic_dissolved_matter_do_bottom

   subroutine process_constituent(self,_ARGUMENTS_DO_BOTTOM_,info)
      class (type_ersem_benthic_column_dissolved_matter),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
      type (type_single_constituent),                    intent(in) :: info

      integer  :: ilayer
      real(rk) :: c_pel,c_top,c_int,c_int_deep
      real(rk) :: sms_per_layer(nlayers),pw_sms_per_layer(nlayers),sms
      real(rk) :: Dm(nlayers)
      real(rk) :: c_bot,c_int_per_layer_eq(nlayers),H_eq,d_top
      real(rk) :: c_int_eq
      real(rk) :: norm_res_int,P_res_int
      real(rk) :: smscorr
      real(rk) :: diff(nlayers),poro,cmix
      real(rk) :: residual_per_layer(nlayers)

      _HORIZONTAL_LOOP_BEGIN_

      ! Retrieve physical properties of sediment column:
      ! porosity, pelagic-bentic transfer coefficient, per-layer diffusivities, per-layer depth of the bottom interface.
      _GET_HORIZONTAL_(self%id_poro,poro)
      _GET_HORIZONTAL_(self%id_cmix,cmix)
      do ilayer=1,nlayers
         _GET_HORIZONTAL_(self%id_diff(ilayer),diff(ilayer))
         _GET_HORIZONTAL_(self%id_Dm(ilayer),Dm(ilayer))
      end do

      ! Retrieve column-integrated mass, lowermost pelagic concentration, and layer-specific depth-integrated sources.
      _GET_HORIZONTAL_(info%id_int,c_int)
      _GET_(info%id_pel,c_pel)
      do ilayer=1,nlayers
         _GET_HORIZONTAL_(info%id_sms(ilayer),sms_per_layer(ilayer))
         _GET_HORIZONTAL_(info%id_pw_sms(ilayer),pw_sms_per_layer(ilayer))
      end do

      ! Combine depth-integrated sources applied to total and pore-water-only matter.
      ! Also convert time units: FABM source terms are given in s-1, and we need d-1.
      sms_per_layer = (sms_per_layer+pw_sms_per_layer)*86400

      ! Column-integrated sources minus sinks (# m-2 d-1)
      sms = sum(sms_per_layer)

      ! Estimate steady state concentration at bed from current pelagic concentration
      ! (typically at centre of lowermost pelagic layer), and column-integrated production/destruction term.
      !
      ! This assumes the lowermost pelagic layer includes a diffusion barrier in the form of a bottom boundary layer (BBL),
      ! characterized by height H, [low] diffusivity D, and no production or destruction of the tracer.
      ! At equilibrium, the vertical flux at any point in the BBL must equal the pelagic-benthic flux, which
      ! in turn equals the production within the benthic column, sms. That implies the equilibrium concentration profile within
      ! the BBL is linear with a slope equal to sms/D. The change between top and bottom of the BBL thus equals sms/D*H.
      ! If the structure (height H, diffusivity D) of the BBL are constant in time and space, this can be rewritten as sms*cmix,
      ! with transfer coefficient cmix = H/D in units d/m. If the lowermost layer excluding BBL is well-mixed, the concentration
      ! at the top of the BBL is equal to the lowermost pelagic concentration c_pel, which is typically vertically positioned at the
      ! centre of the lowermost grid cell. Thus, the concentration at the top of the sediment bed equals:
      !
      !    c_top = c_pel + cmix*sms
      !
      ! For positive definite tracers, this expression presents problems if the sediment column destroys tracer at a fast pace:
      ! if sms<-c_pel/cmix, c_top would become negative. To counter that problem, we apply the Patankar trick
      ! (Patankar 1980; Burchard et al. 2003 ApNumMat Eq 16), which guarantees a positive c_top while preserving first order
      ! accuracy of the (originally linear) relation. JB 24/06/2015, 10/11/2015
      if (sms>=0._rk .or. .not. info%nonnegative) then
         ! Sediment column produces tracer or tracer does not need to be positive definite. Use original linear relation.
         c_top = c_pel + cmix*sms
      else
         ! Sediment column destroys tracer while tracer needs to be positive definite. Apply Patankar trick.
         c_top = c_pel*c_pel/(c_pel - cmix*sms)
      end if
      d_top = 0

      if (self%last_layer/=nlayers) then
         ! This constituent drops to zero in layer number self%last_layer
         ! (and the depth where it equals zero acts as this last layer's lower boundary).

         c_int_per_layer_eq = 0

         ! First compute equilibrium concentration profiles in all layers but the last in order to obtain equilibrium values for
         ! the concentration at their bottom interface, c_bot, and their depth-integrated concentration, c_int_per_layer_eq(ilayer)
         ! NB depth-integrated concentration refers to the integral of the pore water concentration profile from layer top to layer bottom,
         ! without accounting for porosity or adsorption.
         do ilayer=1,self%last_layer-1
            call compute_equilibrium_profile(diff(ilayer),c_top,sms_per_layer(ilayer),sum(sms_per_layer(ilayer+1:)),Dm(ilayer)-d_top,c_bot,c_int_per_layer_eq(ilayer))
            d_top = Dm(ilayer)
            c_top = c_bot
         end do

         ! Last layer: the pore water concentration drops to zero.
         ! Compute equilibrium layer height H_eq and depth-integrated concentration c_int_per_layer_eq(self%last_layer).
         ! The concentration at the bottom interface is zero by definition.
         call compute_final_equilibrium_profile(diff(self%last_layer),c_top,sms_per_layer(self%last_layer),sum(sms_per_layer(self%last_layer+1:)),Dm(nlayers)-d_top,H_eq,c_int_per_layer_eq(self%last_layer))

         d_top = d_top + max(self%minD, H_eq)
         if (H_eq <= self%minD) then
            ! Layer too thin - treat it as completely collapsed (next layer starts with orginal surface concentration c_top)
            c_int_per_layer_eq(self%last_layer) = 0
         else
            ! Layer present - next layer starts with 0 concentration at its surface
            c_top = 0
         end if

         ! Relax depth-integrated mass c_int towards its equilibrium value (sum of depth-integrated equilibrium values of all layers)
         _SET_BOTTOM_ODE_(info%id_int, (poro*sum(self%ads(:self%last_layer)*c_int_per_layer_eq(:self%last_layer))-c_int)/self%relax - sum(sms_per_layer(:self%last_layer)))

         ! Relax the depth of the bottom interface of the last layer towards equilibrium value
         _SET_BOTTOM_ODE_(self%id_layer, (d_top - Dm(self%last_layer)) / self%relax)

         if (.not.info%nonnegative) then
            ! Deeper source terms are allowed to be non-zero [typically negative].
            ! In that case, the concentration gradient at the bottom of the last layer will be non-zero too,
            ! and the deeper layers will contain [typically negative] matter. Integrate this and add it to the column integral.
            do ilayer=self%last_layer+1,nlayers
               call compute_equilibrium_profile(diff(ilayer),c_top,sms_per_layer(ilayer),sum(sms_per_layer(ilayer+1:)),max(Dm(ilayer)-d_top, 0._rk),c_bot,c_int_per_layer_eq(ilayer))
               d_top = Dm(ilayer)
               c_top = c_bot
            end do

            ! Relax depth-integrated mass within deeper layers towards equilibrium value.
            _GET_HORIZONTAL_(info%id_int_deep,c_int_deep)
            _SET_BOTTOM_ODE_(info%id_int_deep,(poro*sum(self%ads(self%last_layer+1:)*c_int_per_layer_eq(self%last_layer+1:))-c_int_deep)/self%relax - sum(sms_per_layer(self%last_layer+1:)))
         else
            c_int_deep = 0
         end if

         ! From layer-specific depth integrals of the pore water concentration profile to column-integrated mass [adsorbed + in pore water]
         c_int_eq = poro*sum(self%ads*c_int_per_layer_eq)

         ! Net change in column-integrated mass must equal column-integrated production - surface exchange.
         ! Thus, surface exchange = column-integrated production - net change (net change = relaxation)
         _SET_BOTTOM_EXCHANGE_(info%id_pel,sms-(c_int_eq-(c_int+c_int_deep))/self%relax)
         _SET_HORIZONTAL_DIAGNOSTIC_(info%id_pbf,sms-(c_int_eq-(c_int+c_int_deep))/self%relax)
      else
         ! Apply a "technical correction" in case flux from the oxygenated
         ! layer is negative by scaling flux using pelagic concentration and
         ! redistributing it between benthic layers. Note, that in current
         ! version pelagic concentration at sediment interface is used
         ! instead of mean pelagic concentration of the older code.
         ! This "technical correction" is a hack that was initially applied to
         ! ammonium as an attempt to preserve positive concentrations and should
         ! be replaced in future.
         if (self%correction) then
            smscorr = sms_per_layer(1)
            if (smscorr .lt. 0._rk) then
               sms_per_layer(1) = smscorr*c_top/(c_top+0.5_rk)
               sms_per_layer(2) = sms_per_layer(2) + (smscorr-sms_per_layer(1))
            end if
         end if

         ! Compute equilibrium concentration profiles in all layers in order to obtain equilibrium values for
         ! the concentration at their bottom interface, c_bot, and their depth-integrated concentration, c_int_per_layer_eq(ilayer)
         do ilayer=1,nlayers
            call compute_equilibrium_profile(diff(ilayer),c_top,sms_per_layer(ilayer),sum(sms_per_layer(ilayer+1:)),Dm(ilayer)-d_top,c_bot,c_int_per_layer_eq(ilayer))
          _SET_HORIZONTAL_DIAGNOSTIC_(self%id_conc_eq(ilayer),c_int_per_layer_eq(ilayer)/(Dm(ilayer)-d_top))
            d_top = Dm(ilayer)
            c_top = c_bot
         end do

         ! From layer-specific depth integrals of the pore water concentration profile to column-integrated mass [adsorbed + in pore water]
         c_int_eq = poro*sum(self%ads*c_int_per_layer_eq)

         ! The equilibrium depth-integrated mass c_int_eq usually differs from the current depth-integrated mass.
         ! We can view the actual [unknown] pore water concentration profile as the sum of the equilibrium profile
         ! and a residual profile. The latter has a vertical integral equal to the difference between actual mass and equilibrium mass.
         ! Our aim is to get at the rate at which the mass is exchanged over the benthic-pelagic interface. This is equal to the product of
         ! the diffusivity and gradient in the residual mass at this interface. We do not know this gradient since we do not
         ! know the shape of the residual profile. Let's make some simple assumptions to infer this.
         ! Constraints: diffusion of the residual across bottom of benthic column must be zero (i.e., zero gradient), and at the surface of the benthic
         ! column the concentration of the residual must equal zero (i.e., equilibrium holds at the very surface of the column).
         ! Since we do not know anything about the processes responsible for the residual, let's assume their contribution
         ! in the past was a constant production or destruction per unit sediment volume thoughout the entire column.
         ! That is, production (#/m^2/d) in the three layers was P_int*d1/d3, P_int(d2-d1)/d3, P_int(d3-d2)/d3.
         ! If we would know P_int, we could supply those rates along with zero surface concentration to "compute_equilibrium_profile"
         ! to derive the residual profile. By checking the equations in compute_equilibrium_profile, we can verify that the resulting bottom concentration
         ! and layer integral are both proportional to P_int. Thus, can can simply supply d1, d2-d1, d3-d2 to
         ! "compute_equilibrium_profile", and find the additional scale factor P_int/d3 by demanding that the sum of layer integrals is
         ! equal to the known residual mass. That is, P_int/d3 equals the ratio of residual mass to the sum of normalized layer integrals
         ! computed for layer production terms d1, d2-d1, d3-d2. As we are assuming the residual profile was previously an equilibrium
         ! profile, the necessary depth-integrated production rate P_int must equal the exchange across the surface, i.e., diffusivity*gradient.
         ! Thus, we can now simply add the P_int as a additional surface exchange term, accounting for the move towards equilibrium.
         d_top = 0
         c_top = 0
         do ilayer=1,nlayers
            call compute_equilibrium_profile(diff(ilayer),c_top,Dm(ilayer)-d_top,Dm(nlayers)-Dm(ilayer),Dm(ilayer)-d_top,c_bot,residual_per_layer(ilayer))
            d_top = Dm(ilayer)
            c_top = c_bot
         end do
         norm_res_int = poro*sum(self%ads*residual_per_layer)
         P_res_int = (c_int-c_int_eq)/norm_res_int*Dm(nlayers)
         _SET_BOTTOM_EXCHANGE_(info%id_pel,sms+P_res_int) ! Equilibrium flux = depth-integrated production sms + residual flux P_res_int
         _SET_HORIZONTAL_DIAGNOSTIC_(info%id_pbf,sms+P_res_int)
         _SET_BOTTOM_ODE_(info%id_int,-P_res_int-sms)

         ! Save final estimates of mean pore water concentration per layer.
         d_top = 0
         do ilayer=1,nlayers
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_conc_tot(ilayer),(P_res_int/Dm(nlayers)*residual_per_layer(ilayer)+c_int_per_layer_eq(ilayer))/(Dm(ilayer)-d_top))
            d_top = Dm(ilayer)
         end do

      end if

      _HORIZONTAL_LOOP_END_
   end subroutine process_constituent

   subroutine compute_equilibrium_profile(sigma,C0,P,P_deep,D,C_bot,C_int)
      real(rk),intent(in)  :: sigma,c0,P,P_deep,D
      real(rk),intent(out) :: C_bot,C_int
      real(rk) :: a_D2,b_D,c
      ! ----------------------------------------------------------------------------------------------------
      ! Determine equilibrium concentration profile in pore water from:
      ! - sigma:    diffusivity (m2/d)
      ! - C0:       concentration at layer surface (#/m3)
      ! - P:        layer-integrated source-sink terms of current layer (#/m2/d)
      ! - P_deep:   depth-integrated source-sink terms of deeper layers (#/m2/d) [equal to -bottom flux at equilibrium]
      ! - D:        layer thickness (m)
      ! Task: find concentration at bottom interface from prescribed layer depth.
      ! Returns:
      ! - C_bot:    equilibrium concentration at bottom interface (#/m3)
      ! - C_int:    depth-integrated concentration in current layer (#/m2)
      ! ----------------------------------------------------------------------------------------------------
      ! Governing equation: $\partial C/\partial t = sigma \partial^2 C/\partial z^2 + sms$
      ! Assumption: diffusivity $\sigma$ and sources-minus-sinks $sms$ are independent of $z$ within the layer.
      ! Thus, $sms$ can be written as layer integrated source-sink terms, divided by layer height: $sms = P/D$
      ! At equilibrium: $\sigma \partial^2 C/\partial z^2 + P/D = 0$
      ! Thus, \partial^2 C/\partial z^2 = -P/D/\sigma
      ! Solution is a quadratic equation: $C(z) = a z^2 + b z + c$
      ! ----------------------------------------------------------------------------------------------------
      ! First we determine the coefficients a,b,c of the quadratic equation:
      !
      ! a)
      ! From second derivative: $\partial^2 C/\partial z^2 = 2a = -P/D/\sigma$.
      ! Thus, $a = -P/D/\sigma/2$.
      !
      ! b)
      ! Inward flux at top interface must balance depth-integrated sinks-sources in deeper layers:
      ! $-P-P_deep$. For consistency, this flux must equal that produced by diffusion at the boundary,
      ! i.e., $\sigma*\partial C/\partial z$. Note: gradient must be positive when deeper layers are a
      ! source (i.e., P+P_deep>0). Thus:
      !    \sigma*\partial C/\partial z(0) = \sigma*b = P+P_deep -> b = (P+P_deep)/\sigma
      !
      ! c)
      ! Concentration at top interface $C(0) = c = C0$ is known.
      !
      ! Thus the pore water concentration is prescribed by the quadratic equation $C(z) = a z^2 + b z + c$
      ! with constants:
      !   $a = -P/D/\sigma/2$
      !   $b = (P+P_deep)/\sigma$
      !   $c = C0$
      ! ----------------------------------------------------------------------------------------------------
      ! Problem 1: find the concentration at the bottom interface of the layer (units: #/m3):
      !
      ! The concentration at the bottom equals the value of the parabola at depth D:
      !    C(D) = a D^2 + b D + c
      ! ----------------------------------------------------------------------------------------------------
      ! Problem 2: find the concentration integrated over the layer (units: #/m2)
      !
      ! Depth-integrated layer contents is found by integrating the parabola between 0 and D:
      !    \int_0^D C(z) dz = \int_0^D a * z^2 + b z + c dz
      !                     = [a/3 z^3 + b/2 z^2 + c z]_0^D
      !                     = a/3 D^3 + b/2 D^2 + c D
      ! ----------------------------------------------------------------------------------------------------
      ! Preventing division by zero (due to D=0):
      ! As $a$ always occurs multiplied with $D^2$, and $b$ always multiplied with $D$ (when solving
      ! problems 1 and 2), we define the combined constants $a_D2 = a D^2$ and $b_D = b D$. This avoids
      ! division by 0 when computing $a$ while $D$ tends to zero.
      ! ----------------------------------------------------------------------------------------------------
      a_D2 = -P/sigma/2*D
      b_D = (P+P_deep)/sigma*D
      c = C0
      C_bot = a_D2 + b_D + c
      C_int = (a_D2/3 + b_D/2 + c)*D
   end subroutine compute_equilibrium_profile

   subroutine compute_final_equilibrium_profile(sigma,c0,P,P_deep,Dmax,D,c_int)
      real(rk),intent(in)  :: sigma,c0,P,P_deep,Dmax
      real(rk),intent(out) :: D,c_int
      real(rk) :: a_D2,b_D,c,c_bot
      ! ----------------------------------------------------------------------------------------------------
      ! Determine layer depth and concentration profile in pore water at equilibrium from:
      ! - sigma:    diffusivity (m2/d)
      ! - C0:       concentration at layer surface (#/m3)
      ! - P:        layer-integrated source-sink terms (#/m2/d)
      ! - P_deep:   depth-integrated source-sink terms of deeper layers (#/m2/d) [equal to -bottom flux at equilibrium]
      ! Task: find layer depth from prescribed concentration at bottom interface (zero).
      ! Returns:
      ! - D:        layer thickness (m)
      ! - C_int:    depth integral of concentration (#/m2)
      ! ----------------------------------------------------------------------------------------------------
      ! Governing equation: $\partial C/\partial t = sigma \partial^2 C/\partial z^2 + sms$
      ! Assumption: diffusivity $\sigma$ and sources-minus-sinks $sms$ are independent of $z$ within the layer.
      ! Thus, $sms$ can be written as layer integrated source-sink terms, divided by layer height: $sms = P/D$
      ! At equilibrium: $\sigma \partial^2 C/\partial z^2 + P/D = 0$
      ! Thus, \partial^2 C/\partial z^2 = -P/D/\sigma
      ! Solution is a quadratic equation: $C(z) = a z^2 + b z + c$
      ! ----------------------------------------------------------------------------------------------------
      ! First we determine the coefficients a,b,c of the quadratic equation:
      !
      ! a)
      ! From second derivative: $\partial^2 C/\partial z^2 = 2a = -P/D/\sigma$.
      ! Thus, $a = -P/D/\sigma/2$.
      !
      ! b)
      ! Inward flux at top interface must balance depth-integrated sinks-sources in deeper layers:
      ! $-P-P_deep$. For consistency, this flux must equal that produced by diffusion at the boundary,
      ! i.e., $\sigma*\partial C/\partial z$. Note: gradient must be positive when deeper layers are a
      ! source (i.e., P+P_deep>0). Thus:
      !    \sigma*\partial C/\partial z(0) = \sigma*b = P+P_deep -> b = (P+P_deep)/\sigma
      !
      ! c)
      ! Concentration at top interface $C(0) = c = C0$ is known.
      !
      ! Thus the pore water concentration is prescribed by the quadratic equation $C(z) = a z^2 + b z + c$
      ! with constants:
      !   $a = -P/D/\sigma/2$
      !   $b = (P+P_deep)/\sigma$
      !   $c = C0$
      ! ----------------------------------------------------------------------------------------------------
      ! Problem 1: find layer depth $D$, defined as the depth at which the pore water concentration is zero:
      !
      !   a D^2 + b D + c = 0
      !
      ! Inserting coefficients $a,b,c$:
      !
      !   -P/D/\sigma/2 D^2 + (P + P_deep)/\sigma D + C0 = (P/2 + P_deep)/\sigma D + C0 = 0
      !
      ! Isolating $D$:
      !
      !   D = -C0 \sigma/(P/2 + P_deep) = -2 C0 \sigma/(P + 2 P_deep)
      !
      ! Note that this returns a positive depth D only if P+2 P_deep<0!
      ! ----------------------------------------------------------------------------------------------------
      ! Problem 2: find the concentration integrated over the layer (units: #/m2)
      !
      ! Depth-integrated layer contents is found by integrating the parabola between 0 and D:
      !    \int_0^D C(z) dz = \int_0^D a * z^2 + b z + c dz
      !                     = [a/3 z^3 + b/2 z^2 + c z]_0^D
      !                     = a/3 D^3 + b/2 D^2 + c D
      ! ----------------------------------------------------------------------------------------------------
      ! Preventing division by zero (due to D=0):
      ! As $a$ always occurs multiplied with $D^2$, and $b$ always multiplied with $D$,
      ! we define the combined constants $a_D2 = a D^2$ and $b_D = b D$. This avoids division by 0 when
      ! computing $a$ while $D$ tends to zero.
      ! ----------------------------------------------------------------------------------------------------

      if (Dmax*(P+2*P_deep)>=-2*sigma*c0) then   ! Dmax<-2 sigma c0/(P+2*P_deep), rearranged for P+2*P_deep<0 close to 0. Result also picks up P+2*P_deep>0.
         ! Destruction within layer is too low (or layer experiences net production, i.e., P>0).
         ! If we would impose zero concentration at the layer bottom, the layer would extend beyond maximum depth.
         ! Fix depth at maximum depth and use parabola with non-zero concentration at bottom interface.
         D = Dmax
         call compute_equilibrium_profile(sigma,c0,P,P_deep,D,c_bot,c_int)
      else
         D = -2*sigma*C0/(P+2*P_deep)
         a_D2 = -P/sigma/2*D
         b_D = (P+P_deep)/sigma*D
         c = C0
         C_int = (a_D2/3 + b_D/2 + c)*D
      end if
   end subroutine compute_final_equilibrium_profile

   subroutine dissolved_matter_per_layer_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_dissolved_matter_per_layer),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      integer  :: iconstituent
      integer  :: ilayer
      real(rk) :: c_int
      real(rk) :: pwconc
      real(rk) :: d(nlayers),poro

      _HORIZONTAL_LOOP_BEGIN_

         ! Layer depths (compute from depth of bottom interfaces) and porosity.
         do ilayer=1,nlayers
            _GET_HORIZONTAL_(self%id_Dm(ilayer),d(ilayer))
         end do
         d(2:nlayers) = d(2:nlayers) - d(1:nlayers-1)
         _GET_HORIZONTAL_(self%id_poro,poro)

         do iconstituent=1,size(self%constituents)
            ! Get depth-integrated concentration.
            _GET_HORIZONTAL_(self%constituents(iconstituent)%id_int,c_int)

            if (self%last_layer/=nlayers) then
               ! Vertically homogeneous in top layers, quadratically decreasing in last layer (zero concentration at bottom interface)

               ! Compute pore water concentration in upper layers (above self%last_layer), in matter/m3
               ! Note that concentration in last layer is 1/3 of that of the top layers,
               ! due to the fact that it has a quadratic profile decreasing from the top concentration to zero.
               pwconc = c_int/(poro*(sum(self%ads(:self%last_layer-1)*d(:self%last_layer-1)) + self%ads(self%last_layer)*d(self%last_layer)/3))

               ! Top layers: homogeneous pore water concentration.
               do ilayer=1,self%last_layer-1
                  ! Total depth-integrated layer contents [pore water + adsorbed]
                  _SET_HORIZONTAL_DIAGNOSTIC_(self%constituents(iconstituent)%id_per_layer_total(ilayer),poro*d(ilayer)*pwconc*self%ads(ilayer))

                  ! Pore water contents: layer contents divided by adsorption [total:dissolved]
                  _SET_HORIZONTAL_DIAGNOSTIC_(self%constituents(iconstituent)%id_per_layer_pw_total(ilayer),poro*d(ilayer)*pwconc)
               end do

               ! Last layer: 1/3 of top concentration.
               _SET_HORIZONTAL_DIAGNOSTIC_(self%constituents(iconstituent)%id_per_layer_total   (self%last_layer),poro*d(self%last_layer)/3*pwconc*self%ads(self%last_layer))
               _SET_HORIZONTAL_DIAGNOSTIC_(self%constituents(iconstituent)%id_per_layer_pw_total(self%last_layer),poro*d(self%last_layer)/3*pwconc)

               ! Deeper layers: pore water concentration is zero.
               do ilayer=self%last_layer+1,nlayers
                  _SET_HORIZONTAL_DIAGNOSTIC_(self%constituents(iconstituent)%id_per_layer_total   (ilayer),0.0_rk)
                  _SET_HORIZONTAL_DIAGNOSTIC_(self%constituents(iconstituent)%id_per_layer_pw_total(ilayer),0.0_rk)
               end do
            else
               ! Vertically homogeneous in all layers.

               ! Compute pore water concentration in matter/m3.
               pwconc = c_int/sum(poro*self%ads*d)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%constituents(iconstituent)%id_conc,pwconc)

               do ilayer=1,nlayers
                  ! Total depth-integrated layer contents [pore water + adsorbed]
                  _SET_HORIZONTAL_DIAGNOSTIC_(self%constituents(iconstituent)%id_per_layer_total(ilayer),poro*self%ads(ilayer)*d(ilayer)*pwconc)

                  ! Pore water contents: layer contents divided by adsorption [total:dissolved]
                  _SET_HORIZONTAL_DIAGNOSTIC_(self%constituents(iconstituent)%id_per_layer_pw_total(ilayer),poro*d(ilayer)*pwconc)
               end do
            end if
         end do

      _HORIZONTAL_LOOP_END_

   end subroutine dissolved_matter_per_layer_do_bottom

end module
