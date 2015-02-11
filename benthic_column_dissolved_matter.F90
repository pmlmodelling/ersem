#include "fabm_driver.h"

module ersem_benthic_column_dissolved_matter

   use fabm_types

   use ersem_shared

   implicit none

   private

   ! Model for dissolved matter in sediment, using idealized equilibrium profiles
   ! to determine pelagic-benthic diffusive flux, and to determine steady state depth of
   ! first and second layer (layer depth will be relaxed to steady state value).
   type,extends(type_base_model),public :: type_ersem_benthic_column_dissolved_matter
      type (type_bottom_state_variable_id) :: id_tot
      type (type_state_variable_id)        :: id_pel
      type (type_bottom_state_variable_id) :: id_D1m,id_D2m
      type (type_horizontal_dependency_id) :: id_sms(3),id_pw_sms(3),id_Dtot,id_poro,id_cmix,id_diff(3)
      type (type_horizontal_diagnostic_variable_id) :: id_pbf
      real(rk) :: ads(3)
      real(rk) :: relax, minD
      integer :: last_layer
   contains
      procedure :: initialize => benthic_dissolved_matter_initialize
      procedure :: do_bottom  => benthic_dissolved_matter_do_bottom
   end type

   ! Model for dissolved matter for a single sediment layer
   type,extends(type_base_model) :: type_dissolved_matter_per_layer
      type (type_horizontal_dependency_id)          :: id_tot
      type (type_horizontal_dependency_id)          :: id_D1m,id_D2m,id_Dtot,id_poro
      type (type_horizontal_diagnostic_variable_id) :: id_layers(3)    ! mass density (mass/m2) within the layer
      type (type_horizontal_diagnostic_variable_id) :: id_layers_pw(3) ! mass density (mass/m2) within layer pore water
      real(rk) :: ads(3)
      integer :: last_layer
   contains
      procedure :: do_bottom => dissolved_matter_per_layer_do_bottom
   end type

contains   

   subroutine benthic_dissolved_matter_initialize(self,configunit)
      class (type_ersem_benthic_column_dissolved_matter),intent(inout),target :: self
      integer,                                           intent(in)           :: configunit

      class (type_dissolved_matter_per_layer), pointer :: profile
      character(len=10) :: composition
      character(len=attribute_length) :: long_name

      ! Set time unit to d-1. This implies that all rates (sink/source terms) are given in d-1.
      self%dt = 86400._rk

      ! Create state variable for depth-integrated concentration.
      call self%get_parameter(composition,'composition','','composition')
      select case (composition)
         case ('c'); long_name = 'carbon'
         case ('n'); long_name = 'nitrogen'
         case ('p'); long_name = 'phosphorus'
         case ('s'); long_name = 'silicate'
         case ('o'); long_name = 'oxygen'
         case default
            call self%fatal_error('benthic_dissolved_matter_initialize','Invalid value for parameter "composition". Permitted: c,n,p,s,o.')
      end select
      call self%register_state_variable(self%id_tot,trim(composition),'mmol/m^2',long_name)
      select case (composition)
         case ('c'); call self%add_to_aggregate_variable(standard_variables%total_carbon,    self%id_tot)
         case ('n'); call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_tot)
         case ('p'); call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_tot)
         case ('s'); call self%add_to_aggregate_variable(standard_variables%total_silicate,  self%id_tot)
      end select

      ! Obtain parameter values
      call self%get_parameter(self%last_layer,  'last_layer','','deepest sediment layer',default=3)
      self%ads = 1.0_rk
      call self%get_parameter(self%ads(1),'ads1','-','adsorption in oxygenated layer (total:dissolved)',default=1.0_rk)
      if (self%last_layer>1) call self%get_parameter(self%ads(2),'ads2','-','adsorption in oxidized layer (total:dissolved)',default=1.0_rk)
      if (self%last_layer>2) call self%get_parameter(self%ads(3),'ads3','-','adsorption in anoxic layer (total:dissolved)',default=1.0_rk)
      if (self%last_layer/=3) then
         call self%get_parameter(self%relax,'relax','1/d','rate of relaxation towards equilibrium concentration profile')
         call self%get_parameter(self%minD, 'minD','m',  'minimum depth of bottom interface of deepest layer')
      end if

      ! Register state variable for bottom-most pelagic concentration.
      call self%register_state_dependency(self%id_pel,trim(composition)//'_pel','mmol/m^3','pelagic '//trim(long_name))
      call self%register_diagnostic_variable(self%id_pbf,'pb_flux','mmol/m^2/day','pelagic-benthic flux')      
      ! Dependencies
      call self%register_state_dependency(self%id_D1m, 'D1m', 'm','depth of bottom interface of oxygenated layer',standard_variable=depth_of_bottom_interface_of_layer_1)
      call self%register_state_dependency(self%id_D2m, 'D2m', 'm','depth of bottom interface of oxidized layer',standard_variable=depth_of_bottom_interface_of_layer_2)
      call self%register_dependency(self%id_Dtot,depth_of_sediment_column)
      call self%register_dependency(self%id_poro,sediment_porosity)
      call self%register_dependency(self%id_diff(1),diffusivity_in_sediment_layer_1)
      call self%register_dependency(self%id_diff(2),diffusivity_in_sediment_layer_2)
      call self%register_dependency(self%id_diff(3),diffusivity_in_sediment_layer_3)
      call self%register_dependency(self%id_cmix,pelagic_benthic_transfer_constant)

      ! Create model that computes concentrations per benthic layer.
      allocate(profile)
      call self%add_child(profile,'per_layer',configunit=configunit)
      profile%ads = self%ads
      profile%last_layer = self%last_layer
      call profile%register_diagnostic_variable(profile%id_layers(1),trim(composition)//'1','mmol/m^2','total '//trim(long_name)//' in oxygenated layer (absorbed + dissolved)',act_as_state_variable=.true.,domain=domain_bottom,output=output_none)
      call profile%register_diagnostic_variable(profile%id_layers(2),trim(composition)//'2','mmol/m^2','total '//trim(long_name)//' in oxidized layer (absorbed + dissolved)',act_as_state_variable=.true.,domain=domain_bottom,output=output_none)
      call profile%register_diagnostic_variable(profile%id_layers(3),trim(composition)//'3','mmol/m^2','total '//trim(long_name)//' in anoxic layer (absorbed + dissolved)',act_as_state_variable=.true.,domain=domain_bottom,output=output_none)
      call profile%register_diagnostic_variable(profile%id_layers_pw(1),trim(composition)//'1_pw','mmol/m^2','dissolved '//trim(long_name)//' in oxygenated layer',act_as_state_variable=.true.,domain=domain_bottom,output=output_none)
      call profile%register_diagnostic_variable(profile%id_layers_pw(2),trim(composition)//'2_pw','mmol/m^2','dissolved '//trim(long_name)//' in oxidized layer',act_as_state_variable=.true.,domain=domain_bottom,output=output_none)
      call profile%register_diagnostic_variable(profile%id_layers_pw(3),trim(composition)//'3_pw','mmol/m^2','dissolved '//trim(long_name)//' in anoxic layer',act_as_state_variable=.true.,domain=domain_bottom,output=output_none)
      call profile%register_dependency(profile%id_D1m, 'D1m', 'm','depth of bottom interface of oxygenated layer')
      call profile%register_dependency(profile%id_D2m, 'D2m', 'm','depth of bottom interface of oxidized layer')
      call profile%register_dependency(profile%id_Dtot,depth_of_sediment_column)
      call profile%register_dependency(profile%id_poro,sediment_porosity)
      call profile%register_dependency(profile%id_tot,trim(composition)//'_int','mmol/m^2',trim(long_name)//', depth-integrated')
      call profile%request_coupling(profile%id_D1m,'D1m')
      call profile%request_coupling(profile%id_D2m,'D2m')
      call profile%request_coupling(profile%id_tot,composition)

      ! Make sure that sources-sinks of layer-specific mass are counted in mass budgets.
      select case (composition)
         case ('c')
            call profile%add_to_aggregate_variable(standard_variables%total_carbon,profile%id_layers(1))
            call profile%add_to_aggregate_variable(standard_variables%total_carbon,profile%id_layers(2))
            call profile%add_to_aggregate_variable(standard_variables%total_carbon,profile%id_layers(3))
            call profile%add_to_aggregate_variable(standard_variables%total_carbon,profile%id_layers_pw(1))
            call profile%add_to_aggregate_variable(standard_variables%total_carbon,profile%id_layers_pw(2))
            call profile%add_to_aggregate_variable(standard_variables%total_carbon,profile%id_layers_pw(3))
         case ('n'); long_name = 'nitrogen'
            call profile%add_to_aggregate_variable(standard_variables%total_nitrogen,profile%id_layers(1))
            call profile%add_to_aggregate_variable(standard_variables%total_nitrogen,profile%id_layers(2))
            call profile%add_to_aggregate_variable(standard_variables%total_nitrogen,profile%id_layers(3))
            call profile%add_to_aggregate_variable(standard_variables%total_nitrogen,profile%id_layers_pw(1))
            call profile%add_to_aggregate_variable(standard_variables%total_nitrogen,profile%id_layers_pw(2))
            call profile%add_to_aggregate_variable(standard_variables%total_nitrogen,profile%id_layers_pw(3))
         case ('p'); long_name = 'phosphorus'
            call profile%add_to_aggregate_variable(standard_variables%total_phosphorus,profile%id_layers(1))
            call profile%add_to_aggregate_variable(standard_variables%total_phosphorus,profile%id_layers(2))
            call profile%add_to_aggregate_variable(standard_variables%total_phosphorus,profile%id_layers(3))
            call profile%add_to_aggregate_variable(standard_variables%total_phosphorus,profile%id_layers_pw(1))
            call profile%add_to_aggregate_variable(standard_variables%total_phosphorus,profile%id_layers_pw(2))
            call profile%add_to_aggregate_variable(standard_variables%total_phosphorus,profile%id_layers_pw(3))
         case ('s'); long_name = 'silicate'
            call profile%add_to_aggregate_variable(standard_variables%total_silicate,profile%id_layers(1))
            call profile%add_to_aggregate_variable(standard_variables%total_silicate,profile%id_layers(2))
            call profile%add_to_aggregate_variable(standard_variables%total_silicate,profile%id_layers(3))
            call profile%add_to_aggregate_variable(standard_variables%total_silicate,profile%id_layers_pw(1))
            call profile%add_to_aggregate_variable(standard_variables%total_silicate,profile%id_layers_pw(2))
            call profile%add_to_aggregate_variable(standard_variables%total_silicate,profile%id_layers_pw(3))
      end select

      ! Couple to layer-specific "sinks minus sources".
      call self%register_dependency(self%id_sms(1),'sms_l1','mmol/m^2/s','sources-sinks of total '//trim(long_name)//' in oxygenated layer')
      call self%register_dependency(self%id_sms(2),'sms_l2','mmol/m^2/s','sources-sinks of total '//trim(long_name)//' in oxidized layer')
      call self%register_dependency(self%id_sms(3),'sms_l3','mmol/m^2/s','sources-sinks of total '//trim(long_name)//' in anoxic layer')
      call self%request_coupling('sms_l1','per_layer/'//trim(composition)//'1_sms_tot')
      call self%request_coupling('sms_l2','per_layer/'//trim(composition)//'2_sms_tot')
      call self%request_coupling('sms_l3','per_layer/'//trim(composition)//'3_sms_tot')

      ! Couple to layer-specific "sinks minus sources" for matter in pore water.
      call self%register_dependency(self%id_pw_sms(1),'pw_sms_l1','mmol/m^2/s','sources-sinks of dissolved '//trim(long_name)//' in oxygenated layer')
      call self%register_dependency(self%id_pw_sms(2),'pw_sms_l2','mmol/m^2/s','sources-sinks of dissolved '//trim(long_name)//' in oxidized layer')
      call self%register_dependency(self%id_pw_sms(3),'pw_sms_l3','mmol/m^2/s','sources-sinks of dissolved '//trim(long_name)//' in anoxic layer')
      call self%request_coupling('pw_sms_l1','per_layer/'//trim(composition)//'1_pw_sms_tot')
      call self%request_coupling('pw_sms_l2','per_layer/'//trim(composition)//'2_pw_sms_tot')
      call self%request_coupling('pw_sms_l3','per_layer/'//trim(composition)//'3_pw_sms_tot')
   end subroutine benthic_dissolved_matter_initialize

   subroutine benthic_dissolved_matter_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_benthic_column_dissolved_matter),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: c_pel,c_int,sms_l1,sms_l2,sms,pw_sms_l1,pw_sms_l2
      real(rk) :: d1,d2,d3
      real(rk) :: c_bot1_eq,c_int1_eq,H1_eq
      real(rk) :: c_bot2_eq,c_int2_eq,H2_eq
      real(rk) :: c_bot3_eq,c_int3_eq
      real(rk) :: c_int_eq
      real(rk) :: norm_res_int,P_res_int

      real(rk) :: diff1,diff2,diff3,poro,cmix

      _HORIZONTAL_LOOP_BEGIN_

      ! Retrieve column-integrated mass, layer-specific source-sink terms, and the lowermost pelagic concentration.
      _GET_HORIZONTAL_(self%id_tot,c_int)
      _GET_HORIZONTAL_(self%id_sms(1),sms_l1)
      _GET_HORIZONTAL_(self%id_sms(2),sms_l2)
      _GET_HORIZONTAL_(self%id_pw_sms(1),pw_sms_l1)
      _GET_HORIZONTAL_(self%id_pw_sms(2),pw_sms_l2)
      _GET_(self%id_pel,c_pel)

      ! Sink-source terms are always /s, and we need /d.
      sms_l1 = (sms_l1+pw_sms_l1)*86400._rk
      sms_l2 = (sms_l2+pw_sms_l2)*86400._rk

      ! Retrieve physical properties of sediment column:
      ! layer depths, porosity, pelagic-benthic transfer coefficient, diffusivities.
      _GET_HORIZONTAL_(self%id_D1m,d1)
      _GET_HORIZONTAL_(self%id_D2m,d2)
      _GET_HORIZONTAL_(self%id_Dtot,d3)
      _GET_HORIZONTAL_(self%id_poro,poro)
      _GET_HORIZONTAL_(self%id_cmix,cmix)
      _GET_HORIZONTAL_(self%id_diff(1),diff1)
      _GET_HORIZONTAL_(self%id_diff(2),diff2)
      _GET_HORIZONTAL_(self%id_diff(3),diff3)

      ! Column-integrated source-sink terms
      sms = sms_l1 + sms_l2

      ! Estimate steady state concentration at sediment interface from current pelagic concentration
      ! (typically at centre of lowermost pelagic layer), and column-integrated production/destruction term.
      ! JB: the logic behind the expressions below is UNKNOWN (in original ERSEM, this was handled by the modconc subroutine)
      if (sms>0._rk) then
         ! Sediment column produces tracer. Steady state concentration at sediment interface > lowermost pelagic concentration.
         ! One way to arrive at this formulation is to assume no production or destruction of the tracer
         ! within the lowermost pelagic layer. In that case, the equilibrium concentration profile near the bed has a slope
         ! equal to sms/diffusivity. The change between centre of the pelagic layer and the sediment interface equals
         ! sms/D*H/2, with D representing diffusivity within the lowermost pelagic layer and H representing the height of the 
         ! lowermost pelagic layer. Thus, cmix must equal H/2/D. In practice, however, H is defined by the vertical grid and
         ! and D is defined by the level of turbulence - it then seems unreasonable to treat cmix=H/2/D as constant - JB 4/2/2015.
         c_pel = c_pel + cmix*sms
      else
         ! Sediment column destroys tracer. Steady state concentration at sediment interface < lowermost pelagic concentration.
         ! Note cmix*sms<0, thus c_pel/(c_pel - cmix*sms)<1
         c_pel = c_pel*c_pel/(c_pel - cmix*sms)
      end if

      if (self%last_layer==1) then
         ! Oxygenated layer: compute steady-state layer height H1_eq and layer integral c_int1_eq
         call compute_final_equilibrium_profile(diff1,c_pel,sms_l1,sms_l2,d3,H1_eq,c_int1_eq)

         ! Benthic dynamics: relax depth-integrated mass towards equilibrium value
         c_int_eq = poro*self%ads(1)*c_int1_eq
         _SET_BOTTOM_ODE_(self%id_tot,(c_int_eq-c_int)/self%relax)

         ! Net change in benthos must equal local production - surface exchange.
         ! Thus, surface exchange = local production - net change (net change = relaxation)
         _SET_BOTTOM_EXCHANGE_(self%id_pel,sms-(c_int_eq-c_int)/self%relax)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_pbf,sms-(c_int_eq-c_int)/self%relax)
      
         ! Relax depth of first/oxic layer towards equilibrium value (H1_eq)
         _SET_BOTTOM_ODE_(self%id_D1m,(max(self%minD,H1_eq)-d1)/self%relax)
      elseif (self%last_layer==2) then
         ! Oxygenated layer: compute steady-state concentration at bottom interface c_bot1_eq and layer integral c_int1_eq
         call compute_equilibrium_profile(diff1,c_pel,sms_l1,sms_l2,d1,c_bot1_eq,c_int1_eq)
         ! Oxidized layer: compute steady-state layer height H2_eq and layer integral c_int2_eq
         call compute_final_equilibrium_profile(diff2,c_bot1_eq,sms_l2,0.0_rk,d3-d1,H2_eq,c_int2_eq)

         ! Benthic dynamics: relax depth-integrated mass towards equilibrium value
         c_int_eq = poro*(self%ads(1)*c_int1_eq+self%ads(2)*c_int2_eq)
         _SET_BOTTOM_ODE_(self%id_tot,(c_int_eq-c_int)/self%relax)

         ! Net change in benthos must equal local production - surface exchange.
         ! Thus, surface exchange = local production - net change (net change = relaxation)
         _SET_BOTTOM_EXCHANGE_(self%id_pel,sms-(c_int_eq-c_int)/self%relax)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_pbf,sms-(c_int_eq-c_int)/self%relax)

         ! Relax depth of bottom interface of second/oxidised layer towards equilibrium value (d1+H2_eq)
         _SET_BOTTOM_ODE_(self%id_D2m,(max(self%minD,d1+H2_eq)-d2)/self%relax)
      else
         ! Oxygenated layer: compute steady-state concentration at bottom interface c_bot1_eq and layer integral c_int1_eq
         call compute_equilibrium_profile(diff1,c_pel,    sms_l1,sms_l2,d1,   c_bot1_eq,c_int1_eq)
         ! Oxidized layer: compute steady-state concentration at bottom interface c_bot2_eq and layer integral c_int2_eq
         call compute_equilibrium_profile(diff2,c_bot1_eq,sms_l2,0.0_rk,d2-d1,c_bot2_eq,c_int2_eq)
         ! Anoxic layer: no sources or sinks: homogeneous equilibrium concentration c_bot2_eq
         c_int3_eq = (d3-d2)*c_bot2_eq
         c_int_eq = poro*(self%ads(1)*c_int1_eq+self%ads(2)*c_int2_eq+self%ads(3)*c_int3_eq)

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
         call compute_equilibrium_profile(diff1,0.0_rk,   d1,   d3-d1, d1,   c_bot1_eq,c_int1_eq)
         call compute_equilibrium_profile(diff2,c_bot1_eq,d2-d1,d3-d2, d2-d1,c_bot2_eq,c_int2_eq)
         call compute_equilibrium_profile(diff3,c_bot2_eq,d3-d2,0.0_rk,d3-d2,c_bot3_eq,c_int3_eq)
         norm_res_int = poro*(self%ads(1)*c_int1_eq+self%ads(2)*c_int2_eq+self%ads(3)*c_int3_eq)
         P_res_int = (c_int-c_int_eq)/norm_res_int*d3
         _SET_BOTTOM_EXCHANGE_(self%id_pel,sms+P_res_int) ! Equilibrium flux = sms, residual flux = P_res_int
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_pbf,sms+P_res_int)
         _SET_BOTTOM_ODE_(self%id_tot,-P_res_int)         ! Local sources-sinks (sms) minus surface flux (sms+P_res_int)

      end if

      _HORIZONTAL_LOOP_END_
   end subroutine benthic_dissolved_matter_do_bottom

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
      ! The concentration at the bottom equal the value of the parabola at depth D:
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
      ! As $a$ always occurs multiplied with $D^2$, and $b$ always multiplied with $D$,
      ! we define the combined constants $a_D2 = a D^2$ and $b_D = b D$. This avoids division by 0 when
      ! computing $a$ while $D$ tends to zero.
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
      
      if (Dmax*(P+2*P_deep)>-2*sigma*c0) then   ! Dmax<-2 sigma c0/(P+2*P_deep), rearranged for P+2*P_deep<0 close to 0. Result also picks up P+2*P_deep>0.
         ! Destruction within layer is too low (or layer experiences net production, i.e., P>0).
         ! If we would impose zero concentration at the layer bottom, the layer would extend beyond maximum depth.
         ! Fix depth at maximum depth and use parabola with zero flux but non-zero concentration at bottom interface.
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

      real(rk) :: c_int
      real(rk) :: factor
      real(rk) :: d(3),D2m,d_totX,poro

      _HORIZONTAL_LOOP_BEGIN_

         ! Get depth-integrated concentration.
         _GET_HORIZONTAL_(self%id_tot,c_int)

         ! Layer depths and porosity.
         _GET_HORIZONTAL_(self%id_D1m,d(1))
         _GET_HORIZONTAL_(self%id_D2m,D2m)
         _GET_HORIZONTAL_(self%id_Dtot,d_totX)
         _GET_HORIZONTAL_(self%id_poro,poro)

         ! Compute layer thicknesses from depth of bottom interface of individual layers.
         d(2) = D2m-d(1)
         d(3) = d_totX-D2m

         if (self%last_layer==1) then
            ! All in oxygenated layer
            ! Typically used for oxygen
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers(1),c_int)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers(2),0.0_rk)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers(3),0.0_rk)

            ! Pore water contents: layer contents divided by adsorption [total:dissolved]
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers_pw(1),c_int/self%ads(1))
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers_pw(2),0.0_rk)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers_pw(3),0.0_rk)
         elseif (self%last_layer==2) then
            ! Vertically homogeneous in oxygenated layer, quadratically decreasing in oxidized layer (zero concentration at bottom interface)
            ! Typically used for nitrate

            ! Mean concentration in pore water (matter/m3)
            factor = c_int/(poro*(self%ads(1)*d(1)+self%ads(2)*d(2)/3.0_rk))

            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers(1),poro*self%ads(1)*d(1)*factor)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers(2),poro*self%ads(2)*d(2)/3.0_rk*factor)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers(3),0.0_rk)
            
            ! Pore water contents: layer contents divided by adsorption [total:dissolved]
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers_pw(1),poro*d(1)*factor)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers_pw(2),poro*d(2)/3.0_rk*factor)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers_pw(3),0.0_rk)
         else
            ! Vertically homogeneous in all individual layers (oxygenated, oxidized, anoxic).
            ! Typically used for all tracers but oxygen and nitrate.

            ! Mean concentration in pore water (matter/m3)
            factor = c_int/sum(poro*self%ads*d)

            ! Layer contents (matter/m2)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers(1),poro*self%ads(1)*d(1)*factor)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers(2),poro*self%ads(2)*d(2)*factor)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers(3),poro*self%ads(3)*d(3)*factor)

            ! Pore water contents: layer contents divided by adsorption [total:dissolved]
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers_pw(1),poro*d(1)*factor)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers_pw(2),poro*d(2)*factor)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers_pw(3),poro*d(3)*factor)
         end if

      _HORIZONTAL_LOOP_END_

   end subroutine dissolved_matter_per_layer_do_bottom

end module
