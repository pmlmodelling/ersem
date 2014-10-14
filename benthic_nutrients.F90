#include "fabm_driver.h"

module pml_ersem_benthic_nutrients

   use fabm_types
   use fabm_standard_variables

   use pml_ersem_shared

   implicit none

   private

   ! Define standard variables that enable automatic coupling between variables across models related to 3-layer benthos.
   ! (variables that are assigned the same "standard variable" identity will be coupled by FABM)
   ! The alternative would be to do manual coupling for each of these, for each model instance, but that would be very verbose.
   type (type_horizontal_standard_variable),parameter :: diffusivity_in_sediment_layer_1 = type_horizontal_standard_variable(name='diffusivity_in_sediment_layer_1',units='m^2/d')
   type (type_horizontal_standard_variable),parameter :: diffusivity_in_sediment_layer_2 = type_horizontal_standard_variable(name='diffusivity_in_sediment_layer_2',units='m^2/d')
   type (type_horizontal_standard_variable),parameter :: diffusivity_in_sediment_layer_3 = type_horizontal_standard_variable(name='diffusivity_in_sediment_layer_3',units='m^2/d')
   type (type_horizontal_standard_variable),parameter :: sediment_porosity = type_horizontal_standard_variable(name='sediment_porosity',units='-')
   type (type_horizontal_standard_variable),parameter :: depth_of_bottom_interface_of_layer_1 = type_horizontal_standard_variable(name='depth_of_bottom_interface_of_layer_1',units='m')
   type (type_horizontal_standard_variable),parameter :: depth_of_bottom_interface_of_layer_2 = type_horizontal_standard_variable(name='depth_of_bottom_interface_of_layer_2',units='m')
   type (type_horizontal_standard_variable),parameter :: depth_of_sediment_column = type_horizontal_standard_variable(name='depth_of_sediment_column',units='m')
   type (type_horizontal_standard_variable),parameter :: pelagic_benthic_transfer_constant = type_horizontal_standard_variable(name='pelagic_benthic_transfer_constant',units='d/m')

   ! Model that specifies the structure of the 3-layer sediment column
   type,extends(type_base_model),public :: type_pml_ersem_benthic_3layer_column
      type (type_bottom_state_variable_id)          :: id_D1m,id_D2m
      type (type_horizontal_diagnostic_variable_id) :: id_poro,id_Dtot,id_diff(3),id_EDZ_mixX

      real(rk) :: EDZ_1X,EDZ_2X,EDZ_3X,d_totX
      real(rk) :: qPWX,EDZ_mixX
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

   ! Model for dissolved matter in sediment, using idealized equilibrium profiles
   ! to determine pelagic-benthic diffusive flux, and to determine target depth of
   ! first and second layer.
   type,extends(type_base_model),public :: type_pml_ersem_benthic_dissolved_matter
      type (type_bottom_state_variable_id) :: id_tot
      type (type_state_variable_id)        :: id_pel
      type (type_bottom_state_variable_id) :: id_D1m,id_D2m
      type (type_horizontal_dependency_id) :: id_sms(3),id_Dtot,id_poro,id_cmix,id_diff(3)

      real(rk) :: ads(3)
      real(rk) :: relax, minD
      integer :: type
   contains
      procedure :: initialize => benthic_dissolved_matter_initialize
      procedure :: do_bottom  => benthic_dissolved_matter_do_bottom
   end type

   ! 
   type,extends(type_base_model) :: type_dissolved_matter_per_layer
      type (type_horizontal_dependency_id)          :: id_tot
      type (type_horizontal_dependency_id)          :: id_D1m,id_D2m,id_Dtot,id_poro
      type (type_horizontal_diagnostic_variable_id) :: id_layers(3)

      real(rk) :: ads(3)
      integer :: type
   contains
      procedure :: do_bottom => dissolved_matter_per_layer_do_bottom
   end type

   type,extends(type_base_model),public :: type_pml_ersem_benthic_nitrification
      type (type_bottom_state_variable_id) :: id_K3n,id_K4n,id_G2o
      type (type_state_variable_id)        :: id_N4n
      type (type_dependency_id)            :: id_ETW,id_phx
      type (type_horizontal_dependency_id) :: id_D1m

      real(rk) :: q10nitX,hM4M3X,sM4M3X,xno3X
      integer :: ISWphx
   contains
      procedure :: initialize => benthic_nitrification_initialize
      procedure :: do_bottom  => benthic_nitrification_do_bottom
   end type

contains   

   subroutine benthic_dissolved_matter_initialize(self,configunit)
      class (type_pml_ersem_benthic_dissolved_matter),intent(inout),target :: self
      integer,                                        intent(in)           :: configunit

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
      call self%get_parameter(self%type,  'type','-','layer presence',default=3)
      self%ads = 1.0_rk
      call self%get_parameter(self%ads(1),'ads1','-','adsorption in layer 1 (total:dissolved)',default=1.0_rk)
      if (self%type>1) call self%get_parameter(self%ads(2),'ads2','-','adsorption in layer 2 (total:dissolved)',default=1.0_rk)
      if (self%type>2) call self%get_parameter(self%ads(3),'ads3','-','adsorption in layer 3 (total:dissolved)',default=1.0_rk)
      if (self%type/=3) then
         call self%get_parameter(self%relax,'relax','/d','diffusion time scale for benthic/pelagic interface')
         call self%get_parameter(self%minD, 'minD','m',  'minimum depth of bottom interface of controlled layer')
      end if

      ! Register state variable for bottom-most pelagic concentration.
      call self%register_state_dependency(self%id_pel,trim(composition)//'_pel','mmol/m^3','pelagic '//trim(long_name))

      ! Dependencies
      call self%register_state_dependency(self%id_D1m, 'D1m', 'm','depth of bottom interface of layer 1',standard_variable=depth_of_bottom_interface_of_layer_1)
      call self%register_state_dependency(self%id_D2m, 'D2m', 'm','depth of bottom interface of layer 2',standard_variable=depth_of_bottom_interface_of_layer_2)
      call self%register_dependency(self%id_Dtot,'Dtot','m','depth of sediment column',standard_variable=depth_of_sediment_column)
      call self%register_dependency(self%id_poro,'poro','-','porosity',standard_variable=sediment_porosity)
      call self%register_dependency(self%id_diff(1),'diff1','m^2/d','diffusivity in layer 1',standard_variable=diffusivity_in_sediment_layer_1)
      call self%register_dependency(self%id_diff(2),'diff2','m^2/d','diffusivity in layer 2',standard_variable=diffusivity_in_sediment_layer_2)
      call self%register_dependency(self%id_diff(3),'diff3','m^2/d','diffusivity in layer 3',standard_variable=diffusivity_in_sediment_layer_3)
      call self%register_dependency(self%id_cmix,'cmix','d/m','equilibrium diffusive speed between sediment surface water',standard_variable=pelagic_benthic_transfer_constant)
      
      ! Create model that computes concentrations per benthic layer.
      allocate(profile)
      call self%add_child(profile,'per_layer',configunit=configunit)
      profile%ads = self%ads
      profile%type = self%type
      call profile%register_diagnostic_variable(profile%id_layers(1),trim(composition)//'1','mmol/m^2',trim(long_name)//' in layer 1')
      call profile%register_diagnostic_variable(profile%id_layers(2),trim(composition)//'2','mmol/m^2',trim(long_name)//' in layer 2')
      call profile%register_diagnostic_variable(profile%id_layers(3),trim(composition)//'3','mmol/m^2',trim(long_name)//' in layer 3')
      call profile%act_as_state_variable(profile%id_layers(1))
      call profile%act_as_state_variable(profile%id_layers(2))
      call profile%act_as_state_variable(profile%id_layers(3))
      call profile%register_dependency(profile%id_D1m, 'D1m', 'm','depth of bottom interface of 1st layer')
      call profile%register_dependency(profile%id_D2m, 'D2m', 'm','depth of bottom interface of 2nd layer')
      call profile%register_dependency(profile%id_Dtot,'Dtot','m','depth of sediment column')
      call profile%register_dependency(profile%id_poro,'poro','-','porosity')
      call profile%register_dependency(profile%id_tot,trim(composition)//'_int','mmol/m^2',trim(long_name)//', depth-integrated')
      call profile%request_coupling(profile%id_D1m,'D1m')
      call profile%request_coupling(profile%id_D2m,'D2m')
      call profile%request_coupling(profile%id_Dtot,'Dtot')
      call profile%request_coupling(profile%id_poro,'poro')
      call profile%request_coupling(profile%id_tot,composition)

      ! Couple to layer-specific "sinks minus sources".
      call self%register_dependency(self%id_sms(1),'sms_l1','mmol/m^2/s',trim(long_name)//' sinks-sources in layer 1')
      call self%register_dependency(self%id_sms(2),'sms_l2','mmol/m^2/s',trim(long_name)//' sinks-sources in layer 2')
      call self%register_dependency(self%id_sms(3),'sms_l3','mmol/m^2/s',trim(long_name)//' sinks-sources in layer 3')
      call self%request_coupling('sms_l1','per_layer/'//trim(composition)//'1_sms')
      call self%request_coupling('sms_l2','per_layer/'//trim(composition)//'2_sms')
      call self%request_coupling('sms_l3','per_layer/'//trim(composition)//'3_sms')
   end subroutine benthic_dissolved_matter_initialize

   subroutine benthic_dissolved_matter_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_pml_ersem_benthic_dissolved_matter),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: c_pel,c_int,sms_l1,sms_l2
      real(rk) :: d1,d2,d3
      real(rk) :: c_bot1_eq,c_int1_eq,H1_eq
      real(rk) :: c_bot2_eq,c_int2_eq,H2_eq
      real(rk) :: c_bot3_eq,c_int3_eq
      real(rk) :: c_int_eq
      real(rk) :: norm_res_int,P_res_int

      real(rk) :: diff1,diff2,diff3,poro,cmix

      _HORIZONTAL_LOOP_BEGIN_

      _GET_(self%id_pel,c_pel)
      _GET_HORIZONTAL_(self%id_tot,c_int)
      _GET_HORIZONTAL_(self%id_sms(1),sms_l1)
      _GET_HORIZONTAL_(self%id_sms(2),sms_l2)

      ! Sink-source terms are always /s, and we need /d.
      sms_l1 = sms_l1*86400._rk
      sms_l2 = sms_l2*86400._rk

      _GET_HORIZONTAL_(self%id_D1m,d1)
      _GET_HORIZONTAL_(self%id_D2m,d2)
      _GET_HORIZONTAL_(self%id_Dtot,d3)
      _GET_HORIZONTAL_(self%id_poro,poro)
      _GET_HORIZONTAL_(self%id_cmix,cmix)
      _GET_HORIZONTAL_(self%id_diff(1),diff1)
      _GET_HORIZONTAL_(self%id_diff(2),diff2)
      _GET_HORIZONTAL_(self%id_diff(3),diff3)
      
      if (self%type==1) then
         ! Layer 1: compute steady-state layer height H1_eq and layer integral c_int1_eq
         call compute_final_equilibrium_profile(diff1,modconc(c_pel,sms_l1,cmix),sms_l1,d3,H1_eq,c_int1_eq)

         ! Benthic dynamics: relax towards equilibrium value
         c_int_eq = poro*self%ads(1)*c_int1_eq
         _SET_BOTTOM_ODE_(self%id_tot,(c_int_eq-c_int)/self%relax)

         ! Net change in benthos must equal local production - surface exchange.
         ! Thus, surface exchange = local production - net change (net change = relaxation)
         _SET_BOTTOM_EXCHANGE_(self%id_pel,sms_l1-(c_int_eq-c_int)/self%relax)

         ! Relax depth of layer towards equilibrium value (H1_eq)
         _SET_BOTTOM_ODE_(self%id_D1m,(max(self%minD,H1_eq)-d1)/self%relax)
      elseif (self%type==2) then
         ! Layer 1: compute steady-state concentration at bottom interface c_bot1_eq and layer integral c_int1_eq
         call compute_equilibrium_profile(d1,diff1,modconc(c_pel,sms_l1+sms_l2,cmix),sms_l1,sms_l2,c_bot1_eq,c_int1_eq)
         ! Layer 2: compute steady-state layer height H2_eq and layer integral c_int2_eq
         call compute_final_equilibrium_profile(diff2,c_bot1_eq,sms_l2,d3-d1,H2_eq,c_int2_eq)

         ! Benthic nitrate dynamics: relax towards equilibrium value
         c_int_eq = poro*(self%ads(1)*c_int1_eq+self%ads(2)*c_int2_eq)
         _SET_BOTTOM_ODE_(self%id_tot,(c_int_eq-c_int)/self%relax)

         ! Net change in benthos must equal local production - surface exchange.
         ! Thus, surface exchange = local production - net change (net change = relaxation)
         _SET_BOTTOM_EXCHANGE_(self%id_pel,sms_l1+sms_l2-(c_int_eq-c_int)/self%relax)
      
         ! Relax depth of bottom interface of second/oxidised layer towards equilibrium value (d1+H2_eq)
         _SET_BOTTOM_ODE_(self%id_D2m,(max(self%minD,d1+H2_eq)-d2)/self%relax)
      else
         ! Layer 1: compute steady-state concentration at bottom interface c_bot1_eq and layer integral c_int1_eq
         call compute_equilibrium_profile(d1,diff1,modconc(c_pel,sms_l1+sms_l2,cmix),sms_l1,sms_l2,c_bot1_eq,c_int1_eq)
         ! Layer 2: compute steady-state concentration at bottom interface c_bot2_eq and layer integral c_int2_eq
         call compute_equilibrium_profile(d2-d1,diff2,c_bot1_eq,sms_l2,0.0_rk,c_bot2_eq,c_int2_eq)
         ! Layer 3: no sources or sinks: homogeneous equilibrium concentration c_bot2_eq
         c_int3_eq = (d3-d2)*c_bot2_eq
         c_int_eq = poro*(self%ads(1)*c_int1_eq+self%ads(2)*c_int2_eq+self%ads(3)*c_int3_eq)

         ! The equilibrium depth-integrated mass c_int_eq usually differs from the current depth-integrated mass K4nP.
         ! We can view the actual [unknown] pore water concentration profile as the sum of the equilibrium profile
         ! and a residual profile. The latter has a vertical integral equal to the difference between actual mass and equilibrium mass.
         ! The rate at which the residual mass is exchanged over the benthic-pelagic interface is equal to the product of
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
         call compute_equilibrium_profile(d1,   diff1,0.0_rk,   d1,   d3-d1, c_bot1_eq,c_int1_eq)
         call compute_equilibrium_profile(d2-d1,diff2,c_bot1_eq,d2-d1,d3-d2, c_bot2_eq,c_int2_eq)
         call compute_equilibrium_profile(d3-d2,diff3,c_bot2_eq,d3-d2,0.0_rk,c_bot3_eq,c_int3_eq)
         norm_res_int = poro*(self%ads(1)*c_int1_eq+self%ads(2)*c_int2_eq+self%ads(3)*c_int3_eq)
         P_res_int = (c_int-c_int_eq)/norm_res_int*d3
         _SET_BOTTOM_EXCHANGE_(self%id_pel,sms_l1+sms_l2+P_res_int) ! Equilibrium flux = sms_l1+sms_l2, residual flux = P_res_int
         _SET_BOTTOM_ODE_(self%id_tot,-P_res_int)
      end if

      _HORIZONTAL_LOOP_END_
   end subroutine benthic_dissolved_matter_do_bottom

   subroutine benthic_nitrification_initialize(self,configunit)
      class (type_pml_ersem_benthic_nitrification),intent(inout),target :: self
      integer,                                     intent(in)           :: configunit

      ! Set time unit to d-1. This implies that all rates (sink/source terms) are given in d-1.
      self%dt = 86400._rk

      call self%get_parameter(self%q10nitX,'q10nitX','-',            'Q_10 temperature coefficient mfor nitrification')
      call self%get_parameter(self%hM4M3X, 'hM4M3X', 'mmol/m^3',     'Michaelis-Menten constant for nitrate limitation of nitrification')
      call self%get_parameter(self%ISWphx, 'ISWphx', '',             'pH influence on nitrification',default=0)
      call self%get_parameter(self%sM4M3X, 'sM4M3X', '1/d',          'maximum nitrification rate')
      call self%get_parameter(self%xno3X,  'xno3X',  'mol O_2/mol N','oxygen consumed per nitrate produced in nitrification')

      call self%register_state_dependency(self%id_K3n,'K3n','mmol/m^2','nitrate')
      call self%register_state_dependency(self%id_K4n,'K4n','mmol/m^2','ammonium')
      call self%register_state_dependency(self%id_G2o,'G2o','mmol/m^2','oxygen')
      call self%register_state_dependency(self%id_N4n,'N4n','mmol N/m^3','pelagic ammonium')
      call self%register_dependency(self%id_D1m,'D1m','m','depth of bottom interface of 1st layer',standard_variable=depth_of_bottom_interface_of_layer_1)

      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      call self%register_dependency(self%id_phx,standard_variables%ph_reported_on_total_scale)
   end subroutine benthic_nitrification_initialize

   subroutine benthic_nitrification_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !INPUT PARAMETERS:
   class (type_pml_ersem_benthic_nitrification),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_

   real(rk) :: K3n,K3nP,K4nP,N4n
   real(rk) :: ETW,phx
   real(rk) :: Mu_m,eT,eN,Fph,jM4M3n,D1m

   _HORIZONTAL_LOOP_BEGIN_
      _GET_HORIZONTAL_(self%id_K3n,K3n) !TODO background
      _GET_HORIZONTAL_(self%id_K3n,K3nP)
      _GET_HORIZONTAL_(self%id_K4n,K4nP)
      _GET_(self%id_N4n,N4n)

      _GET_HORIZONTAL_(self%id_D1m,D1m)
      
      _GET_(self%id_ETW,ETW)

      !* Nitrification (layer 1 only):
      ! Reaction: (NH4+,OH-) + 2*O2 -> (NO3-,H+) + 2*H2O
      ! treated as first order reaction for NH4 in the oxygenated layer
      ! Average nitrate density (mmol/m3): MU_m

      ! Average nitrate density in first/oxygenated layer
      MU_m = max(0.0_rk,K3n)/D1m

      ! Limitation factor for temperature (Q_10)
      eT = self%q10nitX**((ETW-10._rk)/10._rk)

      ! Limitation factor for nitrate (hyperbolic, 1 at zero nitrate, dropping to 0 at high nitrate).
      eN = self%hM4M3X/(self%hM4M3X+MU_m)

      ! Ph influence on nitrification - empirical equation
      if(self%ISWphx==1) then
         _GET_(self%id_phx,phx)
        Fph = min(2._rk,max(0._rk,0.6111_rk*phx-3.8889_rk))
        jM4M3n = Fph * self%sM4M3X * K4nP * eT * eN
      else
        jM4M3n = self%sM4M3X*K4nP*eT*eN
      end if

      ! Impose nitrification source-sink terms for benthic ammonium, nitrate, oxygen.
      _SET_BOTTOM_ODE_(self%id_K3n,jM4M3n)
      _SET_BOTTOM_ODE_(self%id_K4n,-jM4M3n)
      _SET_BOTTOM_ODE_(self%id_G2o,-self%xno3X*jM4M3n)
   _HORIZONTAL_LOOP_END_
   end subroutine benthic_nitrification_do_bottom

   subroutine initialize(self,configunit)
      class (type_pml_ersem_benthic_3layer_column),intent(inout),target :: self
      integer,                                     intent(in)           :: configunit
      
      ! Set time unit to d-1. This implies that all rates (sink/source terms) are given in d-1.
      self%dt = 86400._rk

      call self%get_parameter(self%EDZ_1X,'EDZ_1X','m^2/d','diffusivity in 1st (oxygenated) layer')
      call self%get_parameter(self%EDZ_2X,'EDZ_2X','m^2/d','diffusivity in 2nd (oxidized) layer')
      call self%get_parameter(self%EDZ_3X,'EDZ_3X','m^2/d','diffusivity in 3rd (anoxic) layer')
      call self%get_parameter(self%qPWX,'qPWX','-','fraction of pore water in the sediment')
      call self%get_parameter(self%EDZ_mixX,'EDZ_mixX','d/m','equilibrium diffusive speed between sediment surface water')
      call self%get_parameter(self%d_totX,'d_totX','m','depth of sediment column')

      call self%register_state_variable(self%id_D1m,'D1m','m','depth of bottom interface of 1st layer',standard_variable=depth_of_bottom_interface_of_layer_1)
      call self%register_state_variable(self%id_D2m,'D2m','m','depth of bottom interface of 2nd layer',standard_variable=depth_of_bottom_interface_of_layer_2)
      call self%register_diagnostic_variable(self%id_poro,'poro','-','porosity',standard_variable=sediment_porosity,missing_value=self%qPWX)
      call self%register_diagnostic_variable(self%id_Dtot,'Dtot','m','depth of sediment column',missing_value=self%d_totX,standard_variable=depth_of_sediment_column)
      call self%register_diagnostic_variable(self%id_diff(1),'diff1','m^2/s','diffusivity in layer 1',standard_variable=diffusivity_in_sediment_layer_1,missing_value=self%EDZ_1X)
      call self%register_diagnostic_variable(self%id_diff(2),'diff2','m^2/s','diffusivity in layer 2',standard_variable=diffusivity_in_sediment_layer_2,missing_value=self%EDZ_2X)
      call self%register_diagnostic_variable(self%id_diff(3),'diff3','m^2/s','diffusivity in layer 3',standard_variable=diffusivity_in_sediment_layer_3,missing_value=self%EDZ_3X)
      call self%register_diagnostic_variable(self%id_EDZ_mixX,'cmix','s/m','equilibrium diffusive speed between sediment surface water',standard_variable=pelagic_benthic_transfer_constant,missing_value=self%EDZ_mixX)
   end subroutine

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !INPUT PARAMETERS:
   class (type_pml_ersem_benthic_3layer_column),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
   
      real(rk) :: ETW,phx
      real(rk) :: D1m,D2m
      real(rk) :: irrenh
      real(rk) :: MU_m,eT,eN,Fph,jM4M3n
      real(rk) :: diff1,diff2,diff3,cmix
      real(rk) :: d1,d2,d3
      real(rk) :: poro,M1ads
      real(rk) :: vG2,vM3,vM4,vM1,vM11,vM5,vG3
      real(rk) :: jmun
      real(rk) :: jmu(7),jmi(7),jmn(7)
      real(rk) :: profO2(15),profNO3(15),profN(15),profP(15)
      real(rk) :: profS(15),profCO2(15),profD(15)
      real(rk) :: D1m0,D2m0
      real(rk) :: dum
      real(rk) :: H1_eq,H2_eq
      real(rk) :: c_bot1_eq, c_bot2_eq, c_bot3_eq, c_int1_eq, c_int2_eq, c_int3_eq, c_int_eq
      real(rk) :: norm_res_int,P_res_int
   
      _HORIZONTAL_LOOP_BEGIN_

      _GET_HORIZONTAL_(self%id_D1m,D1m)
      _GET_HORIZONTAL_(self%id_D2m,D2m)
!   
!      jmu = 0.0_rk
!      jmi = 0.0_rk
!      jmn = 0.0_rk
!      
!      irrenh = 1._rk !TODO: retrieve from benthic fauna
!      
!      D1m0 = 0.0_rk
!      D2m0 = 0.0_rk
!
!!!* Denitrification:
!!! The part pammonX (maximal) of anaerobic respiration reduces  NO3.
!!! From this NO3 reduction the part pdenitX goes to N2-gas.
!!
!!      MU_m = K3nP(K)/(D1m(K)+(D2m(K)-D1m(K))/3._rk)
!!      eN = eMM(MU_m/3._rk,hM3G4X)
!!! "borrowed" oxygen consumption in anaerobic layer:
!!
!!! to avoid numerical problems, use the recalculated value from the record
!!        jMIo2=jMI(4)
!!
!!! "expression in nitrate reduction (100%):
!!      jMIno3 = -jMIo2/(xno3X - xn2X*pdenitX)
!!      jM3M4n = jMIno3*eN*pammonX*(1._rk-pdenitX)
!!!      jjM3M4n(i) = jMIno3*eN*pammonX*(1.0d0-pdenitX)
!!      jM3G4n = jMIno3*eN*pammonX*     pdenitX
!!
!!      SK14n = SK14n +       jM3M4n
!!      SK3n(K)  = SK3n(K)  -       jM3M4n -      jM3G4n
!!      SG4n(K)  = SG4n(K)  +       jM3G4n
!!      SK6e  = SK6e  + xno3X*jM3M4n + (xno3X-xn2X)*jM3G4n
!!
!!      jMI(1)   = jMI(1)   +       jM3M4n
!!      jMI(6) = jMI(6) -       jM3M4n -      jM3G4n
!!      jMI(7)  = jMI(7)  +       jM3G4n
!!      jMI(4)  = jMI(4)  + xno3X*jM3M4n + (xno3X-xn2X)*jM3G4n
!
!!!--------------------------------------------
!!! Silicateregeneration - preliminary
!!!--------------------------------------------
!!
!!!aerobic layer: 0 .. D1m(K)
!!      rsil = sQ6M5X*partq(D9m(K),0._rk,D1m(K),d_totX)
!!
!!      SQ6s(K) = SQ6s(K) - rsil*Q6sP(K)
!!      SK5s(K) = SK5s(K) + rsil*Q6sP(K)
!!
!!      jMU(3) = jMU(3) + rsil*Q6sP(K)
!!
!!      SD9m(K) = SD9m(K) + (D1m(K)/2._rk-D9m(K))*rsil
!!#ifdef ERSEMDEBUG
!!      if(ersem_debugger.gt.0 .and. k .eq. kp_dbg) then
!!       PPWRITEDBGALLPRCS 'sil regen:',SD9m(K),D1m(K),D9m(K),rsil
!!      endif
!!#endif
!!
!!!lower layer: D1m(K) .. d_totX
!!      rsil = sQ6M5X*partq(D9m(K),D1m(K),d_totX,d_totX)
!!
!!      SQ6s(K) = SQ6s(K) - rsil*Q6sP(K)
!!      SK5s(K) = SK5s(K) + rsil*Q6sP(K)
!!
!!      jMI(3) = jMI(3) + rsil*Q6sP(K)
!!
!!      SD9m(K) = SD9m(K) + D1m(K)*rsil
!!------------------------------------------
!! Nutrient profiles and surface gradients
!!------------------------------------------
!
!! Preparing the profile:
!      diff1  = self%EDZ_1X *irrenh
!      diff2  = self%EDZ_2X *irrenh
!      diff3  = self%EDZ_3X *irrenh
!      cmix   = 0._rk
!      d1  = D1m
!      d2  = D2m
!      d3  = self%d_totX
!
!! Volume factors:
!      poro=self%qPWX
!      M1ads=self%M1adsX
!
!      !IF (N_COMP-N_UPPERX.GT.0) THEN
!      !   poro=benthic_morfology(1,K) ! porosity factor
!      !   M1ads=benthic_morfology(2,K) ! adsorption factor
!      !ENDIF
!
!      vG2  = poro
!      vM3  = poro
!      vM4  = poro*self%M4adsX
!      vM1  = poro*M1ads
!      vM11 = poro*self%M11adsX
!      vM5  = poro
!      vG3  = poro
!
!! profile parameters:
!
!      jmun = jmu(1)
!      IF (jmun<0._rk) THEN
!         ! To prevent pelagic ammonium in layer 1 from becoming negative under high nitrification rates,
!         ! fulfil part of the ammonium demand in nitrification by taking ammonium from sediment layer 2.
!         jmu(1) = jmun*N4n/(N4n+0.5_rk)
!         jmi(1) = jmi(1) + (jmun - jmu(1))
!      ENDIF
!      jmu(5) = 0.0_rk !TODO: Add flux of benthic carbonate SG3c(K) here
!      jmi(5) = 0._rk
!
!      CALL Prof_Parameter(profO2,  jMU(4), jMI(4), 0._rk, vG2,vG2,vG2)
!      CALL Prof_Parameter(profNO3, jMU(6),jMI(6),0._rk, vM3,vM3,vM3)
!      CALL Prof_Parameter(profN,   jMU(1),  jMI(1),  0._rk, vM4,vM4,vM4)
!      CALL Prof_Parameter(profP,   jMU(2),  jMI(2),  0._rk, vM1,vM1,vM11)
!      CALL Prof_Parameter(profS,   jMU(3),  jMI(3),  0._rk, vM5,vM5,vM5)
!      CALL Prof_Parameter(profCO2, jMU(5),jMI(5),0._rk, vG3,vG3,vG3)
!      CALL Prof_Parameter(profD,   d1,d2-d1,d3-d2, 1._rk,1._rk,1._rk)
!
!      CALL EquProfile(0._rk,1._rk,profD,dum,cmix,d1,d2,d3,diff1,diff2,diff3)
!
!      cmix   = self%EDZ_mixX

!      ! -----------------------------------------------------------------------------------
!      ! Oxygen
!      ! -----------------------------------------------------------------------------------
!      ! Layer 1: compute steady-state layer height H1_eq and layer integral c_int1_eq
!      call compute_final_equilibrium_profile(diff1,modconc(O2oP,jMU(4),cmix),jMU(4),self%d_totX,H1_eq,c_int1_eq)
!
!      !CALL EndProfile(O2oP,G2oP,profO2,jMN(4),cmix,d1,d2,d3,diff1,diff2,diff3)
!      !jMN(4) = jMN(4) - profO2(14)/self%relax_oX
!      !_SET_BOTTOM_ODE_(self%id_D1m,(max(D1m0,profO2(15))-D1m)/self%relax_oX)
!
!      ! Benthic oxygen dynamics: relax towards equilibrium value
!      c_int_eq = c_int1_eq*poro
!      _SET_BOTTOM_ODE_(self%id_G2o,(c_int_eq-G2oP)/self%relax_oX)
!
!      ! Net change in benthos must equal local production - surface exchange.
!      ! Thus, surface exchange = local production - net change (net change = relaxation)
!      _SET_BOTTOM_EXCHANGE_(self%id_O2o,jMU(4)-(c_int_eq-G2oP)/self%relax_oX)
!
!      ! Relax depth of first/oxygenated layer towards equilibrium value (H1_eq)
!      _SET_BOTTOM_ODE_(self%id_D1m,(max(D1m0,H1_eq)-D1m)/self%relax_oX)
!
!      ! -----------------------------------------------------------------------------------
!      ! Nitrate
!      ! -----------------------------------------------------------------------------------
!      ! Layer 1: compute steady-state concentration at bottom interface c_bot1_eq and layer integral c_int1_eq
!      call compute_equilibrium_profile(d1,diff1,modconc(N3nP,jMU(6)+jMI(6),cmix),jMU(6),jMI(6),c_bot1_eq,c_int1_eq)
!      ! Layer 2: compute steady-state layer height H2_eq and layer integral c_int2_eq
!      call compute_final_equilibrium_profile(diff2,c_bot1_eq,jMI(4),self%d_totX-d1,H2_eq,c_int2_eq)
!
!      !CALL EndProfile(N3nP,K3nP,profNO3,jMN(6),cmix,d1,d2,d3,diff1,diff2,diff3)
!      !jMN(6) = jMN(6) - profNO3(14)/self%relax_mX
!      !_SET_BOTTOM_ODE_(self%id_D2m,(max(D2m0,profNO3(15))-D2m)/self%relax_mX)
!
!      ! Benthic nitrate dynamics: relax towards equilibrium value
!      c_int_eq = (c_int1_eq+c_int2_eq)*poro
!      _SET_BOTTOM_ODE_(self%id_K3n,(c_int_eq-K3nP)/self%relax_mX)
!
!      ! Net change in benthos must equal local production - surface exchange.
!      ! Thus, surface exchange = local production - net change (net change = relaxation)
!      _SET_BOTTOM_EXCHANGE_(self%id_N3n,jMU(6)+jMI(6)-(c_int_eq-K3nP)/self%relax_mX)
!      
!      ! Relax depth of bottom interface of second/oxidised layer towards equilibrium value (d1+H2_eq)
!      _SET_BOTTOM_ODE_(self%id_D2m,(max(D2m0,d1+H2_eq)-D2m)/self%relax_mX)
!
!      ! -----------------------------------------------------------------------------------
!      ! Ammonium
!      ! -----------------------------------------------------------------------------------
!
!     !!!GL CALL EquProfile(N4nP,K4nP,profN,jMN(1),cmix,d1,d2,d3,diff1,diff2,diff3)
!      ! Non-equilibrium correction:
!     !!!GL CALL NonEquFlux(profN,profD,jMN(1))
!     !!!GL jMN(1) = profN(1) + profN(2) + profN(3)
!
!      ! Layer 1: compute steady-state concentration at bottom interface c_bot1_eq and layer integral c_int1_eq
!     call compute_equilibrium_profile(d1,diff1,modconc(N4nP,jMU(1)+jMI(1),cmix),jMU(1),jMI(1),c_bot1_eq,c_int1_eq)
!      ! Layer 2: compute steady-state concentration at bottom interface c_bot2_eq and layer integral c_int2_eq
!     call compute_equilibrium_profile(d2-d1,diff2,c_bot1_eq,jMI(1),0.0_rk,c_bot2_eq,c_int2_eq)
!      ! Layer 3: no sources or sinks: homogeneous equilibrium concentration c_bot2_eq
!     c_int3_eq = (d3-d2)*c_bot2_eq
!     c_int_eq = poro*self%M4adsX*(c_int1_eq+c_int2_eq+c_int3_eq)
!
!      ! The equilibrium depth-integrated mass c_int_eq usually differs from the current depth-integrated mass K4nP.
!      ! We can view the actual [unknown] pore water concentration profile as the sum of the equilibrium profile
!      ! and a residual profile. The latter has a vertical integral equal to the difference between actual mass and equilibrium mass.
!      ! The rate at which the residual mass is exchanged over the benthic-pelagic interface is equal to the product of
!      ! the diffusivity and gradient in the residual mass at this interface. We do not know this gradient since we do not
!      ! know the shape of the residual profile. Let's make some simple assumptions to infer this.
!      ! Constraints: diffusion of the residual across bottom of benthic column must be zero (i.e., zero gradient), and at the surface of the benthic
!      ! column the concentration of the residual must equal zero (i.e., equilibrium holds at the very surface of the column).
!      ! Since we do not know anything about the processes responsible for the residual, let's assume their contribution
!      ! in the past was a constant production or destruction per unit sediment volume thoughout the entire column.
!      ! That is, production (#/m^2/d) in the three layers was P_int*d1/d3, P_int(d2-d1)/d3, P_int(d3-d2)/d3.
!      ! If we would know P_int, we could supply those rates along with zero surface concentration to "compute_equilibrium_profile"
!      ! to derive the residual profile. By checking the equations in compute_equilibrium_profile, we can verify that the resulting bottom concentration
!      ! and layer integral are both proportional to P_int. Thus, can can simply supply d1, d2-d1, d3-d2 to
!      ! "compute_equilibrium_profile", and find the additional scale factor P_int/d3 by demanding that the sum of layer integrals is
!      ! equal to the known residual mass. That is, P_int/d3 equals the ratio of residual mass to the sum of normalized layer integrals
!      ! computed for layer production terms d1, d2-d1, d3-d2. As we are assuming the residual profile was previously an equilibrium
!      ! profile, the necessary depth-integrated production rate P_int must equal the exchange across the surface, i.e., diffusivity*gradient.
!      ! Thus, we can now simply add the P_int as a additional surface exchange term, accounting for the move towards equilibrium.
!      call compute_equilibrium_profile(d1,   diff1,0.0_rk,   d1,   d3-d1, c_bot1_eq,c_int1_eq)
!      call compute_equilibrium_profile(d2-d1,diff2,c_bot1_eq,d2-d1,d3-d2, c_bot2_eq,c_int2_eq)
!      call compute_equilibrium_profile(d3-d2,diff3,c_bot2_eq,d3-d2,0.0_rk,c_bot3_eq,c_int3_eq)
!      norm_res_int = poro*self%M4adsX*(c_int1_eq+c_int2_eq+c_int3_eq)
!      P_res_int = (K4nP-c_int_eq)/norm_res_int*d3
!      _SET_BOTTOM_EXCHANGE_(self%id_N4n,jMU(1)+jMI(1)+P_res_int) ! Equilibrium flux = jMU(1)+jMI(1), residual flux = P_res_int
!      _SET_BOTTOM_ODE_(self%id_K4n,-P_res_int)
!
!!------------------------------------------------------------------------------
!!! phospate
!!------------------------------------------------------------------------------
!!
!!      CALL EquProfile(N1pP(I),K1pP(K),profP,jMN(2),k)
!!      ! Non-equilibrium correction:
!!      CALL NonEquFlux(profP,profD,jMN(2))
!!      jMN(2) = profP(1) + profP(2) + profP(3)
!     ! Layer 1: compute steady-state concentration at bottom interface c_bot1_eq
!     ! and layer integral c_int1_eq
!     call compute_equilibrium_profile(d1,diff1,modconc(N1pP,jMU(2)+jMI(2),cmix),jMU(2),jMI(2),c_bot1_eq,c_int1_eq)
!      ! Layer 2: compute steady-state concentration at bottom interface
!      ! c_bot2_eq and layer integral c_int2_eq
!     call compute_equilibrium_profile(d2-d1,diff2,c_bot1_eq,jMI(2),0.0_rk,c_bot2_eq,c_int2_eq)
!      ! Layer 3: no sources or sinks: homogeneous equilibrium concentration
!      ! c_bot2_eq
!     c_int3_eq = (d3-d2)*c_bot2_eq
!     c_int_eq = poro*(self%M1adsX*(c_int1_eq+c_int2_eq)+self%M11adsX*c_int3_eq)
!     call compute_equilibrium_profile(d1,   diff1,0.0_rk,   d1,   d3-d1, c_bot1_eq,c_int1_eq)
!     call compute_equilibrium_profile(d2-d1,diff2,c_bot1_eq,d2-d1,d3-d2, c_bot2_eq,c_int2_eq)
!     call compute_equilibrium_profile(d3-d2,diff3,c_bot2_eq,d3-d2,0.0_rk,c_bot3_eq,c_int3_eq)
!      norm_res_int = poro*(self%M1adsX*(c_int1_eq+c_int2_eq)+self%M11adsX*c_int3_eq)
!      P_res_int = (K1pP-c_int_eq)/norm_res_int*d3
!      _SET_BOTTOM_EXCHANGE_(self%id_N1p,jMU(2)+jMI(2)+P_res_int) ! Equilibrium flux = jMU(2)+jMI(2), residual flux = P_res_int
!      _SET_BOTTOM_ODE_(self%id_K1p,-P_res_int)
!
!!
!!      endif
!!#endif
!!
!!! silicate
!!      
!!      CALL EquProfile(N5sP(I),K5sP(K),profS,jMN(3),k)
!!      ! Non-equilibrium correction:
!!      CALL NonEquFlux(profS,profD,jMN(3))
!!      jMN(3) = profS(1) + profS(2) + profS(3)
!!
!!-------------------------------------------------------------------------
!!! carbon dioxide
!!-------------------------------------------------------------------------
!!     
!!      CALL EquProfile(O3cP(I),G3cP(K),profCO2,jMN(5),k)
!!      ! Non-equilibrium correction:
!!      CALL NonEquFlux(profCO2,profD,jMN(5))
!!      jMN(5) = profCO2(1) + profCO2(2) + profCO2(3)
!
! ! Layer 1: compute steady-state concentration at bottom interface c_bot1_eq
!     ! and layer integral c_int1_eq
!     call compute_equilibrium_profile(d1,diff1,modconc(O3cP,jMU(5)+jMI(5),cmix),jMU(5),jMI(5),c_bot1_eq,c_int1_eq)
!      ! Layer 2: compute steady-state concentration at bottom interface
!      ! c_bot2_eq and layer integral c_int2_eq
!     call compute_equilibrium_profile(d2-d1,diff2,c_bot1_eq,jMI(5),0.0_rk,c_bot2_eq,c_int2_eq)
!      ! Layer 3: no sources or sinks: homogeneous equilibrium concentration
!      ! c_bot2_eq
!     c_int3_eq = (d3-d2)*c_bot2_eq
!     c_int_eq = poro*(c_int1_eq+c_int2_eq+c_int3_eq)
!     call compute_equilibrium_profile(d1,   diff1,0.0_rk,   d1,   d3-d1, c_bot1_eq,c_int1_eq)
!     call compute_equilibrium_profile(d2-d1,diff2,c_bot1_eq,d2-d1,d3-d2, c_bot2_eq,c_int2_eq)
!     call compute_equilibrium_profile(d3-d2,diff3,c_bot2_eq,d3-d2,0.0_rk,c_bot3_eq,c_int3_eq)
!      norm_res_int = poro*(c_int1_eq+c_int2_eq+c_int3_eq)
!      P_res_int = (G3cP-c_int_eq)/norm_res_int*d3
!      _SET_BOTTOM_EXCHANGE_(self%id_O3c,jMU(5)+jMI(5)+P_res_int) ! Equilibrium flux = jMU(2)+jMI(2), residual flux = P_res_int
!      _SET_BOTTOM_ODE_(self%id_G3c,-P_res_int)

!
!! nitrate gas:
!      
!      jmn(7) = 0._rk
!
!!-----------
!! Diffusion
!!-----------
!
      !_SET_BOTTOM_ODE_(self%id_K4n, - jMN(1))
      !_SET_BOTTOM_EXCHANGE_(self%id_N4n, + jMN(1))
      !_SET_BOTTOM_ODE_(self%id_K1p, - jMN(2))
      !_SET_BOTTOM_EXCHANGE_(self%id_N1p, + jMN(2))
      !_SET_BOTTOM_ODE_(self%id_K5s, - jMN(3))
      !_SET_BOTTOM_EXCHANGE_(self%id_N5s, + jMN(3))
      !_SET_BOTTOM_ODE_(self%id_G2o, - jMN(4))
      !_SET_BOTTOM_EXCHANGE_(self%id_O2o, + jMN(4))
      !_SET_BOTTOM_ODE_(self%id_G3c, - jMN(5))
      !_SET_BOTTOM_EXCHANGE_(self%id_O3c, + jMN(5))
      !_SET_BOTTOM_ODE_(self%id_K3n, - jMN(6))
      !_SET_BOTTOM_EXCHANGE_(self%id_N3n,jMN(6))
      !_SET_BOTTOM_ODE_(self%id_G4n, - jMN(7))

      _HORIZONTAL_LOOP_END_

   end subroutine

   subroutine compute_equilibrium_profile(D,sigma,c0,P,flux_bot,c_bot,c_int)
      real(rk),intent(in)  :: D,sigma,c0,P,flux_bot
      real(rk),intent(out) :: c_bot,c_int
      real(rk) :: a,b,c
      ! ----------------------------------------------------------------------------------------------------
      ! Determine equilibrium concentration profile in pore water from:
      ! - D:        layer thickness (m)
      ! - sigma:    diffusivity (m2/d)
      ! - c0:       concentration at layer surface (#/m3)
      ! - P:        layer-integrated source-sink terms (#/m2/d)
      ! - flux_bot: bottom flux (#/m2/d)
      ! Returns:
      ! - c_bot:    equilibrium concentration at bottom interface (#/m3)
      ! - c_int:    depth integral of concentration (#/m2)
      ! ----------------------------------------------------------------------------------------------------
      ! Governing equation: dy/dt = sigma d2y/dz2 + sms
      ! Assumption: diffusivity sigma and sources-minus-sinks sms are independent of z within the layer.
      ! Thus, sms can be written as layer integrated source-sink terms, divided by layer height: sms = P/D
      ! At equilibrium: sigma d^2y/dz^2 + P/D = 0
      ! Thus, d^2y/dz^2 = -P/D/sigma
      ! Solution is a quadratic equation: c(z) = a z^2 + b z + c
      ! From second derivative: d^2y/dz^2 = 2a = -P/D/sigma. Thus, a = -P/D/sigma/2.
      ! ----------------------------------------------------------------------------------------------------
      ! c(z) = -P/D/sigma/2 z^2 + b z + c
      !
      ! Lets adopt a bottom-to-top coordinate system, with the bottom interface at z=0 and the top interface at z=D.
      !
      ! Constraint 1: flux over bottom interface is known: flux_bot
      ! (typically chosen to balance demand depth-integrated sinks-sources in deeper layers)
      ! For consistency, this flux must equal that produced by diffusion at the boundary, i.e., sigma*dy/dz
      ! Note: gradient (bottom to top!) must be positive when deeper layers are a sink (i.e., sms_int_deep<0)
      ! sigma*dy/dz(0) = sigma*b = -flux_bot -> b = -flux_bot/sigma
      !
      ! Constraint 2: top concentration c(D)=c0 is known.
      ! c(D) = -P/sigma/2 D - flux_bot/sigma D + c = c0
      ! -> c = c0 + (P/2 + flux_bot)/sigma D
      ! This is also the concentration at bottom interface: c(0) = c
      !
      ! Depth-integrated layer contents is found by integrating the parabola between 0 and D:
      !    \int{a * z^2 + b z + c} = [a/3 z^3 + b/2 z^2 + c z]_0^D = a/3 D^3 + b/2 D^2 + c Z
      ! ----------------------------------------------------------------------------------------------------
      a = -P/D/sigma/2
      b = -flux_bot/sigma
      c = c0 + D*(P/2 + flux_bot)/sigma

      c_bot = c
      c_int = (a/3*D*D + b/2*D + c)*D
   end subroutine

   subroutine compute_final_equilibrium_profile(sigma,c0,P,Dmax,D,c_int)
      real(rk),intent(in)  :: sigma,c0,P,Dmax
      real(rk),intent(out) :: D,c_int
      real(rk) :: c_bot
      ! ----------------------------------------------------------------------------------------------------
      ! Determine layer depth and concentration profile in pore water at equilibrium from:
      ! - sigma:    diffusivity (m2/d)
      ! - c0:       concentration at layer surface (#/m3)
      ! - P:        layer-integrated source-sink terms (#/m2/d)
      ! Constraint: concentration and flux at bottom interface must be zero.
      ! Returns:
      ! - D:        layer thickness (m)
      ! - c_int:    depth integral of concentration (#/m2)
      ! ----------------------------------------------------------------------------------------------------
      ! Governing equation in 1st layer: dy/dt = sigma d2y/dz2 + sms
      ! Assumption: diffusivity sigma and sink-minus-source sms are independent of z within the layer.
      ! Thus, sms can be written as layer integral divided by layer height: sms = P/D
      ! At equilibrium: sigma d^2y/dz^2 + P/D = 0
      ! Thus, d^2y/dz^2 = -P/D/sigma
      ! Solution is a quadratic equation: c(z) = a z^2 + b z + c
      ! From second derivative: d^2y/dz^2 = 2a = -P/D/sigma. Thus, a = -P/D/sigma/2.
      ! ----------------------------------------------------------------------------------------------------
      ! Task 1: find layer height
      !   Bottom-to-top coordinate system, bottom interface at z=0, top interface at z=D
      !   Bottom constraint: concentration c(0) = 0, flux dy/dz = 0 (i.e., minimum of parabola)
      !   Thus, b=0 and c=0, and c(z) = -P/D/sigma/2 z^2
      !   Top constraint: concentration c(D) is known: c0
      !   Thus c(D) = -P/sigma/2 D = c0 -> layer height D = -2 sigma c0/P
      !   Note that D>=0 for P<0 only: this solution is valid only if the layer is a sink.
      !   Also, if the layer is neither sink nor source (P=0), D->infinity. That is,
      !   the solution is a flat profile, with surface concentration c0 extending forever downward.
      ! Task 2: compute depth-integrated layer contents:
      !   \int{-P/D/sigma/2 z^2} = [-P/D/sigma/6 z^3]_0^D = -P/sigma/6 D^2
      ! Check consistency: flux at top interface is sigma*dy/dz = -P
      !    [OK: in steady state, surface exchange compensates internal loss]
      ! ----------------------------------------------------------------------------------------------------
      if (Dmax*P>-2*sigma*c0) then   ! Dmax<-2 sigma c0/P, rearranged for P<0 close to 0. Result also picks up P>0.
         ! Loss rate within layer is too low (or layer experiences net production, i.e., P>0).
         ! If we would impose zero concentration at the layer bottom, the layer would extend beyond maximum depth.
         ! Fix depth at maximum depth and use parabola with zero flux but non-zero concentration at bottom interface.
         D = Dmax
         call compute_equilibrium_profile(D,sigma,c0,P,0.0_rk,c_bot,c_int)
      else
         D = -2*sigma*c0/P
         c_int = -P/sigma/6*D*D
      end if
   end subroutine

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Prof_Parameter \label{sec:profParameter}
!
! !DESCRIPTION:
!  TODO - CHECK THIS
!
!  Load layer productions and volume factors into profile.
!\\
!\\
! !INTERFACE:
   SUBROUTINE Prof_Parameter(prof,s1,s2,s3,v1,v2,v3)
!
! !INPUT PARAMETERS:
! ! TODO - document these.
  real(rk),intent(in) :: s1,s2,s3,v1,v2,v3
!
! !OUTPUT PARAMETERS:
! ! TODO - document these
  real(rk),intent(out) :: prof(15)
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
!  code................................................

      prof(1) = s1  !layer fluxes
      prof(2) = s2
      prof(3) = s3
      prof(8) = v1 !layer volume factors
      prof(9) = v2
      prof(10) = v3

      RETURN
      END SUBROUTINE Prof_Parameter
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: EquProfile \label{sec:equProfile}
!
! !DESCRIPTION:
!  TODO - CHECK THIS.
!
!  computes diffusion-production equilibrium profiles
!    from layer productions, layer depths and benthic/pelagic interface 
!    concentration.
! 
!    a\_water: concentration in water box above the seafloor (per $m^3$)
!    mass: conent in benthic section (per $m^2$)
!    p: working array with
!       1: first layer production
!       2: second layer production
!       3: third layer production
!       4: concentration at seafloor (benthic/pelagic interface)
!       5: concentration at 1st/2nd layer interface
!       6: concentration at 2nd/3rd layer interface
!       7: concentraion at lower boundary of 3rd layer
!       8: 1st layer volume factor (porosity)
!       9: 2nd layer volume factor (porosity)
!      10: 3rd layer volume factor (porosity)
!      11: 1st layer content (per $m^2$)
!      12: 2nd layer content (per $m^2$)
!      13: 3rd layer content (per $m^2$)
!      14: computed difference in total benthic layer content
!      15: lower boundary of 3rd layer
!    f\_0: flux from benthic to pelagic layer 
!\\
!\\
! !INTERFACE:
   SUBROUTINE EquProfile(a_water,mass,p,f0,cmix,d1,d2,d3,diff1,diff2,diff3)
!
! !INPUT PARAMETERS:
   real(rk),intent(in) :: a_water,mass,cmix,d1,d2,d3,diff1,diff2,diff3
!
! !INPUT/OUTPUT PARAMETERS:
!  ! TODO - document these.
   real(rk) :: f0,p(15)
!
! !LOCAL VARIABLES:
!  ! TODO - Document these.
   real(rk) :: f1,f2
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
! fluxes from Sources:
      f0 = p(1)+ p(2)+ p(3)
      f1 =       p(2)+ p(3)
      f2 =             p(3)

! Surface concentration a0:
      p(4) = modconc(a_water,f0,cmix)

! Profile in form of three parabular pieces:
      CALL arc   (p(4),p(5), f0,f1, diff1, 0._rk,d1, p(11))
      CALL arc   (p(5),p(6), f1,f2, diff2, d1,d2, p(12))
      CALL arc   (p(6),p(7), f2,0._rk, diff3, d2,d3, p(13)) ! 0 bottom flux
      p(15) = d3

! Difference of equilibrium masses to real mass:
      p(14) = p(11)*p(8) + p(12)*p(9) + p(13)*p(10) - mass

   END SUBROUTINE EquProfile

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: modconc \label{sec:modconc}
!
! !DESCRIPTION:
!  TODO - CHECK THIS
!
!  Modification of surface concentration by the diffusion flux
!\\
!\\
! !INTERFACE:
   real(rk) FUNCTION modconc(conc,flux,cmix)
!
! !LOCAL VARIABLES:
!  ! TODO - document these.
   real(rk), intent(in) :: conc,flux,cmix
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
! conc mg/m3, flux mg/m2day, cmix day/m

      IF (flux.GE.0._rk) THEN
         modconc = conc + cmix*flux
      ELSE
         modconc = conc*conc/(conc - cmix*flux)
      ENDIF

      RETURN

   END FUNCTION modconc
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: arc \label{sec:arc}
!
! !DESCRIPTION:
!  TODO - CHECK THIS
!
!  Computes for a layer the lower face concentration and the layer content 
!  from upper face concentration and lower and upper face fluxes, assuming
!  equilibrium between inner production and diffusion.
!  Solves: \verb+a''=(f1-f0)/(diff*(d1-d0)+ 
!     with \verb+diff*a'=f0 at face 0 and diff*a'=f1\verb+ at face 1
!    \verb+=> a(d,t) = (f1-f0)/(diff*(d1-d0))*.5*(d-d0)^2+f0/diff*(d-d0)+a0+
!
!  a0: upper face concentration (per $m^3$)
!  a1: lower face concentration (per $m^3$)
!  f0: upper flace flux (outwards)
!  f1: lower face flux (inwards)
!  d0: depth level upper face
!  d1: depth level lower face
!  m1: total benthic content (per $m^2$)
!\\
!\\
! !INTERFACE:
   SUBROUTINE arc(a0,a1,f0,f1,diff,d0,d1,m1)
!
! !INPUT PARAMETERS:
! ! TODO - document these
  real(rk),intent(in) :: a0,d0,d1,f0,f1,diff
!
! !INPUT/OUTPUT PARAMETERS:
! ! TODO - document these.
  real(rk),intent(out) :: a1,m1
!
! !LOCAL VARIABLES:
! ! TODO - document these.
  real(rk) :: g0,g1
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
!  code................................................

      g0 = f0/diff
      g1 = f1/diff

      a1 = a0 + (g0+g1)*(d1-d0)/2._rk
      m1 = (a0 + (g0+g1/2._rk)*(d1-d0)/3._rk)*(d1-d0)

      RETURN
  END SUBROUTINE arc
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: endarc \label{sec:endarc}
!
! !DESCRIPTION:
!  TODO - CHECK THIS
!
!  For profiles with negative concentrations at the bottom interface cut 
!  the layer at the intersection. Works only for f1<=0.
!  a0: upper face concentration (per $m^3$)
!  a1: lower face concentration (per $m^3$)
!  f0: upper flace flux (outwards)
!  f1: lower face flux (inwards)
!  d0: depth level upper face
!  d1: depth level lower face
!  m1: total benthic content (per $m^2$)
!\\
!\\
! !INTERFACE:
  SUBROUTINE endarc(a0,a1,f0,f1,diff,d0,d1,dt,m1)
!
! !INPUT PARAMETERS:
!     ! TODO - document these.
      real(rk),intent(in) :: diff,f0,f1,a0,d0,dt
!
! !OUTPUT PARAMETERS:
!     ! TODO - document these.
      real(rk),intent(out) :: a1,d1,m1
!
! !LOCAL VARIABLES:
!     ! TODO - document these.
      real(rk) :: df,pf,vlow
      DATA vlow / 1.e-10_rk /
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOC
!
!  code................................................

         pf = f0 + f1


         if ( abs(pf).gt.vlow) then
            !df = (f0 + 2._rk*f1)/pf
            df = 1._rk + f1/pf
         else
            pf = vlow * SIGN(1._rk,pf)
            if ((abs(f0).lt.vlow).AND.(abs(f1).lt.vlow)) then
               df = 1.5_rk
            else
               if ( abs(f0).lt.vlow ) df = 2._rk
               if ( abs(f1).lt.vlow ) df = 1._rk
            endif
         endif

         d1 = d0 - 2._rk*a0/pf*diff

         IF (d1.GE.d0.AND.d1.LE.dt) THEN ! cut profile
            m1 = a0*(d1-d0)*df/3._rk
            a1 = 0._rk
         ELSE ! standard case
            CALL arc(a0,a1,f0,f1,diff,d0,dt,m1)
            d1=dt
         ENDIF

      RETURN
   END SUBROUTINE endarc
!
!EOC
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: EndProfile \label{sec:endProfile}
!
! !DESCRIPTION:
!  TODO - CHECK THIS
!  computes diffusion-production equilibrium profiles
!    from layer productions, layer depths and benthic/pelagic interface 
!    concentration. Clips negativity.
! 
!    a\_water: concentration in water box above the seafloor (per $m^3$)
!    mass: conent in benthic section (per $m^2$)
!    p: working array with
!       1: first layer production
!       2: second layer production
!       3: third layer production
!       4: concentration at seafloor (benthic/pelagic interface)
!       5: concentration at 1st/2nd layer interface
!       6: concentration at 2nd/3rd layer interface
!       7: concentraion at lower boundary of 3rd layer
!       8: 1st layer volume factor (porosity)
!       9: 2nd layer volume factor (porosity)
!      10: 3rd layer volume factor (porosity)
!      11: 1st layer content (per $m^2$)
!      12: 2nd layer content (per $m^2$)
!      13: 3rd layer content (per $m^2$)
!      14: computed difference in total benthic layer content
!      15: lower boundary of 3rd layer
!    f\_0: flux from benthic to pelagic layer
!\\
!\\
! !INTERFACE:
   SUBROUTINE EndProfile(a_water,mass,p,f0,cmix,d1,d2,d3,diff1,diff2,diff3)
!
! !INPUT PARAMETERS:
! ! TODO - document these
  real(rk),intent(in) :: a_water,mass,cmix,d1,d2,d3,diff1,diff2,diff3
!
! !INPUT/OUTPUT PARAMETERS:
! ! TODO - document this.
  real(rk) :: p(15),f0
!
! !LOCAL VARIABLES:
! ! TODO - document these
  real(rk) :: f1,f2
  LOGICAL case1,case2,case3,case4
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
! fluxes from Sources:
      f0 = p(1)+ p(2)+ p(3)
      f1 =       p(2)+ p(3)
      f2 =             p(3)

! Surface concentration a0:
       p(4) = modconc(a_water,f0,cmix)
! Case differentiation:
! case1: oxigen
! case2: nitrate
! case4: ammonium,phosphate,silicate,co2

      case1 = ((p(1).LT.0._rk) .AND. (p(2).LE.0._rk) .AND. (p(3).EQ.0._rk))
      case2 = ((p(1).GE.0._rk) .AND. (p(2).LT.0._rk) .AND. (p(3).LE.0._rk))
      case3 = ((p(1).GE.0._rk) .AND. (p(2).GE.0._rk) .AND. (p(3).LT.0._rk))
      case4 = ((p(1).GE.0._rk) .AND. (p(2).GE.0._rk) .AND. (p(3).EQ.0._rk))

! Profile in form of up to three parabular pieces:

      IF (case1) THEN
         CALL endarc(p(4),p(5), f0,f1, diff1, 0._rk,p(15),d3,p(11))
         p(12) = 0._rk
         p(13) = 0._rk
      ELSEIF (case2) THEN
         CALL arc   (p(4),p(5), f0,f1, diff1, 0._rk,d1,      p(11))
         CALL endarc(p(5),p(6), f1,f2, diff2, d1,p(15),d3, p(12))
         p(13) = 0._rk
      ELSEIF (case3) THEN
         CALL arc   (p(4),p(5), f0,f1, diff1, 0._rk,d1,      p(11))
         CALL arc   (p(5),p(6), f1,f2, diff2, d1,d2,      p(12))
         CALL endarc(p(6),p(7), f2,0._rk, diff3, d2,p(15),d3, p(13))
      ELSEIF (case4) THEN
         CALL arc   (p(4),p(5), f0,f1, diff1, 0._rk,d1,      p(11))
         CALL arc   (p(5),p(6), f1,f2, diff2, d1,d2,      p(12))
         CALL arc   (p(6),p(7), f2,0._rk, diff3, d2,d3,      p(13))
         p(15) = d3
      ELSE
         CALL arc   (p(4),p(5), f0,f1, diff1, 0._rk,d1,      p(11))
         CALL arc   (p(5),p(6), f1,f2, diff2, d1,d2,      p(12))
         CALL arc   (p(6),p(7), f2,0._rk, diff3, d2,d3,      p(13))
         p(15) = d3
      ENDIF

      p(14) = p(11)*p(8) + p(12)*p(9) + p(13)*p(10) - mass

      RETURN

      END SUBROUTINE EndProfile
!
!EOC
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: NonEquFlux \label{sec:nonEquFlux}
!
! !DESCRIPTION:
!  TODO - CHECK THIS
!
!  Adds non-equilibrium modification to equilibrium flux.  
!              assumes homogenous distribution of flux correction
!   prof:    (1-3): layer fluxes
!            (8-10): layer volume-factors (porosity+adsorption)
!            (14): Difference of eq. content to real content
!            (15): flux to pelagic layer
!   profD:   (1-3): layer thicknesses
!            (11:13): layer equilibrium contents (per $m^2$)
!   j\_diff:  outward flux on benthic surface 
!\\
!\\
! !INTERFACE:
   SUBROUTINE NonEquFlux(prof,profD,j_diff)
!
! !OUTPUT PARAMETERS:
!  ! TODO - document this.
   real(rk),intent(out) :: j_diff
!
! !INPUT/OUTPUT PARAMETERS:
!  ! TODO - docuemnt these
   real(rk) :: prof(15),profD(15)
!
! !LOCAL VARIABLES:
!  ! TODO - docuemnt these.
   real(rk) :: md,factor
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
!  code.................................................................
!
!  projects content difference p(14) on standard parabolic with 
!    0 surface concentration and 0 bottom flux
      md = prof(8)*profD(11) + prof(9)*profD(12) +prof(10)*profD(13)
      factor = prof(14)/md 

      prof(1) = prof(1) - profD(1)*factor
      prof(2) = prof(2) - profD(2)*factor
      prof(3) = prof(3) - profD(3)*factor

      j_diff = prof(1) + prof(2) + prof(3)

      RETURN

  END SUBROUTINE NonEquFlux
!
!EOC
!-----------------------------------------------------------------------

   subroutine dissolved_matter_per_layer_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !INPUT PARAMETERS:
   class (type_dissolved_matter_per_layer),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_

   real(rk) :: c_int
   real(rk) :: factor
   real(rk) :: d(3),D2m,d_totX,poro

   _HORIZONTAL_LOOP_BEGIN_

      ! Get depth-integrated concentration.
      _GET_HORIZONTAL_(self%id_tot,c_int)

      ! Layer dpeths and porosity.
      _GET_HORIZONTAL_(self%id_D1m,d(1))
      _GET_HORIZONTAL_(self%id_D2m,D2m)
      _GET_HORIZONTAL_(self%id_Dtot,d_totX)
      _GET_HORIZONTAL_(self%id_poro,poro)

      ! Compute layer ticknesses from depth of bottom interface of individual layers.
      d(2) = D2m-d(1)
      d(3) = d_totX-D2m

      if (self%type==1) then
         ! All in layer 1
         ! Typically used for oxygen
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers(1),c_int)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers(2),0.0_rk)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers(3),0.0_rk)
      elseif (self%type==2) then
         ! Vertically homogeneous in layer 1, quadratically decreasing in layer 2 (zero concentration at bottom interface)
         ! Typically used for nitrate

         ! Mean concentration in pore water (matter/m3)
         factor = c_int/(poro*(self%ads(1)*d(1)+self%ads(2)*d(2)/3.0_rk))

         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers(1),poro*self%ads(1)*d(1)*factor)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers(2),poro*self%ads(2)*d(2)/3.0_rk*factor)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers(3),0.0_rk)
      else
         ! Vertically homogeneous in layers 1,2,3
         ! Typically used for all tracers but oxygen and nitrate.

         ! Mean concentration in pore water (matter/m3)
         factor = c_int/sum(poro*self%ads*d)

         ! Layer contents (matter/m2)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers(1),poro*self%ads(1)*d(1)*factor)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers(2),poro*self%ads(2)*d(2)*factor)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layers(3),poro*self%ads(3)*d(3)*factor)
      end if
      
   _HORIZONTAL_LOOP_END_
   
   end subroutine dissolved_matter_per_layer_do_bottom
   
end module
