.. _api::

#########
ERSEM api
#########



benthic_carbonate module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine do_bottom(self,e__arguments_do_bottom_e) 

light_iop module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine get_light(self,e__arguments_vertical_e) 

light_iop_ady module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine do_bn(self,e__arguments_do_e) 
.. code-block:: bash 

      subroutine get_light(self,e__arguments_vertical_e) 

bacteria_docdyn module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine do_bn(self,e__arguments_do_e) 
.. code-block:: bash 

      subroutine remineralization(self,e__arguments_do_e) 

oxygen module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine do_bn(self,e__arguments_do_e) 
.. code-block:: bash 

      subroutine do_surface(self,e__arguments_do_surface_e) 
.. code-block:: bash 

      function oxygen_saturation_concentration(self,etw,x1x) result (osat) 
          real(kind=rk), intent(in) :: etw
          real(kind=rk), intent(in) :: x1x
          real(kind=rk) :: osat

light module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine get_light(self,e__arguments_vertical_e) 

pelagic_base module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine initialize_ersem_base(self,rm,sedimentation) 
          real(kind=rk), optional,intent(in) :: rm
          logical, optional,intent(in) :: sedimentation
.. code-block:: bash 

      subroutine add_constituent(self,name,initial_value,background_value,qn,qp) 
          character(len=*), intent(in) :: name
          real(kind=rk), intent(in) :: initial_value
          real(kind=rk), optional,intent(in) :: background_value
          real(kind=rk), optional,intent(in) :: qn
          real(kind=rk), optional,intent(in) :: qp
          subroutine register_bn(variable_id,name,base_units,long_name,aggregate_variable,qx,id_xdep,id_dep,scale_factor) ! in ../../src/pelagic_base.F90:ersem_pelagic_base:add_constituent
              type(type_state_variable_id), intent(inout),target :: variable_id
              character(len=*), intent(in) :: name
              character(len=*), intent(in) :: base_units
              character(len=*), intent(in) :: long_name
              type(type_bulk_standard_variable), intent(in) :: aggregate_variable
              real(kind=rk), intent(inout),allocatable,optional,dimension(:) :: qx
              type(type_horizontal_diagnostic_variable_id), intent(inout),allocatable,optional,dimension(:) :: id_xdep
              type(type_bottom_state_variable_id), intent(inout),allocatable,optional,dimension(:) :: id_dep
              real(kind=rk), intent(in),optional :: scale_factor
              integer :: idep
              character(len=16) :: num
          end subroutine register_bn
.. code-block:: bash 

      function get_sinking_rate(self,e__arguments_local_e) result (rm) 
          real(kind=rk) :: rm
.. code-block:: bash 

      subroutine do_bottom(self,e__arguments_do_bottom_e) 

primary_producer module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine do_bn(self,e__arguments_do_e) 
.. code-block:: bash 

      function get_sinking_rate(self,e__arguments_local_e) result (sd) 
          real(kind=rk) :: sd
.. code-block:: bash 

      subroutine get_vertical_movement(self,e__arguments_get_vertical_movement_e) 

nitrification module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine do_bn(self,e__arguments_do_e) 

calcification module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine do_bn(self,e__arguments_do_e) 

ersem_model_library module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self) 
.. code-block:: bash 

      subroutine create(self,name,model) 
          character(len=*), intent(in) :: name

zenith_angle module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine do_surface(self,e__arguments_do_surface_e) 

benthic_base module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine initialize_ersem_benthic_base(self) 
.. code-block:: bash 

      subroutine benthic_base_add_constituent(self,name,initial_value,background_value,qn,qp) 
          character(len=*), intent(in) :: name
          real(kind=rk), intent(in) :: initial_value
          real(kind=rk), optional,intent(in) :: background_value
          real(kind=rk), optional,intent(in) :: qn
          real(kind=rk), optional,intent(in) :: qp
          subroutine register_bn(variable_id,resuspended_id,resuspended_flux_id,name,base_units,long_name,aggregate_variable,scale_factor) ! in ../../src/benthic_base.F90:ersem_benthic_base:benthic_base_add_constituent
              type(type_bottom_state_variable_id), intent(inout),target :: variable_id
              type(type_state_variable_id), intent(inout),target :: resuspended_id
              type(type_horizontal_diagnostic_variable_id), intent(inout),target :: resuspended_flux_id
              character(len=*), intent(in) :: name
              character(len=*), intent(in) :: base_units
              character(len=*), intent(in) :: long_name
              type(type_bulk_standard_variable), intent(in) :: aggregate_variable
              real(kind=rk), intent(in),optional :: scale_factor
          end subroutine register_bn
.. code-block:: bash 

      subroutine do_bottom(self,e__arguments_do_bottom_e) 

benthic_column module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine do_bottom(self,e__arguments_do_bottom_e) 
.. code-block:: bash 

      subroutine bioturbation_initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine bioturbation_do_bottom(self,e__arguments_do_bottom_e) 

benthic_bacteria module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine do_bottom(self,e__arguments_do_bottom_e) 

mesozooplankton module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine do_bn(self,e__arguments_do_e) 

bacteria module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine do_bn(self,e__arguments_do_e) 
.. code-block:: bash 

      subroutine remineralization(self,e__arguments_do_e) 

microzooplankton module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine do_bn(self,e__arguments_do_e) 

benthic_fauna module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine do_bottom(self,e__arguments_do_bottom_e) 

benthic_calcite module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine do_bottom(self,e__arguments_do_bottom_e) 

benthic_column_particulate_matter module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine do_bottom(self,e__arguments_do_bottom_e) 
.. code-block:: bash 

      subroutine layer_initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine layer_initialize_constituent(self,info,name,units,long_name,remin,q10,source_depth_distribution,aggregate_target,aggregate_scale_factor) 
          character(len=*), intent(in) :: name
          character(len=*), intent(in) :: units
          character(len=*), intent(in) :: long_name
          real(kind=rk), intent(in) :: remin
          real(kind=rk), intent(in) :: q10
          integer, intent(in) :: source_depth_distribution
          type(type_bulk_standard_variable), intent(in) :: aggregate_target
          real(kind=rk), optional,intent(in) :: aggregate_scale_factor
.. code-block:: bash 

      subroutine layer_do_bottom(self,e__arguments_do_bottom_e) 
.. code-block:: bash 

      subroutine layer_process_constituent(self,e__arguments_do_bottom_e,info) 
          type(type_constituent_in_single_layer), intent(in) :: info
.. code-block:: bash 

      subroutine constituent_for_single_layer_change_do_bottom(self,e__arguments_do_bottom_e) 
.. code-block:: bash 

      function partq(d_pen,d_top,d_bot,d_max) 
          real(kind=rk), intent(in) :: d_pen
          real(kind=rk), intent(in) :: d_top
          real(kind=rk), intent(in) :: d_bot
          real(kind=rk), intent(in) :: d_max

benthic_erosion module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine do_bottom(self,e__arguments_do_bottom_e) 

fluff module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine do_bottom(self,e__arguments_do_bottom_e) 

benthic_nitrogen_cycle module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine do_bottom(self,e__arguments_do_bottom_e) 

benthic_column_dissolved_matter module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine benthic_dissolved_matter_initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      subroutine initialize_constituent(self,info,profile,profile_info,name,units,long_name,aggregate_target,background_value,nonnegative) 
          type(type_single_constituent), intent(inout),target :: info
          type(type_single_constituent_estimates), intent(inout),target :: profile_info
          character(len=*), intent(in) :: name
          character(len=*), intent(in) :: units
          character(len=*), intent(in) :: long_name
          type(type_bulk_standard_variable), optional,intent(in) :: aggregate_target
          real(kind=rk), optional,intent(in) :: background_value
          logical, optional,intent(in) :: nonnegative
.. code-block:: bash 

      subroutine benthic_dissolved_matter_do_bottom(self,e__arguments_do_bottom_e) 
.. code-block:: bash 

      subroutine process_constituent(self,e__arguments_do_bottom_e,info) 
          type(type_single_constituent), intent(in) :: info
.. code-block:: bash 

      subroutine compute_equilibrium_profile(sigma,c0,p,p_deep,d,c_bot,c_int) 
          real(kind=rk), intent(in) :: sigma
          real(kind=rk), intent(in) :: c0
          real(kind=rk), intent(in) :: p
          real(kind=rk), intent(in) :: p_deep
          real(kind=rk), intent(in) :: d
          real(kind=rk), intent(out) :: c_bot
          real(kind=rk), intent(out) :: c_int
.. code-block:: bash 

      subroutine compute_final_equilibrium_profile(sigma,c0,p,p_deep,dmax,d,c_int) 
          real(kind=rk), intent(in) :: sigma
          real(kind=rk), intent(in) :: c0
          real(kind=rk), intent(in) :: p
          real(kind=rk), intent(in) :: p_deep
          real(kind=rk), intent(in) :: dmax
          real(kind=rk), intent(out) :: d
          real(kind=rk), intent(out) :: c_int
.. code-block:: bash 

      subroutine dissolved_matter_per_layer_do_bottom(self,e__arguments_do_bottom_e) 

carbonate module with the following subroutines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash 

      subroutine initialize(self,configunit) 
          integer, intent(in) :: configunit
.. code-block:: bash 

      function approximate_alkalinity(iswtalk,t,s) result (ta) 
          integer, intent(in) :: iswtalk
          real(kind=rk), intent(in) :: t
          real(kind=rk), intent(in) :: s
          real(kind=rk) :: ta
.. code-block:: bash 

      subroutine do_bn(self,e__arguments_do_e) 
.. code-block:: bash 

      subroutine do_surface(self,e__arguments_do_surface_e) 
.. code-block:: bash 

      subroutine co2dyn(t,s,prss,ctot,ta,ph,pco2,h2co3,hco3,co3,k0co2,success,hscale) 
          real(kind=rk), intent(in) :: t
          real(kind=rk), intent(in) :: s
          real(kind=rk), intent(in) :: prss
          real(kind=rk), intent(inout) :: ctot
          real(kind=rk), intent(inout) :: ta
          real(kind=rk), intent(out) :: ph
          real(kind=rk), intent(out) :: pco2
          real(kind=rk), intent(out) :: h2co3
          real(kind=rk), intent(out) :: hco3
          real(kind=rk), intent(out) :: co3
          real(kind=rk), intent(out) :: k0co2
          logical, intent(out) :: success
          integer, intent(in) :: hscale
.. code-block:: bash 

      subroutine co2set(p,t,s,k0co2,k1co2,k2co2,kb,hscale) 
          real(kind=rk), intent(in) :: p
          real(kind=rk), intent(in) :: t
          real(kind=rk), intent(in) :: s
          real(kind=rk), intent(out) :: k0co2
          real(kind=rk), intent(out) :: k1co2
          real(kind=rk), intent(out) :: k2co2
          real(kind=rk), intent(out) :: kb
          integer, intent (in) :: hscale
.. code-block:: bash 

      subroutine co2clc(k0co2,k1co2,k2co2,kb,icalc,boron,btot,ctot,ta,ph,pco2,h2co3,hco3,co3,success) 
          real(kind=rk), intent(in) :: k0co2
          real(kind=rk), intent(in) :: k1co2
          real(kind=rk), intent(in) :: k2co2
          real(kind=rk), intent(in) :: kb
          integer :: icalc
          logical :: boron
          real(kind=rk) :: btot
          real(kind=rk), intent(inout) :: ctot
          real(kind=rk), intent(inout) :: ta
          real(kind=rk), intent(inout) :: ph
          real(kind=rk), intent(inout) :: pco2
          real(kind=rk), intent(inout) :: h2co3
          real(kind=rk), intent(inout) :: hco3
          real(kind=rk), intent(inout) :: co3
          logical, intent(out) :: success
.. code-block:: bash 

      subroutine caco3_saturation(tc,s,pr,co3,om_cal,om_arg) 
          real(kind=rk), intent(in) :: tc
          real(kind=rk), intent(in) :: s
          real(kind=rk), intent(in) :: pr
          real(kind=rk), intent(in) :: co3
          real(kind=rk), intent(out) :: om_cal
          real(kind=rk), intent(out) :: om_arg