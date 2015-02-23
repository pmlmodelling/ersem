module ersem_shared

   use fabm_types
   use fabm_standard_variables

   implicit none

   public

   ! Parameter to reproduce pre-SSB-v1 ERSEM behaviour, including several obvious flaws/inconsistencies.
   logical, parameter :: legacy_ersem_compatibility = .true.

   real(rk),parameter :: CMass   = 12.011_rk
   real(rk),parameter :: qnRPIcX = 1.26E-02_rk
   real(rk),parameter :: qpRPIcX = 7.86E-04_rk
   real(rk),parameter :: qsRPIcX = 15._rk/106._rk/CMass
   real(rk),parameter :: ZeroX   = 1e-8_rk
   real(rk), parameter :: pi=acos(-1._rk)
   real(rk), parameter :: deg2rad=pi/180._rk

   logical :: use_light_iop = .false.
   logical :: use_light = .false.
#ifdef IRON
   logical,parameter :: use_iron = .true.
#else
   logical,parameter :: use_iron = .false.
#endif

   ! Aggregate diagnostics for e.g., carbon budgets.
   type (type_bulk_standard_variable),parameter :: photosynthesis_rate = type_bulk_standard_variable(name='photosynthesis_rate',units='mg C/m^3/d',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: phytoplankton_respiration_rate = type_bulk_standard_variable(name='phytoplankton_respiration_rate',units='mg C/m^3/d',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: zooplankton_respiration_rate = type_bulk_standard_variable(name='zooplankton_respiration_rate',units='mg C/m^3/d',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: bacterial_respiration_rate = type_bulk_standard_variable(name='bacterial_respiration_rate',units='mg C/m^3/d',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: total_chlorophyll = type_bulk_standard_variable(name='total_chlorophyll',units='mg/m^3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: total_calcite_in_biota = type_bulk_standard_variable(name='total_calcite_in_biota',units='mg C/m^3',aggregate_variable=.true.)

   ! Aggregate variables for benthic bioturbation and bioirrigation (summed over all fauna).
   type (type_bulk_standard_variable),parameter :: total_bioturbation_activity = type_bulk_standard_variable(name='total_bioturbation_activity',units='mg C/m^2/d',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: total_bioirrigation_activity = type_bulk_standard_variable(name='total_bioirrigation_activity',units='mg C/m^2/d',aggregate_variable=.true.)

   ! Standard benthic variables used to make implicit based on matching standard names coupling possible.
   type (type_horizontal_standard_variable),parameter :: depth_of_sediment_column = type_horizontal_standard_variable(name='depth_of_sediment_column',units='m')
   type (type_horizontal_standard_variable),parameter :: diffusivity_in_sediment_layer_1 = type_horizontal_standard_variable(name='diffusivity_in_sediment_layer_1',units='m^2/d')
   type (type_horizontal_standard_variable),parameter :: diffusivity_in_sediment_layer_2 = type_horizontal_standard_variable(name='diffusivity_in_sediment_layer_2',units='m^2/d')
   type (type_horizontal_standard_variable),parameter :: diffusivity_in_sediment_layer_3 = type_horizontal_standard_variable(name='diffusivity_in_sediment_layer_3',units='m^2/d')
   type (type_horizontal_standard_variable),parameter :: particulate_diffusivity_due_to_bioturbation = type_horizontal_standard_variable(name='particulate_diffusivity_due_to_bioturbation',units='m^2/d')
   type (type_horizontal_standard_variable),parameter :: bioturbation_depth = type_horizontal_standard_variable(name='bioturbation_depth',units='m')
   type (type_horizontal_standard_variable),parameter :: sediment_porosity = type_horizontal_standard_variable(name='sediment_porosity',units='-')
   type (type_horizontal_standard_variable),parameter :: depth_of_bottom_interface_of_layer_1 = type_horizontal_standard_variable(name='depth_of_bottom_interface_of_layer_1',units='m')
   type (type_horizontal_standard_variable),parameter :: depth_of_bottom_interface_of_layer_2 = type_horizontal_standard_variable(name='depth_of_bottom_interface_of_layer_2',units='m')
   type (type_horizontal_standard_variable),parameter :: pelagic_benthic_transfer_constant = type_horizontal_standard_variable(name='pelagic_benthic_transfer_constant',units='d/m')

   ! Aggregate absorption and backscatter.
   type (type_bulk_standard_variable),parameter :: particulate_organic_absorption_coefficient = type_bulk_standard_variable(name='particulate_organic_absorption_coefficient',units='1/m',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: particulate_organic_backscatter_coefficient = type_bulk_standard_variable(name='particulate_organic_backscatter_coefficient',units='1/m',aggregate_variable=.true.)

   ! Zenith angle.
   type (type_horizontal_standard_variable),parameter :: zenith_angle = type_horizontal_standard_variable(name='zenith_angle',units='degrees')

end module
