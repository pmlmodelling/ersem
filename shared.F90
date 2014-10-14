module pml_ersem_shared

   use fabm_types
   use fabm_standard_variables

   implicit none

   public
   
   real(rk),parameter :: CMass   = 12.011_rk
   real(rk),parameter :: qnRPIcX = 1.26E-02_rk
   real(rk),parameter :: qpRPIcX = 7.86E-04_rk
   real(rk),parameter :: qsRPIcX = 15._rk/106._rk/CMass
   real(rk),parameter :: ZeroX   = 1e-8_rk

   type (type_bulk_standard_variable),parameter :: photosynthesis_rate = type_bulk_standard_variable(name='photosynthesis_rate',units='mg C/m^3/d',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: phytoplankton_respiration_rate = type_bulk_standard_variable(name='phytoplankton_respiration_rate',units='mg C/m^3/d',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: zooplankton_respiration_rate = type_bulk_standard_variable(name='zooplankton_respiration_rate',units='mg C/m^3/d',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: bacterial_respiration_rate = type_bulk_standard_variable(name='bacterial_respiration_rate',units='mg C/m^3/d',aggregate_variable=.true.)

end module