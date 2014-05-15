module pml_ersem_shared

   use fabm_types

   implicit none

   private
   
   real(rk),parameter,public :: CMass   = 12.011_rk
   real(rk),parameter,public :: qnRPIcX = 1.26E-02_rk
   real(rk),parameter,public :: qpRPIcX = 7.86E-04_rk
   real(rk),parameter,public :: qsRPIcX = 15._rk/106._rk/CMass
   real(rk),parameter,public :: ZeroX   = 1e-8_rk
end module