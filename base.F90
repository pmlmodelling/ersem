#include "fabm_driver.h"
#define IRON
module pml_ersem_base

   use fabm_types

   implicit none

!  default: all is private.
   private

   type,extends(type_base_model),public :: type_ersem_base_model
      type (type_state_variable_id)      :: id_c,id_n,id_p,id_f,id_s,id_chl
      type (type_diagnostic_variable_id) :: id_cD,id_nD,id_pD,id_fD,id_sD,id_chlD
      type (type_conserved_quantity_id)  :: id_totc,id_totn,id_totp,id_tots,id_totf
   contains
      procedure :: initialize_ersem_base
   end type

   type type_component
      character(len=attribute_length) :: name   = ''
      real(rk)                        :: weight = 1._rk
      type (type_dependency_id)       :: id
      type (type_component),pointer   :: next   => null()
   end type

   type,extends(type_base_model) :: type_weighted_sum
      character(len=attribute_length) :: output_name      = ''
      character(len=attribute_length) :: output_long_name = ''
      character(len=attribute_length) :: output_units     = ''
      type (type_diagnostic_variable_id) :: id_output
      type (type_component),pointer   :: first => null()
   contains
      procedure :: initialize    => weighted_sum_initialize
      procedure :: do            => weighted_sum_do
      procedure :: add_component => weighted_sum_add_component
   end type
      
contains
      
   subroutine initialize_ersem_base(self,c_ini,n_ini,p_ini,s_ini,f_ini,chl_ini,qn,qp)
      class (type_ersem_base_model), intent(inout), target :: self
      real(rk),optional,             intent(in)            :: c_ini,n_ini,p_ini,s_ini,f_ini,chl_ini
      real(rk),optional,             intent(in)            :: qn,qp
      
      type (type_weighted_sum),pointer :: child

      call self%register_conserved_quantity(self%id_totc,standard_variables%total_carbon)
      call self%register_conserved_quantity(self%id_totn,standard_variables%total_nitrogen)
      call self%register_conserved_quantity(self%id_totp,standard_variables%total_phosphorus)
      call self%register_conserved_quantity(self%id_tots,standard_variables%total_silicate)
      call self%register_conserved_quantity(self%id_totf,standard_variables%total_iron)

      if (present(c_ini)) then
         call self%register_state_variable(self%id_c,'c','mg m-3','carbon',c_ini,minimum=0._rk)
         call self%add_conserved_quantity_component(self%id_totc,self%id_c,scale_factor=1._rk/12._rk)
      else
         call self%register_diagnostic_variable(self%id_cD,'c','mg m-3','carbon',missing_value=0._rk)
      end if

      if (present(s_ini)) then
         call self%register_state_variable(self%id_s,'s','mmol Si m-3','silicate',s_ini,minimum=0._rk)
         call self%add_conserved_quantity_component(self%id_tots,self%id_s)
      else
         call self%register_diagnostic_variable(self%id_sD,'s','mmol Si m-3','silicate',missing_value=0._rk)
      end if

      if (present(chl_ini)) then
         call self%register_state_variable(self%id_chl,'Chl','mg C m-3','chlorophyll-a',chl_ini,minimum=0._rk)
      else
         call self%register_diagnostic_variable(self%id_chlD,'Chl','mg C m-3','chlorophyll-a',missing_value=0._rk)
      end if

#ifdef IRON   
      if (present(f_ini)) then
         ! Iron as state variable
         call self%register_state_variable(self%id_f,'f','umol Fe m-3','iron',f_ini,minimum=0._rk)
         call self%add_conserved_quantity_component(self%id_totf,self%id_f)
      else
         ! No iron
         call self%register_diagnostic_variable(self%id_fD,'f','umol Fe m-3','iron',missing_value=0._rk)
      end if
#endif

      if (present(p_ini)) then
         ! Phosphorous as state variable
         call self%register_state_variable(self%id_p, 'p', 'mmol P m-3', 'phosphorus', p_ini, minimum=0._rk)
         call self%add_conserved_quantity_component(self%id_totp,self%id_p)
      elseif (present(qp)) then
         ! Phosphorous derived from carbon (constant C:P)
         ! Create an expression that computes P content, so external models (e.g., predators) can access it.
         allocate(child)
         child%output_name = 'p'
         child%output_units = 'mmol P m-3'
         call child%add_component('c',qp)
         call self%add_child(child,'dynamic_p',configunit=-1)
         call self%add_conserved_quantity_component(self%id_totp,self%id_c,scale_factor=qp)
      else
         ! No phosphorous
         call self%register_diagnostic_variable(self%id_pD, 'p', 'mmol P m-3', 'phosphorus', missing_value=0._rk)
      end if

      if (present(n_ini)) then
         ! Nitrogen as state variable
         call self%register_state_variable(self%id_n, 'n', 'mmol N m-3', 'nitrogen', n_ini, minimum=0._rk)
         call self%add_conserved_quantity_component(self%id_totn,self%id_n)
      elseif (present(qn)) then
         ! Nitrogen derived from carbon (constant C:P)
         ! Create an expression that computes N content, so external models (e.g., predators) can access it.
         allocate(child)
         child%output_name = 'n'
         child%output_units = 'mmol N m-3'
         call child%add_component('c',qn)
         call self%add_child(child,'dynamic_n',configunit=-1)
         call self%add_conserved_quantity_component(self%id_totn,self%id_c,scale_factor=qn)
      else
         ! No nitrogen
         call self%register_diagnostic_variable(self%id_nD, 'n', 'mmol N m-3', 'nitrogen', missing_value=0._rk)
      end if
   end subroutine

   subroutine weighted_sum_initialize(self,configunit)
      class (type_weighted_sum),intent(inout),target :: self
      integer,                  intent(in)           :: configunit

      type (type_component),pointer :: component
      component => self%first
      do while (associated(component))
         call self%register_dependency(component%id,trim(component%name))
         component => component%next
      end do
      if (self%output_long_name=='') self%output_long_name = self%output_name
      call self%register_diagnostic_variable(self%id_output,self%output_name,self%output_units,self%output_long_name)
      call self%parent%add_alias(self%id_output,trim(self%output_name))
   end subroutine

   subroutine weighted_sum_add_component(self,name,weight)
      class (type_weighted_sum),intent(inout) :: self
      character(len=*),         intent(in)    :: name
      real(rk),optional,        intent(in)    :: weight
      
      type (type_component),pointer :: component

      if (.not.associated(self%first)) then
         allocate(self%first)
         component => self%first
      else
         component => self%first
         do while (associated(component%next))
            component => component%next
         end do
         allocate(component%next)
         component => component%next
      end if
      component%name = name
      if (present(weight)) component%weight = weight
   end subroutine

   subroutine weighted_sum_do(self,_ARGUMENTS_DO_)
      class (type_weighted_sum),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      
      type (type_component),pointer :: component
      real(rk)                      :: sum,value
      
      _LOOP_BEGIN_
         sum = 0._rk
         component => self%first
         do while (associated(component))
            _GET_(component%id,value)
            sum = sum + component%weight*value
            component => component%next
         end do
         _SET_DIAGNOSTIC_(self%id_output,sum)
      _LOOP_END_
   end subroutine
end module