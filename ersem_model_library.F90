module ersem_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   use ersem_pelagic_base
   use ersem_benthic_base
   use ersem_oxygen
   use ersem_carbonate
   use ersem_primary_producer
   use ersem_microzooplankton
   use ersem_mesozooplankton
   use ersem_bacteria
   use ersem_bacteria_docdyn
   use ersem_nitrification
   use ersem_light
   use ersem_light_iop
   use ersem_calcification
   use ersem_benthic_column
   use ersem_benthic_column_dissolved_matter
   use ersem_benthic_column_particulate_matter
   use ersem_benthic_nitrogen_cycle
   use ersem_benthic_bacteria
   use ersem_benthic_fauna
   use ersem_zenith_angle
   use ersem_shared

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: ersem_factory

contains

   subroutine create(self,name,model)
      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('pelagic_base');                            allocate(type_ersem_pelagic_base::model)
         case ('benthic_base');                            allocate(type_ersem_benthic_base::model)
         case ('oxygen');                                  allocate(type_ersem_oxygen::model)
         case ('carbonate');                               allocate(type_ersem_carbonate::model)
         case ('primary_producer');                        allocate(type_ersem_primary_producer::model)
         case ('microzooplankton');                        allocate(type_ersem_microzooplankton::model)
         case ('mesozooplankton');                         allocate(type_ersem_mesozooplankton::model)
         case ('bacteria');                                allocate(type_ersem_bacteria::model)
         case ('bacteria_docdyn');                         allocate(type_ersem_bacteria_docdyn::model)
         case ('nitrification');                           allocate(type_ersem_nitrification::model)
         case ('light');                                   allocate(type_ersem_light::model); use_light=.true.
         case ('light_iop');                                   allocate(type_ersem_light_iop::model); use_light_iop=.true.
         case ('calcification');                           allocate(type_ersem_calcification::model)
         case ('benthic_column');                          allocate(type_ersem_benthic_column::model)
         case ('benthic_column_dissolved_matter');         allocate(type_ersem_benthic_column_dissolved_matter::model)
         case ('benthic_column_particulate_matter');       allocate(type_ersem_benthic_column_particulate_matter::model)
         case ('benthic_nitrogen_cycle');                  allocate(type_ersem_benthic_nitrogen_cycle::model)
         case ('benthic_bacteria');                        allocate(type_ersem_benthic_bacteria::model)
         case ('benthic_column_particulate_matter_layer'); allocate(type_ersem_benthic_pom_layer::model)
         case ('benthic_fauna');                           allocate(type_ersem_benthic_fauna::model)
         case ('zenith_angle');                           allocate(type_ersem_zenith_angle::model)
         ! Add new models here
         case default
            call self%type_base_model_factory%create(name,model)
      end select
   end subroutine create

end module
