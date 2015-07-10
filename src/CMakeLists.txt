option(ERSEM_USE_IRON "Include support for iron in ERSEM" OFF)
if(ERSEM_USE_IRON)
  add_definitions(-DIRON)
endif()

add_library(fabm_models_ersem OBJECT
            ersem_model_library.F90
            shared.F90
            pelagic_base.F90
            benthic_base.F90
            light.F90
            light_iop.F90
            primary_producer.F90
            microzooplankton.F90
            mesozooplankton.F90
            bacteria.F90
            bacteria_docdyn.F90
            calcification.F90
            carbonate.F90
            oxygen.F90
            nitrification.F90
            benthic_column.F90
            benthic_column_dissolved_matter.F90
            benthic_column_particulate_matter.F90
            benthic_nitrogen_cycle.F90
            benthic_bacteria.F90
            benthic_fauna.F90
            benthic_carbonate.F90
            zenith_angle.F90
           )

add_dependencies(fabm_models_ersem fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2015 - PML
#-----------------------------------------------------------------------
