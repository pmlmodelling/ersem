option(ERSEM_USE_IRON "Include support for iron in ERSEM" OFF)
if(ERSEM_USE_IRON)
  add_definitions(-DIRON)
endif()

# Generate placeholder "ersem_version" module (versions all set to "unknown" - updated at build time if FABM_EMBED_VERSION is set)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.F90.in ${CMAKE_CURRENT_BINARY_DIR}/version.F90)

add_library(fabm_models_ersem OBJECT
            ersem_model_library.F90
            ${CMAKE_CURRENT_BINARY_DIR}/version.F90
            shared.F90
            pelagic_base.F90
            benthic_base.F90
            light.F90
            light_iop.F90
            light_iop_ady.F90
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
            benthic_calcite.F90
            benthic_erosion.F90
            fluff.F90
            zenith_angle.F90
           )

add_dependencies(fabm_models_ersem fabm_base)

if(FABM_EMBED_VERSION)
  # Dynamically generate "ersem_version" module with git commit/branch information at every build.
  add_custom_target(fabm_models_ersem_version
     ${CMAKE_COMMAND} -DINFILE=${CMAKE_CURRENT_SOURCE_DIR}/version.F90.in -DOUTFILE=${CMAKE_CURRENT_BINARY_DIR}/version.F90 -P ${GET_GIT_INFO_SCRIPT}
     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
     COMMENT "Retrieving description of last ERSEM commit..."
     SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/version.F90.in ${GET_GIT_INFO_SCRIPT}
     VERBATIM
  )
  add_dependencies(fabm_models_ersem fabm_models_ersem_version)
endif()

#-----------------------------------------------------------------------
# Copyright (C) 2015 - PML
#-----------------------------------------------------------------------
