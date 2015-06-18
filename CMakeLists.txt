add_subdirectory(src)
set_property(TARGET fabm_models_ersem PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
