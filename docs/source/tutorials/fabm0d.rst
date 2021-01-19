.. _fabm0d:

fabm0d: ERSEM in an aquarium
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To use the driver, you need a configuration file can ``run.nml``, which you
could take from ``<SOURCEDIR>/fabm/testcases/0d/run.nml``. In addition,
you need a file with forcing data for surface shortwave radiation,
temperature and salinity, e.g.,
``<SOURCEDIR>/fabm/testcases/0d/env_nns_annual.dat``. The local name of
this file is configured in ``run.nml``, variable ``env_file``. Finally,
you can provide additional forcing with a ``input.yaml`` file, which can
contain entries such as

::

   wind_speed:
     file: wind.dat
     column: 1
     scale_factor: 1.0
   absorption_of_silt:
     constant_value: 0.07
   mole_fraction_of_carbon_dioxide_in_air:
     constant_value: 411.29
   bottom_stress:
     constant_value: 0.0

This specifies that wind speed must be read from file wind.dat (first
column, no scaling), the light absorption by silt must be set to a
constant value of 0.07 m-1, atmospheric pCO2 to 411.29 ppm, and bottom
stress to 0 Pa. The “column” and “scale_factor” attributes in the case
of wind_speed are optional. They default to 1 and 1.0, respectively.
