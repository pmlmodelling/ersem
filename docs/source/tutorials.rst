.. _tutorials:


ERSEM tutorials
===============


GOTM: ERSEM in a water column
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To use GOTM with ERSEM, copy an ERSEM configuration (recommended:
``testcases/fabm-ersem-15.06-L4-ben-docdyn-iop.yaml``) to your GOTM
setup directory. This file needs to be named ``fabm.yaml``. You also
need to prescribe two additional forcing variables: atmospheric pCO2 (in
ppm) and light absorption by sediment (1/m). This can be done in
``gotm.yaml`` under the ``fabm`` section, by adding something like:

::

   input:
     mole_fraction_of_carbon_dioxide_in_air: 411.29
     absorption_of_silt: 0.07

It is good practice to keep up to date with the latest code from the
ERSEM, FABM and GOTM repositories by regularly running ``git pull`` in
all directories under ``<SOURCEDIR>``.

If the ERSEM, FABM or GOTM source codes change (e.g., because you made
changes to the code yourself or ran ``git pull``), you need to
recompile. This does *not* require rerunning cmake. Instead, return to
the build directory and rerun ``make install``. For instance
``cd gotm && make install``.

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

Python front-end
~~~~~~~~~~~~~~~~

Example Jupyter notebooks that use the Python front-end can be found in
``<FABMDIR>/testcases/python``. More examples can be found on `the FABM
wiki <https://github.com/fabm-model/fabm/wiki/python>`__.

NEMO + FABM + ERSEM
~~~~~~~~~~~~~~~~~~~

This requires a customised NEMO codebase with the FABM coupler
integrated. To obtain access to this code, please
`register <https://pml.ac.uk/Modelling_at_PML/Access_Code>`__.

FABM-ERSEM needs to be compiled separately before the compilation of
NEMO. Usually, the following suffices for this:

::

   mkdir -p ~/build/nemo
   cd ~/build/nemo
   cmake <SOURCEDIR>/fabm -DFABM_HOST=nemo -DFABM_ERSEM_BASE=<SOURCEDIR>/ersem
   make install

This will create the library in the standard folder
``~/local/fabm/nemo/lib`` where NEMO-FABM will look for linking to NEMO.
Module files (.mod) are placed at ``~/local/fabm/nemo/include``.


