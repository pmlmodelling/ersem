.. _gotm:

###############
GOTM-FABM-ERSEM
###############


Building GOTM-FABM-ERSEM
~~~~~~~~~~~~~~~~~~~~~~~~

To run the install script, you will need to have ``netCDF`` installed.
An example of how to do this is here:

.. literalinclude:: ../../../github-actions/gotm-fabm-ersem/gotm-fabm-ersem-dep-debian.sh
    :language: bash
    :linenos:

To install GOTM-FABM-ERSEM, we suggest you use the following script below

.. literalinclude:: ../../../github-actions/gotm-fabm-ersem/gotm-fabm-ersem-build.sh
    :language: bash
    :linenos:


If you experience NetCDF issues when running ``make install``, see `tips
and tricks/troubleshooting <#tips-and-trickstroubleshooting>`__.

Now you should have a GOTM executable with FABM and ERSEM support at
``~/local/gotm/bin/gotm``.


Running GOTM-FABM-ERSEM
~~~~~~~~~~~~~~~~~~~~~~~

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
