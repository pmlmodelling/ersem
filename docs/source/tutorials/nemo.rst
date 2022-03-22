.. _nemo:


#############################
NEMO: Atlantic Margin Model 7
#############################

This tutorial gives an end to end example how to install and then use
NEMO-ERSEM on the Atlantic Margin Model 7 (AMM7) using a high performance computing
(HPC) machine. We have used
`ARCHER2 <https://www.archer2.ac.uk/>`__ together with a 
`singularity container <https://sylabs.io/guides/3.5/user-guide/introduction.html>`__ with
NEMO-ERSEM install on it.

This tutorial is based on utilising the 
`container <https://sylabs.io/guides/3.5/user-guide/introduction.html>`_ with the configuration
scripts container `here <https://github.com/dalepartridge/AMM7-NEMO4-FABM-setup>`_. You will
also need access to the shared folder within the ARCHER2 project id `n01` to obtain the forcing
files.

.. note::
    You will need to get access to the 
    `NEMO-ERSEM code <https://www.pml.ac.uk/Modelling_at_PML/Access_Code>`_.

Obtain NEMO-ERSEM container
###########################

The key packages for that are installed into 
`NEMO-ERSEM container <https://github.com/pmlmodelling/NEMO-container>`_ are:

* `NEMO <https://gitlab.ecosystem-modelling.pml.ac.uk/gle/nemo404>`__ (UK version -- see note)
* `XIOS <http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-2.5>`__
* `FABM <https://github.com/fabm-model/fabm>`__
* `ERSEM <https://github.com/pmlmodelling/ersem>`__

ERSEM and FABM are freely available on GitHub and XIOS is available through the NEMO consotium
svn server, however, for NEMO you will have to register for the code -- see note above. Full instructions 
how to build the container can be found `here <https://github.com/pmlmodelling/NEMO-container>`__.

Running the container on ARCHER2
################################

To generate the `slurm <https://slurm.schedmd.com/documentation.html>`__ scheduling scripts to 
run the model follow the instructions 
`here <https://docs.archer2.ac.uk/research-software/nemo/nemo/#building-a-run-script>`__. 

The slurm script schedulers, 
`mkslurm <https://docs.archer2.ac.uk/research-software/nemo/nemo/#building-a-run-script>`__,
requires two executables, namely `xios_server.exe` and `nemo`. The script below shows an example how
you set up the singulartiy container, `nemo.sif`, to be used as an executable. 

.. code-block:: bash
    :caption: Example script to run container, user is required to change `[RUNING-DIR]` and 
              `[INPUT-DIR]` to the coresponding directories and `[NEMO-XIOS]` to either `nemo`
              or `xios`.

    #! /bin/bash
    BIND_OPTS="-B [RUNING-DIR],"
    BIND_OPTS="${BIND_OPTS}[INPUT-DIR],"
    BIND_OPTS="${BIND_OPTS}/usr/lib64/liblustreapi.so:/opt/lib64/liblustreapi.so,"
    BIND_OPTS="${BIND_OPTS}/usr/lib64/liblustreapi.so.1:/opt/lib64/liblustreapi.so.1,"
    BIND_OPTS="${BIND_OPTS}/usr/lib64/liblustreapi.so.1.0.0:/opt/lib64/liblustreapi.so.1.0.0,"
    BIND_OPTS="${BIND_OPTS}/opt/cray/libfabric/1.11.0.4.71/lib64/:/opt/fabric,"
    BIND_OPTS="${BIND_OPTS}/usr/lib64/libpals.a:/opt/lib64/libpals.a,"
    BIND_OPTS="${BIND_OPTS}/usr/lib64/libpals.so:/opt/lib64/libpals.so,"
    BIND_OPTS="${BIND_OPTS}/usr/lib64/libpals.so.0:/opt/lib64/libpals.so.0,"
    BIND_OPTS="${BIND_OPTS}/usr/lib64/libpals.so.0.0.0:/opt/lib64/libpals.so.0.0.0,"
    BIND_OPTS="${BIND_OPTS}/opt/cray/pe/mpich/8.1.9/ofi/gnu/9.1:/opt/mpi/install/,"
    BIND_OPTS="${BIND_OPTS}/opt/cray/pe/pmi/6.0.13/lib,"
    BIND_OPTS="${BIND_OPTS}/var/spool/slurmd/mpi_cray_shasta"
    
    export SINGULARITYENV_LD_LIBRARY_PATH=/opt/hdf5/install/lib:$SINGULARITYENV_LD_LIBRARY_PATH
    export SINGULARITYENV_LD_LIBRARY_PATH=/opt/mpi/install/lib-abi-mpich:$SINGULARITYENV_LD_LIBRARY_PATH
    export SINGULARITYENV_LD_LIBRARY_PATH=/.singularity.d/libs:$SINGULARITYENV_LD_LIBRARY_PATH
    export SINGULARITYENV_LD_LIBRARY_PATH=/opt/cray/pe/pmi/6.0.13/lib:$SINGULARITYENV_LD_LIBRARY_PATH
    export SINGULARITYENV_LD_LIBRARY_PATH=/opt/lib64:$SINGULARITYENV_LD_LIBRARY_PATH
    export SINGULARITYENV_LD_LIBRARY_PATH=/opt/fabric:$SINGULARITYENV_LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=/opt/hdf5/install/lib
    export OMP_NUM_THREADS=1
    
    singularity run ${BIND_OPTS} --home [RUNING-DIR] nemo.sif [NEMO-XIOS]

The CoNES `documentation <https://cones.readthedocs.io/en/latest/?badge=latest>`__ gives additional 
information how to run NEMO singularity containers.


Example output from NEMO-ERSEM
##############################

To visualise NEMO output, we recomend using `nctoolkit <https://github.com/pmlmodelling/nctoolkit>`__.
A basic example of how to use `nctoolkit` is given below. Within the 
`documentation <https://nctoolkit.readthedocs.io/en/latest/>`__ of `nctoolkit` one will find additional
example uses.

.. code-block:: python
    :caption: Example plotting script using `nctoolkit`, user is required to change `[PATH_TO_NETCDF_FILE]` and 
            `[VARIABLE]` to the exact location of the NEMO output and the variable to be plotted, respectively.

    import nctoolkit as nc

    ff = [PATH_TO_NETCDF_FILE]
    ds = nc.open_data(ff)
    plot = ds.plot([VARIABLE])

The following plots represent a months output of a NEMO-ERSEM simulation on the AMM7 domain.

.. dropdown:: Potential temperature

	.. image:: ../../images/NEMO_temp.gif

.. dropdown::  Salinity

	.. image:: ../../images/NEMO_sal.gif

.. dropdown:: Phosphate phosphorus

	.. image:: ../../images/NEMO_N1_p.gif

.. dropdown::  Nitrate nitrogen

	.. image:: ../../images/NEMO_N3_n.gif

.. dropdown:: Carbonate total dissolved inorganic carbon

	.. image:: ../../images/NEMO_O3_c.gif

.. dropdown:: Diatoms chlorophyll

	.. image:: ../../images/NEMO_P1_Chl.gif

.. dropdown:: Medium-sized POM carbon

	.. image:: ../../images/NEMO_R6_c.gif

.. dropdown:: Oxygen

	.. image:: ../../images/NEMO_O2_o.gif

