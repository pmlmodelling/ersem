
.. _nemo:

#############################
NEMO-ERSEM in a water column 
#############################

This tutorial demonstrates how to configure and run NEMO 1D (water-column)
testcase C1D_PAPA with ERSEM biogeochemistry. It is a demonstration of
concept only, as it combines hydrodynamics of PAPA station in the North
Pacific Ocean with biogeochemistry of L4 station in the Western English
Channel. Users are encouraged to modify this example for their own purposes.

Step 1: Obtaining the code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To get started, you need to follow the instructions for setting up
the `C1D_PAPA configuration <https://forge.ipsl.jussieu.fr/nemo/chrome/site/doc/NEMO/guide/html/cfgs.html#c1d-papa>`__. The input data must be obtained from NEMO Reference configurations inputs repository on `Zenodo <https://zenodo.org/record/1472245#.Yt6_QIzMKEI>`__.

NEMO4 code base with FABM support can be obtained in the corresponding `repository <https://github.com/pmlmodelling/NEMO4.0-FABM>`__. You will also need to download `FABM <https://github.com/fabm-model/fabm>`__ and `ERSEM <https://github.com/pmlmodelling/ersem>`__. Finally, you need to download I/O server `XIOS-2.5 <https://forge.ipsl.jussieu.fr/nemo/chrome/site/doc/NEMO/guide/html/install.html#extract-and-install-xios>`__, and `install it. <https://forge.ipsl.jussieu.fr/ioserver/>`__. Copy xios_server.exe executable into your working directory.


Step 2: Compiling the code
~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, you need to compile FABM with ERSEM support and specifying nemo as a physical host. The following commands can be written in an executable file for easy recompilation at any future point:

  .. code-block:: bash
        
        old = `pwd`                      #remember current directory
        mkdir -p ~/build/nemo-fabm-ersem # create directory for the build
        cd ~/build/nemo-fabm-ersem       # go to the build directory
        cmake <FABM sourcedir> -DFABM_HOST=nemo -DFABM_ERSEM_BASE=<ERSEM sourcedir> -DCMAKE_INSTALL_PREFIX=~/local/fabm/nemo-fabm-ersem
        #replace <FABM sourcedir> and <ERSEM sourcedir> with corresponding directories of FABM and ERSEM code bases.
        make install
        make -j4
        cd $old                         # return to the original directory
        
After that, you need to compile NEMO. C1D_PAPA_FABM configuration has been added to the NEMO4.0-FABM to support this, or you can create your own. The critical point is to provide the necessary compilation keys in cpp_X.fcm file, i.e. key_c1d for compilation in 1D, and key_fabm for FABM support:

  .. code-block:: bash
  
       bld::tool::fppkeys   key_c1d key_mpp_mpi key_iomput key_nosignedzero key_top key_fabm
       
Next, compile the model by executing the following lines:

  .. code-block:: bash
  
    #!/bin/bash

    module load mpi

    NEMO_BUILD_DIR=$<NEMODIR>
    RUNDIR=$<MYRUNDIR>
    export XIOS_HOME=$<XIOSDIR>
    export FABM_HOME=$HOME/local/fabm/nemo-fabm-ersem
    # replace <NEMODIR> and <XIOSDIR> with location of corresponding code bases on your machine, <MYRUNDIR> with your working directory. FABM_HOME in this example corresponds to directory where we installed FABM-ERSEM above.
    
    ARCH=GCC_PMPC

    cd $NEMO_BUILD_DIR
    ./makenemo -m $ARCH -r C1D_PAPA_FABM -n C1D_PAPA_FABM_BLD_SCRATCH | tee compile.log
    mv $NEMO_BUILD_DIR/cfgs/C1D_PAPA_FABM_BLD_SCRATCH/BLD/bin/nemo.exe $RUNDIR/
    echo "Done."
    
The script above will compile the model and move the nemo executable into your working directory.

Note that you might need to edit the architecture file depending on the machine you are compiling and running the model on. This will include compiler version, compiler flags and links to netCDF libraries. Here, we are pointing our compilation to arch-GCC_PMPC.fcm file (available within NEMO4.0-FABM repository) to compile on a typical PML workstation running Fedora Linux distribution and using a GNU Fortran compiler.

Step 2: Getting ready to run the model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now, you need to let NEMO know that you are running the simulation with ERSEM. Within your working directory, create a link to the desired configuration, e.g.:

  .. code-block:: bash
  
     ln -sf <ERSEMDIR>/testcases/fabm-ersem-15.06-L4-noben-docdyn-iop.yaml fabm.yaml
     # replace <ERSEMDIR> with the location of your ERSEM code.
     
ERSEM requires some external inputs, which we must provde. The following lines should be appended to your fabm.yaml file. Note that for simplicity we are using constant values here. Depending on the configuration, the list of external inputs will vary.

  .. code-block:: bash
  
       pco2a:
         model: horizontal_constant
         parameters:
           value: 400.
           standard_name: mole_fraction_of_carbon_dioxide_in_air
       ADY_0:
         model: horizontal_constant
         parameters:
           value: 1.0e-10
           standard_name: gelbstoff_absorption_satellite

    


