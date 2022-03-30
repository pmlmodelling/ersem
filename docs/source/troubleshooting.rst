.. _trouble:

###############################
Tips and tricks/troubleshooting
###############################

General:
~~~~~~~~

-  cmake autodetects a suitable Fortran compiler. If a Fortran compiler
   is not found or you are unhappy with the compiler cmake selected, you
   can override this by providing cmake with
   ``-DCMAKE_Fortran_COMPILER=<COMPILER_EXECUTABLE>``. Note: if you
   change ``CMAKE_Fortran_COMPILER`` (or provide it for the first time,
   while you previously ran cmake without), you need to empty your build
   directory first!
-  When building GOTM or FABM’s 0d driver, cmake will try to auto-detect
   NetCDF using ``nf-config``. If nf-config is not present on your
   system, you’ll need to provide cmake with the path(s) to the NetCDF
   include directories (``-DNetCDF_INCLUDE_DIRS=<PATH>``) and the
   path(s) to the NetCDF libraries (``-DNetCDF_LIBRARIES=<PATH>``). If
   you need to provide multiple paths to these variables, the individual
   paths should be separated by semi-colons. A common reason to use
   this: On some systems, ``nf-config`` does not detect where
   ``netcdf.mod`` is installed, which means you have to tell cmake by
   adding ``-DNetCDF_INCLUDE_DIRS=<netcdf.mod_DIR>``. For instance, when
   using gfortran on Fedora, ``<netcdf.mod_DIR>`` can be
   ``/usr/lib64/gfortran/modules``. In that case you usually do not need
   to provide ``-DNetCDF_LIBRARIES=<PATH>``. Also, on some systems (like
   in the current LTS release of Ubuntu), ``nf-config`` is not included
   in the NetCDF packages. In this case, you can use ``nc-config``:
   specify ``-DNetCDF_CONFIG_EXECUTABLE=<PATH_TO_nc-config>`` when
   calling cmake, or create a link to ``nc-config`` somewhere in your
   default path, to get auto-detection working.
-  By default, cmake will select the “release” build type, which creates
   an executable without debugging information. If you want to compile
   in debug mode, specify ``-DCMAKE_BUILD_TYPE=debug`` in the call to
   cmake. When using gfortran, you may also want to add
   ``-DCMAKE_Fortran_FLAGS_DEBUG=-fcheck=bounds`` to catch out-of-bounds
   array access.
-  To see the settings you specified when you ran cmake, and to
   selectively make changes to these setttings, you can run ``ccmake .``
   in your build directory.
-  To speed up compilation, you can perform a parallel build by
   providing the ``-j N`` switch to ``make`` (not ``cmake``!), with
   ``N`` being the number of cores you want to use.

Specific platforms:

-  `ARCHER <https://www.archer.ac.uk>`__: make sure to first load the
   cmake and cray-netcdf modules (``module load cmake cray-netcdf``).
   These enable cmake and autodetection of NetCDF paths, respectively.
   Provide ``-DCMAKE_Fortran_COMPILER=ftn`` to cmake to make sure cmake
   uses Cray’s compiler wrapper script for Fortran compilation.
   Compilation should succeed with all three programming environments
   (intel, gnu, cray).
-  Ubuntu LTS: ``nf-config`` is not included in the netcdf packages, but
   ``nc-config`` is. Specify
   ``-DNetCDF_CONFIG_EXECUTABLE=<PATH_TO_nc-config>`` when calling
   cmake, or create a link from ``nf-config`` to ``nc-config`` somewhere
   in your default path to get auto-detection of NetCDF paths working.
-  Fedora: ``nf-config`` does not detect where ``netcdf.mod`` is
   installed (typically in ``/usr/lib64/gfortran/modules``). This means
   you have to tell cmake by adding
   ``-DNetCDF_INCLUDE_DIRS=/usr/lib64/gfortran/modules``.


Migrating from an earlier ERSEM version
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you were using an earlier release of ERSEM, you can update your old
ERSEM configuration (``fabm.yaml``) to the latest by running the Python
script ``<SOURCEDIR>/ersem/testcases/update.py`` with one argument: the
path to your old fabm.yaml file.

This script requires a recent version of Python 2 or 3 and the `PyYAML
package <https://pyyaml.org/wiki/PyYAML>`__. Ideally, you also have the
latest version of `the Python front end to
FABM-ERSEM <#python-front-end>`__ installed; this enables the update
script to clean up the yaml file and add documentation to it.
