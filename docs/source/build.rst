.. _build:

Building ERSEM
==============

Linux
~~~~~

First go to a directory where you want to place the source codes. This
will be referred to as ``<SOURCEDIR>`` below.

Now get the ERSEM, FABM and GOTM source codes:

.. code-block:: bash

   git clone https://github.com/pmlmodelling/ersem.git
   git clone https://github.com/fabm-model/fabm.git
   git clone --recurse-submodules https://github.com/gotm-model/code.git gotm

This locally creates ``gotm``, ``fabm``, ``ersem`` directories with source
code.

Note: the above gets you ERSEM’s public stable release. Developers can
check out the the developers’ version by substituting
``git@gitlab.ecosystem-modelling.pml.ac.uk:edge/ersem.git`` for
``https://github.com/pmlmodelling/ersem.git``.

To build the code, you will need: \* `a recent Fortran
compiler <https://github.com/fabm-model/fabm/wiki/Building-and-installing#supported-compilers>`__
\* `cmake <https://www.cmake.org>`__ 3.0 or higher. First check whether
you have that installed: run ``cmake --version`` on the command line.

FABM0d-GOTM-ERSEM
~~~~~~~~~~~~~~~~~

FABM’s 0d driver allows you to run biogeochemical models in a
"well-mixed box" under arbitrary (time-varying) environmental forcing.

To install FABMod-GOTM-ERSEM we suggest you use the following script below

.. literalinclude:: ../../github-actions/fabm0d-gotm-ersem-build.sh
    :language: bash
    :linenos:

PyFABM-ERSEM
~~~~~~~~~~~~

The Python front-end (``pyfabm``) can be used for enumerating model
metadata (e.g., variable information), to add documentation to fabm.yaml
files, or to retrieve sources-sinks and diagnostics while manipulating
the environment.

To install PyFABM-ERSEM we suggest you use the following script below

.. literalinclude:: ../../github-actions/pyfabm-ersem-build.sh
    :language: bash
    :linenos:

GOTM-FABM-ERSEM
~~~~~~~~~~~~~~~

To install GOTM-FABM-ERSEM we suggest you use the following script below

.. literalinclude:: ../../github-actions/gotm-fabm-ersem-build.sh
    :language: bash
    :linenos:


If you experience NetCDF issues when running ``make install``, see `tips
and tricks/troubleshooting <#tips-and-trickstroubleshooting>`__.

Now you should have a GOTM executable with FABM and ERSEM support at
``~/local/gotm/bin/gotm``.

Tips and tricks/troubleshooting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

General:

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

Windows
-------

Note: below are quick-start instructions tailored to Visual Studio.
Further information, including instruction for building with the free
MinGW compiler, can be found on `the FABM
wiki <http://fabm.net/wiki>`__, section “Building and installing”.

NB this was tested with Visual Studio 2017 in combination with Intel
Visual Fortran 19.

To obtain the source code of GOTM, FABM and ERSEM, you need a git
client. First install `the Windows version of git
itself <https://git-scm.com/download/win>`__. For convenience, you can
install the graphical git client
`TortoiseGit <https://tortoisegit.org/>`__ on top; the following
instructions assume you have installed TortoiseGit.

After these two program are installed, you can obtain the code by
right-clicking in Windows Explorer within a directory where you want the
source code directories, and choosing “Git Clone…”. In the window that
appears, set the URL, ensure “recursive” is checked and click OK. Do
this for the following URLs:

-  ERSEM: https://github.com/pmlmodelling/ersem.git
-  FABM: https://github.com/fabm-model/fabm.git
-  GOTM: https://github.com/gotm-model/code.git (also change ``code`` to
   ``gotm`` for the target directory)

To compile the code, you need `CMake <https://www.cmake.org/>`__. If
CMake is installed, open “CMake (cmake-gui)”, specify the directory with
GOTM code for “Where is the source code”, choose an empty directory of
your choice for “Where to build the binaries”, and click Configure.
Choose the generator that matches your Visual Studio version and click
Finish. Several configuration options will appear, among which
``FABM_BASE``. This option must be set to the directory with the FABM
source code. After doing so, click “Configure”. Then a new option
``FABM_ERSEM_BASE`` will appear, which must be set to the path to the
directory with ERSEM source code. Keep clicking “Configure” until there
are no red-coloured options left. Then press “Generate”. This will
create a ``gotm.sln`` Visual Studio solution in the specified build
directory, which you can now open with Visual Studio.

It is good practice to keep up to date with the latest code from the
GOTM, FABM and ERSEM repositories by regularly right-clicking each
repository directory, choosing “Git Sync…”, and clicking the “Pull”
button in the window that then appears.

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

.. |DOI| image:: https://zenodo.org/badge/302390544.svg
   :target: https://zenodo.org/badge/latestdoi/302390544
