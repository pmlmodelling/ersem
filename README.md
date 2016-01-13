# Obtaining the code and building 
## Linux

First obtain the FABM source code:

    git clone git://git.code.sf.net/p/fabm/code <FABMDIR>

(Replace `<FABMDIR>` with the directory with the FABM code, e.g., ~/fabm-git.)

Obtain the ERSEM code for FABM:

    git clone git@gitlab.ecosystem-modelling.pml.ac.uk:edge/ersem.git <ERSEMDIR>

(Replace `<ERSEMDIR>` with the directory with the ERSEM code, e.g., ~/ersem-git.) Note that for this to work, you have to provide [the PML GitLab server](https://gitlab.ecosystem-modelling.pml.ac.uk/profile/keys) with your public SSH key.

FABM and ERSEM use object-oriented Fortran and therefore require a recent Fortran compiler, such as Intel Fortran 12.1 or higher and gfortran 4.7 or higher. Compilation is regularly tested with Intel Fortran 12.1, 13.0 and 14.0. as well as gfortran 4.7, 4.8 and 4.9.

FABM and ERSEM use a platform-independent build system based on [cmake](http://www.cmake.org). You'll need version 2.8.8 or higher. First check whether you have that installed: execute `cmake --version` on the command line.

### GOTM + FABM + ERSEM

First obtain the latest (developers') version of the GOTM code from its git repository:

    git clone git://git.code.sf.net/p/gotm/code gotm-git

To build GOTM with FABM support, create a build directory, call cmake to generate makefiles, and make to compile and install. For instance:

    mkdir -p ~/build/gotm && cd ~/build/gotm
    cmake <GOTMDIR>/src -DGOTM_USE_FABM=ON -DFABM_BASE=<FABMDIR> -DFABM_ERSEM_BASE=<ERSEMDIR>
    make install

In the above, replace `<GOTMDIR>` with the directory with the GOTM source code, e.g., ~/gotm-git if you executed `git clone` in you home directory. Also, replace `<FABMDIR>` with the directory with the FABM code, e.g., ~/fabm-git and `<ERSEMDIR>` with the directory with the ERSEM code, e.g., ~/ersem-git. If you experience issues related to NetCDF, see [tips and tricks/troubleshooting](#tips-and-tricks-troubleshooting).

Now you should have a GOTM executable with FABM and ERSEM support at `~/local/gotm/bin/gotm`.

It is good practice to keep up to date with the latest code from the ERSEM, FABM and GOTM repositories by regularly executing `git pull` in a directory of each repository.

If either the ERSEM, FABM or GOTM source codes change (e.g., because changes you made to the code yourself, or after `git pull`), you will need to recompile. This does NOT require rerunning cmake. Instead, you need to return to the build directory and rerun `make install`. For instance `cd ~/build/gotm && make install`.

### FABM-ERSEM 0d ("aquarium")

The 0d driver allows you to run FABM models in a "well-mixed box", under arbitrary (time-varying) environmental forcing.

To build the 0d driver, you need to create a directory to build the code in, call `cmake` to generate makefiles, and call `make` to compile and install the FABM library. Usually, the following suffices for this:

    mkdir -p ~/build/fabm-0d && cd ~/build/fabm-0d
    cmake <FABMDIR>/src/drivers/0d -DGOTM_BASE=<GOTMDIR> -DFABM_ERSEM_BASE=<ERSEMDIR>
    make install

In the above, replace `<ERSEMDIR>` with the path to directory with the ERSEM source code (e.g., ~/ersem-git), `<FABMDIR>` with the path to directory with the FABM source code (e.g., ~/fabm-git), and `<GOTMDIR>` with the path to directory with the GOTM source code (e.g., ~/gotm-git). The latter is needed because the 0d driver uses GOTM routines for input, output, time integration, etc. If you experience issues related to NetCDF, see [tips and tricks/troubleshooting](#tips-and-tricks-troubleshooting).

This will give you an executable at `~/local/fabm/0d/bin/fabm0d`.

To use the driver, you need a configuration file can run.nml, which you could take from [<FABMDIR>/testcases/0d/run.nml](../blob/master/testcases/0d/run.nml) In addition, you need a file with forcing data for surface shortwave radiation, temperature and salinity, e.g., [<FABMDIR>/testcases/0d/env_nns_annual.dat](../blob/master/testcases/0d/env_nns_annual.dat). The local name of this file is configured in run.nml, variable env_file. Finally, you can provide additional forcing with a `input.yaml` file, which can contain entries such as

    wind_speed:
      file: wind.dat
      column: 1
      scale_factor: 1.0
    mass_concentration_of_silt:
      constant_value: 400.
    mole_fraction_of_carbon_dioxide_in_air:
      constant_value: 401.0
    bottom_stress:
      constant_value: 0.0

This specifies that wind speed must be read from file wind.dat (first column, no scaling), the mass concentration of silt must be set to a constant value of 400 mg/m3, atmospheric pCO2 to 401 ppm, and bottom stress to 0 Pa. The "column" and "scale_factor" attributes in the case of wind_speed are optional. They default to 1 and 1.0, respectively.

### Python front-end

The Python front-end can be used for enumerating model metadata (e.g., variable information), to add documentation to fabm.yaml files, to retrieve sources-sinks and diagnostics while manipulating the environment, and to show a configuration GUI [the last is still under development].

To build the Python driver, you need to create a directory to build the code in, call `cmake` to generate makefiles, and call `make` to compile and install the FABM library. Usually, the following suffices for this:

    mkdir -p ~/build/fabm-python && cd ~/build/fabm-python
    cmake <FABMDIR>/src/drivers/python -DFABM_ERSEM_BASE=<ERSEMDIR>
    make install

In the above, replace `<FABMDIR>` with the directory with the FABM source code, e.g., ~/fabm-git and `<ERSEMMDIR>` with the directory with the ERSEM source code, e.g., ~/ersem-git.

This will install the python-fabm module in a directory that is automatically looked in by Python. Typically, this is ~/.local/lib/python2.X/site-packages (with X being Python's minor version number). This means you can now just open Python and enter `import pyfabm` to load FABM's Python module.

Example scripts that use the Python front-end can be found in `<FABMDIR>/testcases/python`. For instance, use `<FABMDIR>/testcases/python/fabm_complete_yaml.py` to add documentation to a fabm.yaml file.

### NEMO+FABM+ERSEM

FABM needs to be compiled separately before the compilation of NEMO.
Usually, the following suffices for this:

   mkdir -p ~/build/nemo && cd ~/build/nemo
   cmake <FABMDIR>/src/ -DFABM_HOST=nemo -DFABM_ERSEM_BASE=<ERSEMDIR>
   make install

In the above, replace `<FABMDIR>` with the directory with the FABM source code, e.g., ~/fabm-git and `<ERSEMMDIR>` with the directory with the ERSEM source code, e.g., ~/ersem-git.

This will create the library in the standard folder ~/local/fabm/nemo/lib where NEMO-FABM will look for linking to NEMO.

### Tips and tricks/troubleshooting

General:

* cmake autodetects a suitable Fortran compiler. If a Fortran compiler is not found or you are unhappy with the compiler cmake selected, you can override this by providing cmake with `-DCMAKE_Fortran_COMPILER=<COMPILER_EXECUTABLE>`. Note: if you change CMAKE_Fortran_COMPILER (or provide it for the first time, while you previously ran cmake without), you need to clean (remove everything in) your build directory first!
* When building GOTM or FABM's 0d driver, cmake will try to auto-detect NetCDF using `nf-config`. If nf-config is not present on your system, you'll need to provide cmake with the path(s) to the NetCDF include directories (`-DNetCDF_INCLUDE_DIRS=<PATH>`) and the path(s) to the NetCDF libraries (`-DNetCDF_LIBRARIES=<PATH>`). If you need to provide multiple paths to these variables, the individual paths should be separated by semi-colons. A common reason to use this: On some systems, nf-config does not detect where netcdf.mod is installed, which means you have to tell cmake by adding `-DNetCDF_INCLUDE_DIRS=<PATH_TO_netcdf.mod>`. For instance, when using gfortran on Fedora, `<PATH_TO_netcdf.mod>` can be `/usr/lib64/gfortran/modules`. In that case you usually do not need to provide `-DNetCDF_LIBRARIES=<PATH>`. Also, on some systems (like in the current LTS release of Ubuntu), nf-config is not included in the netcdf packages. In these case, specify `-DNetCDF_CONFIG_EXECUTABLE=<PATH_TO_nc-config>` when calling cmake, or create a link to nc-config somewhere in your default path, to get auto-detection working.
* By default, cmake will select the "release" build type, which creates an executable without debugging information. If you want to compile in debug mode, specify `-DCMAKE_BUILD_TYPE=debug` in the call to cmake. When using gfortran, you may also want to add `-DCMAKE_Fortran_FLAGS_DEBUG=-fcheck=bounds` to catch out-of-bounds array access.
* To see the settings you specified when you ran cmake, and to selectively make changes to these setttings, you can run `ccmake .` in your build directory.
* To speed up compilation, you can perform a parallel build by providing the -j N switch to make (not cmake!), with N being the number of cores you want to use.

Specific platforms:

* [ARCHER](http://www.archer.ac.uk): make sure to first load the cmake and cray-netcdf modules (`module load cmake cray-netcdf`). These enable cmake and autodetection of NetCDF paths, respectively. Provide `-DCMAKE_Fortran_COMPILER=ftn` to cmake to make sure cmake uses Cray's compiler wrapper script for Fortran compilation. Compilation should succeed with the GNU (`module load PrgEnv-gnu`) and the Intel programming environments (`module load PrgEnv-intel`). Compilation with the Cray compiler (`module load PrgEnv-cray`) currently fails due to compiler bugs related to Fortran 2003 support.
* Ubuntu LTS: nf-config is not included in the netcdf packages, but nc-config is. Specify `-DNetCDF_CONFIG_EXECUTABLE=<PATH_TO_nc-config>` when calling cmake, or create a link from nf-config to nc-config somewhere in your default path to get auto-detection of NetCDF paths working.
* Fedora: nf-config does not detect where netcdf.mod is installed (typically in /usr/lib64/gfortran/modules). This means you have to tell cmake by adding `-DNetCDF_INCLUDE_DIRS=/usr/lib64/gfortran/modules`.

## Windows

Note: below are quick-start instructions tailored to Visual Studio. Further information, including instruction for building with the free MinGW compiler, can be found on [the FABM wiki](http://fabm.net/wiki), section "Building and installing".

NB tested with Visual Studio 2010 in combination with Intel Visual Fortran (version 12.1 or higher required).

To obtain the source code of GOTM, FABM and ERSEM, you need a git client. First install [the Windows version of git itself](http://git-scm.com/download/win). For convenience, you can install the graphical git client [TortoiseGit](https://tortoisegit.org/) on top; the following instructions assume you have installed TortoiseGit.

After these two program are installed, you can obtain the code by right-clicking in Windows Explorer within a directory where you want the source code directories, and choosing "Git Clone...". In the window that appears, set the URL, check the target directory and click OK. Do this for the following URLs:

* GOTM: git://git.code.sf.net/p/gotm/code (suggested target directory: gotm-git)
* FABM: git://git.code.sf.net/p/fabm/code (suggested target directory: fabm-git)
* ERSEM: git@gitlab.ecosystem-modelling.pml.ac.uk:edge/ersem.git (suggested target directory: ersem-git). For this repository, you also need to provide your private SSH key, which must match the public key that you provided on [the PML GitLab site](https://gitlab.ecosystem-modelling.pml.ac.uk/profile/keys). To do so, check "Load Putty key" and set its path to the file with your private key. For creating private/public keys, we suggest using [PuTTYgen](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html).

To compile the code, you need [CMake](http://www.cmake.org/). If CMake is installed, open "CMake (cmake-gui)", specify GOTM's `src` directory for "Where is the source code", choose a directory of your choice (but outside the source directory!) for "Where to build the binaries", and click Configure. Choose the generator that matches your Visual Studio version (avoid the IA64 and Win64 options) and click Finish. Several configuration options will appear, among which `FABM_BASE`. This option must be set to the directory with the FABM source code (the root directory, not the `src` subdirectory). After doing so, click "Configure". Then a new option `FABM_ERSEM_BASE` will appear, which must be set to the path to the directory with ERSEM source code (the root directory, not the `src` subdirectory). Keep clicking "Configure" until there are no red-coloured options left. Then press "Generate". This will create a `gotm.sln` Visual Studio solution in the specified build directory, which you can now open with Visual Studio.

Note that it is good practice to keep up to date with the latest code from the GOTM, FABM and ERSEM repositories by regularly right-clicking the repository directory, choosing "Git Sync...", and clicking the "Pull" button in the window that then appears.