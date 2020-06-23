# locust_mc

locust_mc is a monte carlo simulation package of the receiver chain for Project 8.

### Instructions For Use:

Once locust is installed and added to your path (see Installing in Section #4), it runs and is typically configured on the command line as:
```
  $ LocustSim config=/path/to/config/LocustConfigFile.json
```
The output of the simulation is a file in egg format.  It is recommended to process the egg files with Katydid with an appropriate config file, for example as
```
  $ Katydid -c /path/to/config/katydid.json
```
After running Katydid, the processed fft spectra will have a span from 0 to the sampling frequency fs where fs/2 corresponds to DC.  For example, if fs=200 MHz, signals that are mixed down to an intermediate frequency of 50 MHz in Locust will appear at 100 MHz + 50 MHz = 150 MHz after Katydid.

### Tutorial

The [Config](https://github.com/project8/locust_mc/tree/develop/Config) directory contains instructions to guide the user through several short examples. 

Issues should be posted via [GitHub](https://github.com/project8/locust_mc/issues).

Installation
---------------

### Dependencies

**External**
 - Boost (www.boost.org) version 1.66 or higher (date_time, filesystem, program_options, system, thread)
 - CMake (www.cmake.org) version 3.15 or higher
 - FFTW3 (3.3 or newer)
 - G++ version 5.0 or higher (if compiling with GCC)
 - GSL (www.gnu.org/software/gsl) version 2.0 or higher
 - HDF5 (required by Monarch3 and for outputing to HDF5 files)
 - ROOT (www.cern.ch/root) version 5.24 or higher (6.x should work too)

**Optional Dependencies**
 - LibXml2 (xmlsoft.org)
 - Log4CXX (logging.apache.org/log4cxx)
 - MPI (www.open-mpi.org or mpich.org)
 - OpenCL (www.khronos.org/opencl), installation details depend on your system
 - OpenSSL (openssl.org) version 0.9.6 or higher
 - PETSc (mcs.anl.gov/petsc)
 - VTK (www.vtk.org) version 5.0 or higher
 - zlib (www.zlib.net)
 - pdflatex (for making the documentation; minimum version not known)
 - doxygen (for making the documentation; minimum version not known)(3.1 or better)

**Submodules** (included with locust; must be fetched via Git)
- [Monarch](https://github.com/project8/monarch)
  - [Scarab](https://github.com/project8/scarab)
    - yaml-cpp
- [Kassiopeia](https://github.com/project8/kassiopeia)


### Operating System Support

* Mac OS X (usually tested on OS X 10.12)
* Linux (usually tested on Ubuntu)


### Installing

The following steps will build locust from scratch.  In the terminal:

1. Clone the repository and make a build directory as recommended above. You will also have to initialize the submodules.
  ```
  $ git clone https://github.com/project8/locust_mc
  $ cd locust_mc
  $ git submodule update --init --recursive
  $ mkdir build
  ```

2. To configure the installation you can use cmake, ccmake, or cmake-gui.

  ```
  $ cd build
  $ ccmake ..
  ```

  You will be prompted to press [c] to configure, and the window will fill up with several options. 

  You should set the CMake variable `CMAKE_BUILD_TYPE` to either `RELEASE`, `STANDARD`, or `DEBUG` (default), in order
  of how much text output you would like (from least to most) and how much compiler optimization
  should be performed (from most to least).

  The install prefix is specified by the CMake variable `CMAKE_INSTALL_PREFIX`.
  The library, binaries, and header files will be installed in the
  lib, bin, and include subdirectories. The default install prefix is the
  build directory.

  After you've finished, if you've changed anything press [c] again to configure.  Then [g] to generate and exit.

3. Build and install.

  ```
  $ make install
  ```

  Or if you want to take advantage of parallel building to get things done faster:
  ```
  $ make -j2 install
  ```

  If the compiler runs into errors during the build, first check that you've updated the submodules and that you have all of the required dependencies installed (many are called "optional" on this page, but if you want to build without them you must also specify this in the cmake window). If you made a change to the dependencies or submodules, you may have to wipe the build directory and start again from step 1; simply writing `make install` again will not always work. 

4. Add the locust binaries to your $PATH to call executables directly from the command line.

    ```
    $ source /path/to/locust_mc/build/bin/kasperenv.sh 
    ```
  If you have another kassiopeia installation on your computer, you should ensure that locust is not using the environmental variables of the independent kassiopeia installation by removing this line from your .bashrc.


Directory Structure
-------------------

*  Config - Example configuration files
*  Documentation - Doxygen-based code documentation.
*  kassiopeia - Submodule for simulation of charged particle trajectories
*  monarch - Submodule library for file I/O in Project 8
*  Scarab - Submodule
*  Source - Locust source code
