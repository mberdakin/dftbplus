#
# Global architecture independent build settings
#

# CMAKE_BUILD_TYPE is commented out in order to allow for multi-configuration builds. It will
# automatically default to RelWithDebInfo if used in a single configuration build. Uncomment or
# override it only if you want a non-default single configuration build.
#
#set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type (Release|RelWithDebInfo|Debug|MinSizeRel)")

option(WITH_OMP "Whether OpenMP thread parallisation should be enabled" TRUE)

option(WITH_MPI "Whether DFTB+ should support MPI-parallelism" FALSE)
# If you build an MPI-parallised binary, consider to set WITH_OMP (OpenMP thread parallelisaton) to
# FALSE unless you want hybrid parallelisation (for experts only).

option(WITH_ELSI "Whether DFTB+ with MPI-parallelism should use the ELSI libraries" FALSE)
# Works only with MPI-parallel build.

option(WITH_GPU "Whether DFTB+ should support GPU-acceleration via the MAGMA-library" FALSE)

option(WITH_TRANSPORT "Whether transport via libNEGF should be included." FALSE)
# Works only when building static libraries (see option BUILD_SHARED_LIBS)

option(WITH_SOCKETS "Whether socket communication should be allowed for" FALSE)

option(WITH_ARPACK "Whether the ARPACK library should be included (needed for TD-DFTB)" FALSE)
# Works only with non-MPI (serial) build, needed for Casida linear response

option(WITH_DFTD3 "Whether the DFTD3 library should be included" FALSE)
# NOTE: Due to the license of the DFTD3 library, the combined code must be distributed under the
# GPLv3 license (as opposed to the LGPLv3 license of the DFTB+ package)

option(WITH_MBD "Whether DFTB+ should be built with many-body-dispersion support" FALSE)

option(WITH_PLUMED "Whether metadynamics via the PLUMED2 library should be allowed for" FALSE)

option(WITH_API "Whether public API should be included and the DFTB+ library installed" TRUE)
# Turn this on, if you want to use the DFTB+ library to integrate DFTB+ into other software
# packages. (Otherwise only a stripped down version of the library without the public API is built.)
# This will also install necessary include and module files and further libraries needed to link the
# DFTB+ library.

option(WITH_PYTHON "Whether the Python components of DFTB+ should be tested and installed" TRUE)

option(BUILD_SHARED_LIBS "Whether the libraries built should be shared" FALSE)
# Turn this on, if the DFTB+ library (and other compiled libraries) should be shared libraries and
# dynamically linked to their applications. This results in smaller applications, but the libraries
# must be present at run-time (and the correct LD_LIBRARY_PATH environment variable must be set, so
# that they can be found by the operating system). If you want use the DFTB+ library from other
# software packages (see WITH_API option above), they may also require a shared library (e.g.
# calling DFTB+ functions from Python or Julia).


#
# Test environment settings
#
set(TEST_MPI_PROCS "1" CACHE STRING "Nr. of MPI processes used for testing")

set(TEST_OMP_THREADS "1" CACHE STRING "Nr. of OpenMP-threads used for testing")

# Command line used to launch the test code.
# The escaped variables (\${VARIABLE}) will be substituted by the corresponding CMake variables.
if(WITH_MPI)
  set(TEST_RUNNER_TEMPLATE "env OMP_NUM_THREADS=\${TEST_OMP_THREADS} mpiexec -n \${TEST_MPI_PROCS}"
    CACHE STRING "How to run the tests")
else()
  set(TEST_RUNNER_TEMPLATE "env OMP_NUM_THREADS=\${TEST_OMP_THREADS}" CACHE STRING
    "How to run the tests")
  set(MODES_RUNNER_TEMPLATE "env OMP_NUM_THREADS=\${TEST_OMP_THREADS}" CACHE STRING
    "How to run the modes code for tests")
endif()


#
# Installation options
#
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/_install" CACHE STRING
  "Directory to install the compiled code into")

set(INSTALL_INCLUDEDIR "dftbplus" CACHE PATH
  "Name of the project specific sub-folder within the install folder for include files")

set(INSTALL_MODULEDIR "${INSTALL_INCLUDEDIR}/modfiles" CACHE PATH
  "Installation directory for Fortran module files (within the install folder for include files)")

set(PKGCONFIG_LANGUAGE "Fortran" CACHE STRING
  "Compiler and Linker language to assume when creating the pkg-config export file (C or Fortran)")
# The pkg-config export file (lib/pkgconfig/dftbplus.pc) contains the compiler and linker options
# needed to link the DFTB+ library to an application. (It can be queried with the pkg-config tool.)
# Depending on the language setting ("C" or "Fortran") you would get the flags for the case of using
# that compiler for the linking.


#
# Whether to build external dependencies instead of finding and using installed versions of those.
# When enabled for a given dependency, its source code must be present in the corresponding
# external/DEPENDENCY/orig/ folder.
#
option(BUILD_MPIFX "Whether to build MpiFx as part of the build process" TRUE)

option(BUILD_SCALAPACKFX "Whether to build ScalapackFx as part of the build process" TRUE)

option(BUILD_MBD "Whether to build LibMbd as part of the build process" TRUE)
