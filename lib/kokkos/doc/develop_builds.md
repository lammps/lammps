
# Places to build options: architecture, device, advanced options, cuda options

These are the files that need to be updated when a new architecture or device is
added:

  + generate_makefile.bash
      * Interface for makefile system
  + cmake/kokkos_options.cmake
      * Interface for cmake system
  + Makefile.kokkos
      * Main logic for build (make and cmake) and defines (KokkosCore_config.h)
  + core/unit_test/UnitTestConfig.make
      * Unit test for Makefile.kokkos

In general, an architecture is going to be from on of these platforms:
  + AMD
  + ARM
  + IBM
  + Intel
  + Intel Xeon Phi
  + NVIDIA
Although not strictly necessary, it is helpful to keep things organized by
grouping by platform.

### generate_makefile.sh

The bash code does not do any error checking on the `--arch=`  or `--device=`
arguments thus strictly speaking you do not *need* to do anything to add a 
device or architecture; however, you should add it to the help menu.  For the
archictectures, please group by one of the platforms listed above.


### cmake/kokkos_options.cmake and cmake/kokkos_settings.cmake

The options for the CMake build system are: `-DKOKKOS_HOST_ARCH:STRING=` and
`-DKOKKOS_ENABLE_<device>:BOOL=`.  Although any string can be passed into
KOKKOS_HOST_ARCH option, it is checked against an accepted list.  Likewise, the
KOKKOS_ENABLE_<device> must have the option added AND it is formed using the
list. Thus: 
  + A new architecture should be added to the KOKKOS_HOST_ARCH_LIST variable.
  + A new device should be added to the KOKKOS_DEVICES_LIST variable **AND** a
    KOKKOS_ENABLE_<newdevice> option specified (see KOKKOS_ENABLE_CUDA for
    example).
  + A new device should be added to the KOKKOS_DEVICES_LIST variable **AND** a

The translation from option to the `KOKKOS_SETTINGS` is done in
`kokkos_settings.cmake`.  This translation is automated for some types if you ad
to the list, but for others, it may need to be hand coded. 


### Makefile.kokkos

This is the main coding used by both the make and cmake system for defining
the sources (generated makefile and cmake snippets by `core/src/Makefile`), for
setting the defines in KokkosCore_config.h, and defining various internal
variables.  To understand how to add to this file, you should work closely with
the Kokkos development team.


### core/unit_test/UnitTestConfig.make

This file is used to check the build system in a platform-independent way.  It
works by looping over available architectures and devices; thus, you should add
your new architecure to KOKKOS_ARCH_OPTIONS and your new device to 
KOKKOS_DEVICE_OPTIONS to be tested.  The build system tests work by grepping the
generated build files (automatically).  The header file tests work by diffing
the generated file with results that are stored in
`core/unit_tests/config/results` (namespaced by ARCH_DEVICE_).  Thus, you will
need to add accepted results to this directory for diffing.

The CMake build system is also tested in `core/unit_tests/config/cmaketest`.
Because it uses cmake/kokkos_options.cmake, it already has the tests to loop
over.  It is diffed with the same files that the build system is tested with.
Thus, if you are consistent in all of the files listed, the unit tests should
pass automatically.
