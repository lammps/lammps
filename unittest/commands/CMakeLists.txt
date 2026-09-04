# -*- CMake -*- file for testing LAMMPS commands

add_executable(test_simple_commands test_simple_commands.cpp)
add_library(noplugin MODULE noplugin.cpp)
set_target_properties(noplugin PROPERTIES PREFIX "" SUFFIX ".so")

# tests for the plugin command require the PLUGIN package.
# With a STATIC LAMMPS library the plugins must be built without embedding a
# full copy of LAMMPS (loading such a plugin corrupts memory and aborts at
# program exit, see the comment in examples/plugins/CMakeLists.txt), which
# requires the COMPILE_ONLY generator expression from CMake 3.27; the test
# executable must then export its symbols so the loaded plugins can resolve
# them.  Windows is not affected and keeps linking the full library.
if(PKG_PLUGIN)
  if(BUILD_SHARED_LIBS OR (CMAKE_VERSION VERSION_GREATER_EQUAL 3.27)
     OR (CMAKE_SYSTEM_NAME STREQUAL "Windows"))
    add_subdirectory(${LAMMPS_DIR}/examples/plugins ${CMAKE_BINARY_DIR}/build-plugins)
    add_dependencies(test_simple_commands plugins noplugin)
    target_compile_definitions(test_simple_commands PRIVATE -DLMP_PLUGIN)
    if((NOT BUILD_SHARED_LIBS) AND (NOT (CMAKE_SYSTEM_NAME STREQUAL "Windows")))
      set_target_properties(test_simple_commands PROPERTIES ENABLE_EXPORTS TRUE)
    endif()
  else()
    message(WARNING "Skipping tests for the plugin command: they require "
      "CMake 3.27 or later or setting BUILD_SHARED_LIBS=on")
  endif()
endif()

target_link_libraries(test_simple_commands PRIVATE lammps GTest::GMock)
add_test(NAME SimpleCommands COMMAND test_simple_commands)

# the geturl test downloads from the network, which takes long and varies with
# the connection, so it is a separate executable with the "slow" label
add_executable(test_geturl_command test_geturl_command.cpp)
target_link_libraries(test_geturl_command PRIVATE lammps GTest::GMock)
add_test(NAME GeturlCommand COMMAND test_geturl_command)
set_tests_properties(GeturlCommand PROPERTIES LABELS "slow")
if(CMAKE_VERSION VERSION_LESS 3.20)
  set(LAMMPS_PLUGIN_BIN_DIR ${CMAKE_BINARY_DIR})
else()
  cmake_path(SET LAMMPS_PLUGIN_BIN_DIR NORMALIZE ${CMAKE_BINARY_DIR})
  if(CMAKE_CONFIG_TYPE)
    cmake_path(APPEND LAMMPS_PLUGIN_BIN_DIR ${CMAKE_CONFIG_TYPE})
  endif()
endif()
set_tests_properties(SimpleCommands PROPERTIES
                     ENVIRONMENT "LAMMPS_PLUGIN_DIR=${LAMMPS_PLUGIN_BIN_DIR}")
add_kokkos_test(NAME SimpleCommands COMMAND test_simple_commands
                ENVIRONMENT "LAMMPS_PLUGIN_DIR=${LAMMPS_PLUGIN_BIN_DIR}")

add_executable(test_lattice_region test_lattice_region.cpp)
target_link_libraries(test_lattice_region PRIVATE lammps GTest::GMock)
add_test(NAME LatticeRegion COMMAND test_lattice_region)
add_kokkos_test(NAME LatticeRegion COMMAND test_lattice_region)

add_executable(test_fix_spring_chunk test_fix_spring_chunk.cpp)
target_link_libraries(test_fix_spring_chunk PRIVATE lammps GTest::GMock)
add_test(NAME FixSpringChunk COMMAND test_fix_spring_chunk)
add_kokkos_test(NAME FixSpringChunk COMMAND test_fix_spring_chunk)

add_executable(test_kspace_auto_slab test_kspace_auto_slab.cpp)
target_compile_definitions(test_kspace_auto_slab PRIVATE -DTEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_kspace_auto_slab PRIVATE lammps GTest::GMock)
add_test(NAME KSpaceAutoSlab COMMAND test_kspace_auto_slab)
add_kokkos_test(NAME KSpaceAutoSlab COMMAND test_kspace_auto_slab)

add_executable(test_kspace_group_group_slab test_kspace_group_group_slab.cpp)
target_compile_definitions(test_kspace_group_group_slab PRIVATE -DTEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_kspace_group_group_slab PRIVATE lammps GTest::GMock)
add_test(NAME KSpaceGroupGroupSlab COMMAND test_kspace_group_group_slab)
add_kokkos_test(NAME KSpaceGroupGroupSlab COMMAND test_kspace_group_group_slab)

if(PKG_QMMM-XTB)
  add_executable(test_kspace_group_potential test_kspace_group_potential.cpp)
  target_include_directories(test_kspace_group_potential PRIVATE ${LAMMPS_SOURCE_DIR}/KSPACE
                                                                  ${LAMMPS_SOURCE_DIR}/QMMM-XTB)
  target_link_libraries(test_kspace_group_potential PRIVATE lammps GTest::GMock)
  if(FFT_SINGLE)
    target_compile_definitions(test_kspace_group_potential PRIVATE -DFFT_SINGLE)
  endif()
  add_test(NAME KSpaceGroupPotential COMMAND test_kspace_group_potential)
endif()

add_executable(test_groups test_groups.cpp)
target_link_libraries(test_groups PRIVATE lammps GTest::GMock)
add_test(NAME Groups COMMAND test_groups)
add_kokkos_test(NAME Groups COMMAND test_groups)

add_executable(test_regions test_regions.cpp)
target_link_libraries(test_regions PRIVATE lammps GTest::GMock)
add_test(NAME Regions COMMAND test_regions)
add_kokkos_test(NAME Regions COMMAND test_regions)

add_executable(test_info test_info.cpp)
target_link_libraries(test_info PRIVATE lammps GTest::GMock)
add_test(NAME InfoCommand COMMAND test_info)
add_kokkos_test(NAME InfoCommand COMMAND test_info)

add_executable(test_thermo test_thermo.cpp)
target_link_libraries(test_thermo PRIVATE lammps GTest::GMock)
add_test(NAME ThermoCommands COMMAND test_thermo)
add_kokkos_test(NAME ThermoCommands COMMAND test_thermo)

add_executable(test_delete_atoms test_delete_atoms.cpp)
target_link_libraries(test_delete_atoms PRIVATE lammps GTest::GMock)
add_test(NAME DeleteAtoms COMMAND test_delete_atoms)
add_kokkos_test(NAME DeleteAtoms COMMAND test_delete_atoms)

add_executable(test_set_property test_set_property.cpp)
target_link_libraries(test_set_property PRIVATE lammps GTest::GMock)
add_test(NAME SetProperty COMMAND test_set_property)
add_kokkos_test(NAME SetProperty COMMAND test_set_property)

add_executable(test_pair_ldd test_pair_ldd.cpp)
target_link_libraries(test_pair_ldd PRIVATE lammps GTest::GMock)
add_test(NAME PairLDD COMMAND test_pair_ldd)
add_kokkos_test(NAME PairLDD COMMAND test_pair_ldd)

add_executable(test_labelmap test_labelmap.cpp)
target_compile_definitions(test_labelmap PRIVATE -DTEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_labelmap PRIVATE lammps GTest::GMock)
add_test(NAME Labelmap COMMAND test_labelmap)
add_kokkos_test(NAME Labelmap COMMAND test_labelmap)

add_executable(test_compute_frenkel test_compute_frenkel.cpp)
target_link_libraries(test_compute_frenkel PRIVATE lammps GTest::GMock)
add_test(NAME ComputeFrenkel COMMAND test_compute_frenkel)
add_kokkos_test(NAME ComputeFrenkel COMMAND test_compute_frenkel)

add_executable(test_variables test_variables.cpp)
target_link_libraries(test_variables PRIVATE lammps GTest::GMock)
add_test(NAME Variables COMMAND test_variables)
add_kokkos_test(NAME Variables COMMAND test_variables)

add_executable(test_kim_commands test_kim_commands.cpp)
if(KIM_EXTRA_UNITTESTS)
  if(CURL_FOUND)
    target_compile_definitions(test_kim_commands PRIVATE -DKIM_EXTRA_UNITTESTS)
  else()
    message(FATAL_ERROR "CURL not found. Enabling KIM extra unit tests requires to have libcurl installed.")
  endif()
endif()
target_link_libraries(test_kim_commands PRIVATE lammps GTest::GMock)
add_test(NAME KimCommands COMMAND test_kim_commands)

add_executable(test_reset_atoms test_reset_atoms.cpp)
target_compile_definitions(test_reset_atoms PRIVATE -DTEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_reset_atoms PRIVATE lammps GTest::GMock)
add_test(NAME ResetAtoms COMMAND test_reset_atoms)
add_kokkos_test(NAME ResetAtoms COMMAND test_reset_atoms)

if(PKG_MOLECULE)
  add_executable(test_compute_global test_compute_global.cpp)
  target_compile_definitions(test_compute_global PRIVATE -DTEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(test_compute_global PRIVATE lammps GTest::GMock)
  add_test(NAME ComputeGlobal COMMAND test_compute_global)
  add_kokkos_test(NAME ComputeGlobal COMMAND test_compute_global)

  add_executable(test_compute_chunk test_compute_chunk.cpp)
  target_compile_definitions(test_compute_chunk PRIVATE -DTEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(test_compute_chunk PRIVATE lammps GTest::GMock)
  add_test(NAME ComputeChunk COMMAND test_compute_chunk)
  add_kokkos_test(NAME ComputeChunk COMMAND test_compute_chunk)

  add_executable(test_fix_halt test_fix_halt.cpp)
  target_compile_definitions(test_fix_halt PRIVATE -DTEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(test_fix_halt PRIVATE lammps GTest::GMock)
  add_test(NAME FixHalt COMMAND test_fix_halt)
  add_kokkos_test(NAME FixHalt COMMAND test_fix_halt)
endif()

if(PKG_PHONON)
  add_executable(test_phonon_commands test_phonon_commands.cpp)
  target_link_libraries(test_phonon_commands PRIVATE lammps GTest::GMock)
  add_test(NAME PhononCommands COMMAND test_phonon_commands)
  add_kokkos_test(NAME PhononCommands COMMAND test_phonon_commands)
endif()

if(PKG_REAXFF)
  add_executable(test_reaxff_output test_reaxff_output.cpp)
  target_link_libraries(test_reaxff_output PRIVATE lammps GTest::GMock)
  add_test(NAME ReaxFFOutput COMMAND test_reaxff_output)
  set_tests_properties(ReaxFFOutput PROPERTIES
                       ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR}")
  add_kokkos_test(NAME ReaxFFOutput COMMAND test_reaxff_output
                  ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR}")
endif()

add_executable(test_mpi_load_balancing test_mpi_load_balancing.cpp)
target_link_libraries(test_mpi_load_balancing PRIVATE lammps GTest::GMock)
add_mpi_test(NAME MPILoadBalancing NUM_PROCS 4 COMMAND $<TARGET_FILE:test_mpi_load_balancing>)

add_executable(test_restart_settings test_restart_settings.cpp)
target_link_libraries(test_restart_settings PRIVATE lammps GTest::GMock)
add_test(NAME RestartSettings COMMAND test_restart_settings)
add_kokkos_test(NAME RestartSettings COMMAND test_restart_settings)

add_executable(test_input_validation test_input_validation.cpp)
target_link_libraries(test_input_validation PRIVATE lammps GTest::GMock)
add_test(NAME InputValidation COMMAND test_input_validation)
add_kokkos_test(NAME InputValidation COMMAND test_input_validation)
