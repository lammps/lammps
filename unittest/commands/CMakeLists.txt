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
if(CMAKE_VERSION VERSION_LESS 3.20)
  set_tests_properties(SimpleCommands PROPERTIES
                       ENVIRONMENT "LAMMPS_PLUGIN_DIR=${CMAKE_BINARY_DIR}")
else()
  cmake_path(SET LAMMPS_PLUGIN_BIN_DIR NORMALIZE ${CMAKE_BINARY_DIR})
  if(CMAKE_CONFIG_TYPE)
    cmake_path(APPEND LAMMPS_PLUGIN_BIN_DIR ${CMAKE_CONFIG_TYPE})
  endif()
  set_tests_properties(SimpleCommands PROPERTIES
                       ENVIRONMENT "LAMMPS_PLUGIN_DIR=${LAMMPS_PLUGIN_BIN_DIR}")
endif()

add_executable(test_lattice_region test_lattice_region.cpp)
target_link_libraries(test_lattice_region PRIVATE lammps GTest::GMock)
add_test(NAME LatticeRegion COMMAND test_lattice_region)

add_executable(test_kspace_auto_slab test_kspace_auto_slab.cpp)
target_compile_definitions(test_kspace_auto_slab PRIVATE -DTEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_kspace_auto_slab PRIVATE lammps GTest::GMock)
add_test(NAME KSpaceAutoSlab COMMAND test_kspace_auto_slab)

add_executable(test_kspace_group_group_slab test_kspace_group_group_slab.cpp)
target_compile_definitions(test_kspace_group_group_slab PRIVATE -DTEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_kspace_group_group_slab PRIVATE lammps GTest::GMock)
add_test(NAME KSpaceGroupGroupSlab COMMAND test_kspace_group_group_slab)

add_executable(test_groups test_groups.cpp)
target_link_libraries(test_groups PRIVATE lammps GTest::GMock)
add_test(NAME Groups COMMAND test_groups)

if(PKG_EXTRA-FIX AND PKG_EXTRA-COMPUTE)
  add_executable(test_fix_uvt test_fix_uvt.cpp)
  target_link_libraries(test_fix_uvt PRIVATE lammps GTest::GMock)
  add_test(NAME FixUVT COMMAND test_fix_uvt)

  add_executable(test_compute_temp_uvt test_compute_temp_uvt.cpp)
  target_link_libraries(test_compute_temp_uvt PRIVATE lammps GTest::GMock)
  add_test(NAME ComputeTempUVT COMMAND test_compute_temp_uvt)
  add_mpi_test(NAME ComputeTempUVTMPI NUM_PROCS 2
               COMMAND $<TARGET_FILE:test_compute_temp_uvt>)
endif()

if(PKG_REPLICA)
  add_executable(test_fix_pimd_nve test_fix_pimd_nve.cpp)
  target_link_libraries(test_fix_pimd_nve PRIVATE lammps GTest::GMock)
  add_test(NAME FixPIMDNVE COMMAND test_fix_pimd_nve
           --gtest_filter=FixPIMDNVESerialTest.*:PIMDStyles/FixPIMDIntegrationCompatibilityTest.*)
  string(CONCAT pimd_nve_partition_filter
         "--gtest_filter=FixPIMDNVEMPI.FreeRingPolymerMatchesAnalyticalMotion:"
         "FixPIMDNVEMPI.PIMDPartitionedRunSupportsOneRankPerBead:"
         "FixPIMDNVEMPI.ConservesTotalEnergyOverShortRun:"
         "FixPIMDNVEMPI.RemoveCOMRespectsRepresentationAndGroup")
  add_mpi_test(NAME FixPIMDNVEPartition NUM_PROCS 2
               COMMAND $<TARGET_FILE:test_fix_pimd_nve> "${pimd_nve_partition_filter}")
  add_mpi_test(NAME FixPIMDNVETagCommunication NUM_PROCS 6
               COMMAND $<TARGET_FILE:test_fix_pimd_nve>
                       --gtest_filter=FixPIMDNVEMPI.CommunicationMatchesTagsAcrossDifferentDecompositions)
  add_mpi_test(NAME FixPIMDNVEMultiRankPerBead NUM_PROCS 4
               COMMAND $<TARGET_FILE:test_fix_pimd_nve>
                       --gtest_filter=FixPIMDNVEMPI.MultiRankPerBeadRunProducesFiniteVector)

  add_executable(test_fix_pimd_nvt test_fix_pimd_nvt.cpp)
  target_link_libraries(test_fix_pimd_nvt PRIVATE lammps GTest::GMock)
  add_test(NAME FixPIMDNVTSerial
           COMMAND test_fix_pimd_nvt --gtest_filter=FixPIMDNVTSerialTest.*)
  add_mpi_test(NAME FixPIMDNVTMPI NUM_PROCS 2
               COMMAND $<TARGET_FILE:test_fix_pimd_nvt>
                       --gtest_filter=FixPIMDNVTMPI.PIMDPartitionedRunSupportsOneRankPerBead:FixPIMDNVTMPI.CMDLeavesCentroidChainInactive)
  add_mpi_test(NAME FixPIMDNVTMultiRankPerBead NUM_PROCS 4
               COMMAND $<TARGET_FILE:test_fix_pimd_nvt>
                       --gtest_filter=FixPIMDNVTMPI.MultiRankPerBeadRunProducesFiniteThermostatState)

  add_executable(test_fix_pimd_uvt test_fix_pimd_uvt.cpp)
  target_link_libraries(test_fix_pimd_uvt PRIVATE lammps GTest::GMock)
  add_test(NAME FixPIMDUVTSerial COMMAND test_fix_pimd_uvt --gtest_filter=FixPIMDUVTSerialTest.*)
  string(CONCAT pimd_uvt_partition_filter
         "--gtest_filter=FixPIMDUVTMPI.InitialVectorDerivativeUsesBeadCoordinates:"
         "FixPIMDUVTMPI.BeadAveragedQuadraticDerivativeAtFixedPoint:"
         "FixPIMDUVTMPI.SharedElectronicMassAndEnergyAreConsistentAcrossBeads:"
         "FixPIMDUVTMPI.ElectronicVelocityUsesMeanChainFriction")
  add_mpi_test(NAME FixPIMDUVTMPI NUM_PROCS 2
               COMMAND $<TARGET_FILE:test_fix_pimd_uvt> "${pimd_uvt_partition_filter}")
  add_mpi_test(NAME FixPIMDUVTMultiRankPerBead NUM_PROCS 4
               COMMAND $<TARGET_FILE:test_fix_pimd_uvt>
                       --gtest_filter=FixPIMDUVTMPI.MultiRankPerBeadAveragesDerivativeOncePerBead:FixPIMDUVTMPI.ElectronicVelocityUsesMeanChainFriction)
  add_mpi_test(NAME FixPIMDUVTSharedElectronicMass NUM_PROCS 4
               COMMAND $<TARGET_FILE:test_fix_pimd_uvt>
                       --gtest_filter=FixPIMDUVTMPI.SharedElectronicMassAndEnergyAreConsistentAcrossBeads)

endif()

add_executable(test_regions test_regions.cpp)
target_link_libraries(test_regions PRIVATE lammps GTest::GMock)
add_test(NAME Regions COMMAND test_regions)

add_executable(test_delete_atoms test_delete_atoms.cpp)
target_link_libraries(test_delete_atoms PRIVATE lammps GTest::GMock)
add_test(NAME DeleteAtoms COMMAND test_delete_atoms)

add_executable(test_set_property test_set_property.cpp)
target_link_libraries(test_set_property PRIVATE lammps GTest::GMock)
add_test(NAME SetProperty COMMAND test_set_property)

add_executable(test_pair_ldd test_pair_ldd.cpp)
target_link_libraries(test_pair_ldd PRIVATE lammps GTest::GMock)
add_test(NAME PairLDD COMMAND test_pair_ldd)

add_executable(test_labelmap test_labelmap.cpp)
target_compile_definitions(test_labelmap PRIVATE -DTEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_labelmap PRIVATE lammps GTest::GMock)
add_test(NAME Labelmap COMMAND test_labelmap)

add_executable(test_variables test_variables.cpp)
target_link_libraries(test_variables PRIVATE lammps GTest::GMock)
add_test(NAME Variables COMMAND test_variables)

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

if(PKG_MOLECULE)
  add_executable(test_compute_global test_compute_global.cpp)
  target_compile_definitions(test_compute_global PRIVATE -DTEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(test_compute_global PRIVATE lammps GTest::GMock)
  add_test(NAME ComputeGlobal COMMAND test_compute_global)

  add_executable(test_compute_chunk test_compute_chunk.cpp)
  target_compile_definitions(test_compute_chunk PRIVATE -DTEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(test_compute_chunk PRIVATE lammps GTest::GMock)
  add_test(NAME ComputeChunk COMMAND test_compute_chunk)
endif()

add_executable(test_mpi_load_balancing test_mpi_load_balancing.cpp)
target_link_libraries(test_mpi_load_balancing PRIVATE lammps GTest::GMock)
target_compile_definitions(test_mpi_load_balancing PRIVATE ${TEST_CONFIG_DEFS})
add_mpi_test(NAME MPILoadBalancing NUM_PROCS 4 COMMAND $<TARGET_FILE:test_mpi_load_balancing>)

if(PKG_REPLICA)
  add_executable(test_fix_pimd_langevin test_fix_pimd_langevin.cpp)
  target_link_libraries(test_fix_pimd_langevin PRIVATE lammps GTest::GMock)
  add_test(NAME FixPIMDLangevin
           COMMAND test_fix_pimd_langevin
                   --gtest_filter=FixPIMDLangevinSerialTest.*)
  string(CONCAT pimd_langevin_partition_filter
         "--gtest_filter=FixPIMDLangevinMPI.NormalModeIntegration:"
         "FixPIMDLangevinMPI.CartesianIntegration:"
         "FixPIMDLangevinRestartTest.NVTStateRestored")
  add_mpi_test(NAME FixPIMDLangevinPartition NUM_PROCS 2
               COMMAND $<TARGET_FILE:test_fix_pimd_langevin>
                       "${pimd_langevin_partition_filter}")
  string(CONCAT pimd_langevin_multirank_filter
         "--gtest_filter=FixPIMDLangevinMPI.NMPIMDNVTOBABOMultiRankPerBead:"
         "FixPIMDLangevinMPI.IsotropicBarostats:"
         "FixPIMDLangevinRestartTest.NMPIMDNPTBZPRestartContinuation:"
         "FixPIMDLangevinExpectedErrorTest.PIMDRejectsMultipleRanksPerBead")
  add_mpi_test(NAME FixPIMDLangevinMultiRank NUM_PROCS 4
               COMMAND $<TARGET_FILE:test_fix_pimd_langevin>
                       "${pimd_langevin_multirank_filter}")

  add_mpi_test(NAME FixPIMDLangevinAnisotropic NUM_PROCS 4
               COMMAND $<TARGET_FILE:test_fix_pimd_langevin>
                       --gtest_filter=FixPIMDLangevinMPI.AnisotropicBarostatStateMatchesAcrossPartitions)

  add_mpi_test(NAME FixPIMDLangevinBosonic NUM_PROCS 4
               COMMAND $<TARGET_FILE:test_fix_pimd_langevin>
                       --gtest_filter=FixPIMDLangevinBosonicTest.*)
endif()

add_executable(test_restart_settings test_restart_settings.cpp)
target_link_libraries(test_restart_settings PRIVATE lammps GTest::GMock)
add_test(NAME RestartSettings COMMAND test_restart_settings)

add_executable(test_input_validation test_input_validation.cpp)
target_link_libraries(test_input_validation PRIVATE lammps GTest::GMock)
add_test(NAME InputValidation COMMAND test_input_validation)

# Partition variables and restart files are local to each test invocation.
get_property(command_tests DIRECTORY PROPERTY TESTS)
foreach(test_name IN LISTS command_tests)
  if(test_name MATCHES "^(FixUVT|ComputeTempUVT|FixPIMD)")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
    set_tests_properties(${test_name} PROPERTIES
                         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
  endif()
endforeach()
