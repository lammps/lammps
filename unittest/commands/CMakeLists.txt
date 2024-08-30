# -*- CMake -*- file for testing LAMMPS commands

add_executable(test_simple_commands test_simple_commands.cpp)

# tests for the plugin command require the PLUGIN package and won't work on windows
if((NOT (CMAKE_SYSTEM_NAME STREQUAL "Windows")) AND PKG_PLUGIN)
  add_subdirectory(${LAMMPS_DIR}/examples/plugins ${CMAKE_BINARY_DIR}/build-plugins)
  add_dependencies(test_simple_commands plugins)
  target_compile_definitions(test_simple_commands PRIVATE -DLMP_PLUGIN)
endif()

target_link_libraries(test_simple_commands PRIVATE lammps GTest::GMock)
add_test(NAME SimpleCommands COMMAND test_simple_commands)
set_tests_properties(SimpleCommands PROPERTIES
          ENVIRONMENT "LAMMPS_PLUGIN_BIN_DIR=${CMAKE_BINARY_DIR}")

add_executable(test_lattice_region test_lattice_region.cpp)
target_link_libraries(test_lattice_region PRIVATE lammps GTest::GMock)
add_test(NAME LatticeRegion COMMAND test_lattice_region)

add_executable(test_groups test_groups.cpp)
target_link_libraries(test_groups PRIVATE lammps GTest::GMock)
add_test(NAME Groups COMMAND test_groups)

add_executable(test_regions test_regions.cpp)
target_link_libraries(test_regions PRIVATE lammps GTest::GMock)
add_test(NAME Regions COMMAND test_regions)

add_executable(test_delete_atoms test_delete_atoms.cpp)
target_link_libraries(test_delete_atoms PRIVATE lammps GTest::GMock)
add_test(NAME DeleteAtoms COMMAND test_delete_atoms)

add_executable(test_set_property test_set_property.cpp)
target_link_libraries(test_set_property PRIVATE lammps GTest::GMock)
add_test(NAME SetProperty COMMAND test_set_property)

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
