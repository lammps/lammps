
add_executable(test_simple_commands test_simple_commands.cpp)

# tests for the plugin command require the PLUGIN package and won't work on windows
if((NOT (CMAKE_SYSTEM_NAME STREQUAL "Windows")) AND PKG_PLUGIN)
  add_subdirectory(${LAMMPS_DIR}/examples/plugins ${CMAKE_BINARY_DIR}/build-plugins)
  add_dependencies(test_simple_commands plugins)
  target_compile_definitions(test_simple_commands PRIVATE -DLMP_PLUGIN)
endif()

target_link_libraries(test_simple_commands PRIVATE lammps gmock)
add_test(NAME SimpleCommands COMMAND test_simple_commands WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
if(BUILD_IS_MULTI_CONFIG)
  set(LAMMPS_PLUGIN_BIN_DIR ${CMAKE_BINARY_DIR}/build-plugins/$<CONFIG>)
else()
  set(LAMMPS_PLUGIN_BIN_DIR ${CMAKE_BINARY_DIR}/build-plugins)
endif()
if(APPLE)
  set_tests_properties(SimpleCommands PROPERTIES
    ENVIRONMENT "DYLD_LIBRARY_PATH=${LAMMPS_PLUGIN_BIN_DIR}/:${DYLD_LIBRARY_PATH}")
else()
  set_tests_properties(SimpleCommands PROPERTIES
    ENVIRONMENT "LD_LIBRARY_PATH=${LAMMPS_PLUGIN_BIN_DIR}/:${LD_LIBRARY_PATH}")
endif()

add_executable(test_lattice_region test_lattice_region.cpp)
target_link_libraries(test_lattice_region PRIVATE lammps gmock)
add_test(NAME LatticeRegion COMMAND test_lattice_region WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_groups test_groups.cpp)
target_link_libraries(test_groups PRIVATE lammps gmock)
add_test(NAME Groups COMMAND test_groups WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_variables test_variables.cpp)
target_link_libraries(test_variables PRIVATE lammps gmock)
add_test(NAME Variables COMMAND test_variables WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_kim_commands test_kim_commands.cpp)
if(KIM_EXTRA_UNITTESTS)
  if(CURL_FOUND)
    target_compile_definitions(test_kim_commands PRIVATE -DKIM_EXTRA_UNITTESTS)
  else()
    message(FATAL_ERROR "CURL not found. Enabling KIM extra unit tests requires to have libcurl installed.")
  endif()
endif()
target_link_libraries(test_kim_commands PRIVATE lammps gmock)
add_test(NAME KimCommands COMMAND test_kim_commands WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_reset_ids test_reset_ids.cpp)
target_compile_definitions(test_reset_ids PRIVATE -DTEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_reset_ids PRIVATE lammps gmock)
add_test(NAME ResetIDs COMMAND test_reset_ids WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

if(BUILD_MPI)
  add_executable(test_mpi_load_balancing test_mpi_load_balancing.cpp)
  target_link_libraries(test_mpi_load_balancing PRIVATE lammps gmock)
  target_compile_definitions(test_mpi_load_balancing PRIVATE ${TEST_CONFIG_DEFS})
  add_mpi_test(NAME MPILoadBalancing NUM_PROCS 4 COMMAND $<TARGET_FILE:test_mpi_load_balancing>)
endif()
