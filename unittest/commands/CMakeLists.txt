
# build LAMMPS plugins, but not on Windows
if(NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
  ExternalProject_Add(plugins
                      SOURCE_DIR      "${LAMMPS_DIR}/examples/plugins"
                      BINARY_DIR      "${CMAKE_BINARY_DIR}/plugins"
                      CMAKE_ARGS      ${CMAKE_REQUEST_PIC}
                                      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                      -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                                      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                                      -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
                                      -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
                     BUILD_BYPRODUCTS <BINARY_DIR>/morse2plugin${CMAKE_STATIC_LIBRARY_SUFFIX}
                                      <BINARY_DIR>/nve2plugin${CMAKE_STATIC_LIBRARY_SUFFIX}
                                      <BINARY_DIR>/helloplugin${CMAKE_STATIC_LIBRARY_SUFFIX}
                     INSTALL_COMMAND ""
                     TEST_COMMAND    "")
endif()

add_executable(test_simple_commands test_simple_commands.cpp)
add_dependencies(test_simple_commands plugins)
target_link_libraries(test_simple_commands PRIVATE lammps GTest::GMock GTest::GTest)
add_test(NAME SimpleCommands COMMAND test_simple_commands WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_lattice_region test_lattice_region.cpp)
target_link_libraries(test_lattice_region PRIVATE lammps GTest::GMock GTest::GTest)
add_test(NAME LatticeRegion COMMAND test_lattice_region WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_kim_commands test_kim_commands.cpp)
if(KIM_EXTRA_UNITTESTS)
  if(CURL_FOUND)
    target_compile_definitions(test_kim_commands PRIVATE -DKIM_EXTRA_UNITTESTS)
  else()
    message(FATAL_ERROR "CURL not found. Enabling KIM extra unit tests requires to have libcurl installed.")
  endif()
endif()
target_link_libraries(test_kim_commands PRIVATE lammps GTest::GMock GTest::GTest)
add_test(NAME KimCommands COMMAND test_kim_commands WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_reset_ids test_reset_ids.cpp)
target_compile_definitions(test_reset_ids PRIVATE -DTEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_reset_ids PRIVATE lammps GTest::GMock GTest::GTest)
add_test(NAME ResetIDs COMMAND test_reset_ids WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
