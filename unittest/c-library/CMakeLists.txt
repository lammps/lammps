# -*- CMake -*- file for testing the c-library interface

add_executable(test_library_open test_library_open.cpp test_main.cpp)
target_link_libraries(test_library_open PRIVATE lammps GTest::GMock)
add_test(NAME LibraryOpen COMMAND test_library_open)

add_executable(test_library_commands test_library_commands.cpp test_main.cpp)
target_link_libraries(test_library_commands PRIVATE lammps GTest::GMock)
add_test(NAME LibraryCommands COMMAND test_library_commands)

add_executable(test_library_external test_library_external.cpp test_main.cpp)
target_link_libraries(test_library_external PRIVATE lammps GTest::GMock)
add_test(NAME LibraryExternal COMMAND test_library_external)

add_executable(test_library_properties test_library_properties.cpp test_main.cpp)
target_link_libraries(test_library_properties PRIVATE lammps GTest::GMock)
target_compile_definitions(test_library_properties PRIVATE -DTEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME LibraryProperties COMMAND test_library_properties)
set_tests_properties(LibraryProperties PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR}")

add_executable(test_library_scatter_gather test_library_scatter_gather.cpp test_main.cpp)
target_link_libraries(test_library_scatter_gather PRIVATE lammps GTest::GMock)
target_compile_definitions(test_library_scatter_gather PRIVATE -DTEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME LibraryScatterGather COMMAND test_library_scatter_gather)
set_tests_properties(LibraryScatterGather PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR}")

set(TEST_CONFIG_DEFS "-DTEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR};-DLAMMPS_${LAMMPS_SIZES}")
set(PKG_COUNT 0)
foreach(PKG ${STANDARD_PACKAGES} ${EXTRA_PACKAGES} ${SUFFIX_PACKAGES})
  if(PKG_${PKG})
    MATH(EXPR PKG_COUNT "${PKG_COUNT}+1")
  endif()
endforeach()
list(APPEND TEST_CONFIG_DEFS -DNUM_LAMMPS_PACKAGES=${PKG_COUNT})

if(PKG_MANYBODY)
  set(HAS_MANYBODY 1)
else()
  set(HAS_MANYBODY 0)
endif()
list(APPEND TEST_CONFIG_DEFS -DLAMMPS_HAS_MANYBODY=${HAS_MANYBODY})

if(BUILD_MPI)
  set(HAS_MPI 1)
else()
  set(HAS_MPI 0)
endif()
list(APPEND TEST_CONFIG_DEFS -DLAMMPS_HAS_MPI=${HAS_MPI})

foreach(WITH "JPEG" "PNG" "GZIP" "FFMPEG")
  if(WITH_${WITH})
    set(HAS_${WITH} 1)
  else()
    set(HAS_${WITH} 0)
  endif()
  list(APPEND TEST_CONFIG_DEFS -DLAMMPS_HAS_${WITH}=${HAS_${WITH}})
endforeach()

add_executable(test_library_config test_library_config.cpp test_main.cpp)
target_link_libraries(test_library_config PRIVATE lammps GTest::GMock)
target_compile_definitions(test_library_config PRIVATE ${TEST_CONFIG_DEFS})
add_test(NAME LibraryConfig COMMAND test_library_config)

add_executable(test_library_mpi test_library_mpi.cpp)
target_link_libraries(test_library_mpi PRIVATE lammps GTest::GMock)
target_compile_definitions(test_library_mpi PRIVATE ${TEST_CONFIG_DEFS})
add_mpi_test(NAME LibraryMPI NUM_PROCS 4 COMMAND $<TARGET_FILE:test_library_mpi>)

# simple run tests for coupling to the LAMMPS library
if(BUILD_MPI)
  if(BUILD_SHARED_LIBS)
    add_subdirectory(${LAMMPS_DIR}/examples/COUPLE/plugin ${CMAKE_BINARY_DIR}/build-couple)
    add_test(NAME RunCoupleSimplePlugin
      COMMAND $<TARGET_FILE:simple-plugin> 1 ${LAMMPS_DIR}/examples/COUPLE/plugin/in.lj $<TARGET_FILE:lammps>)
    set_tests_properties(RunCoupleSimplePlugin PROPERTIES
      ENVIRONMENT "TSAN_OPTIONS=ignore_noninstrumented_modules=1;HWLOC_HIDE_ERRORS=2"
      PASS_REGULAR_EXPRESSION "LAMMPS \\([0-9]+ [A-Za-z]+ 2[0-9][0-9][0-9]( - Update [0-9]+)?( - Development.*)?( - Maintenance.*)?\\)")
  endif()
  add_subdirectory(${LAMMPS_DIR}/examples/COUPLE/simple ${CMAKE_BINARY_DIR}/build-simple)
  add_test(NAME RunCoupleSimpleC
    COMMAND $<TARGET_FILE:simpleC> 1 ${LAMMPS_DIR}/examples/COUPLE/simple/in.lj)
  add_test(NAME RunCoupleSimpleCC
    COMMAND $<TARGET_FILE:simpleCC> 1 ${LAMMPS_DIR}/examples/COUPLE/simple/in.lj)

  if(TARGET simpleF90)
    add_test(NAME RunCoupleSimpleF90
      COMMAND $<TARGET_FILE:simpleF90> 1 ${LAMMPS_DIR}/examples/COUPLE/simple/in.lj)
    set_tests_properties(RunCoupleSimpleF90 PROPERTIES
      ENVIRONMENT "TSAN_OPTIONS=ignore_noninstrumented_modules=1;HWLOC_HIDE_ERRORS=2"
      PASS_REGULAR_EXPRESSION "LAMMPS \\([0-9]+ [A-Za-z]+ 2[0-9][0-9][0-9]( - Update [0-9]+)?( - Development.*)?( - Maintenance.*)?\\)")
  endif()

  set_tests_properties(RunCoupleSimpleC RunCoupleSimpleCC PROPERTIES
    ENVIRONMENT "TSAN_OPTIONS=ignore_noninstrumented_modules=1;HWLOC_HIDE_ERRORS=2"
    PASS_REGULAR_EXPRESSION "LAMMPS \\([0-9]+ [A-Za-z]+ 2[0-9][0-9][0-9]( - Update [0-9]+)?( - Development.*)?( - Maintenance.*)?\\)")
endif()
