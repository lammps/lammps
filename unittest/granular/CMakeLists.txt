# -*- CMake -*- file for tests using DEM / granular models

function(extract_tags out yaml_file)
  file(STRINGS ${yaml_file} TAGS_LINE REGEX "^tags:")
  string(REPLACE "tags:" "" TAGS_LINE "${TAGS_LINE}")
  string(REGEX MATCHALL "[^, \t\r\n]+" TAGS "${TAGS_LINE}")
  set(${out} "${TAGS}" PARENT_SCOPE)
endfunction()

# build support library
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  string(REPLACE "/" "\\\\" TEST_INPUT_FOLDER "${CMAKE_CURRENT_SOURCE_DIR}/tests")
else()
  set(TEST_INPUT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR}/tests)
endif()
add_library(granular_tests STATIC
  ../error_stats.cpp
  ../yaml_writer.cpp
  test_config_reader.cpp
  test_analytic_models.cpp
  test_dem_common.cpp
  test_main.cpp)
if(YAML_FOUND)
  target_compile_definitions(granular_tests PRIVATE TEST_INPUT_FOLDER=${TEST_INPUT_FOLDER})
else()
  # we always use static linkage with local compiled libyaml
  target_compile_definitions(granular_tests PRIVATE TEST_INPUT_FOLDER=${TEST_INPUT_FOLDER} YAML_DECLARE_STATIC)
endif()
target_include_directories(granular_tests PRIVATE "${LAMMPS_SOURCE_DIR};${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(granular_tests PUBLIC GTest::GMock Yaml::Yaml lammps)

# propagate sanitizer options to test tools
if(ENABLE_SANITIZER AND (NOT (ENABLE_SANITIZER STREQUAL "none")))
  target_compile_options(granular_tests PUBLIC -fsanitize=${ENABLE_SANITIZER})
  target_link_options(granular_tests PUBLIC -fsanitize=${ENABLE_SANITIZER})
endif()

# prepare environment for testers
if(WIN32)
  set(GRANULAR_TEST_ENVIRONMENT PYTHONPATH=${TEST_INPUT_FOLDER})
else()
  set(GRANULAR_TEST_ENVIRONMENT PYTHONPATH=${TEST_INPUT_FOLDER}:${LAMMPS_PYTHON_DIR}:$ENV{PYTHONPATH})
endif()
list(APPEND GRANULAR_TEST_ENVIRONMENT "PYTHONUNBUFFERED=1")
list(APPEND GRANULAR_TEST_ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1")
list(APPEND GRANULAR_TEST_ENVIRONMENT "OMP_PROC_BIND=false")
list(APPEND GRANULAR_TEST_ENVIRONMENT "OMP_NUM_THREADS=4")
list(APPEND GRANULAR_TEST_ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR}")
get_property(BUILD_IS_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(BUILD_IS_MULTI_CONFIG)
  set(LAMMPS_LIB_PATH ${CMAKE_BINARY_DIR}/$<CONFIG>)
else()
  set(LAMMPS_LIB_PATH ${CMAKE_BINARY_DIR})
endif()
if(APPLE)
  list(APPEND GRANULAR_TEST_ENVIRONMENT "DYLD_LIBRARY_PATH=${LAMMPS_LIB_PATH}:$ENV{DYLD_LIBRARY_PATH}")
elseif(WIN32)
  list(APPEND GRANULAR_TEST_ENVIRONMENT "LAMMPSDLLPATH=${LAMMPS_LIB_PATH}")
else()
  list(APPEND GRANULAR_TEST_ENVIRONMENT "LD_LIBRARY_PATH=${LAMMPS_LIB_PATH}:$ENV{LD_LIBRARY_PATH}")
endif()

# helper to register all YAML files matching a prefix with a given test driver
function(register_dem_tests driver prefix label)
  file(GLOB DEM_TESTS LIST_DIRECTORIES false CONFIGURE_DEPENDS ${TEST_INPUT_FOLDER}/${prefix}-*.yaml)
  foreach(TEST ${DEM_TESTS})
    string(REGEX REPLACE "^.*${prefix}-(.*)\.yaml" "${label}:\\1" TNAME ${TEST})
    extract_tags(TEST_TAGS ${TEST})
    list(FIND TEST_TAGS "no${CMAKE_SYSTEM_NAME}" _SKIPME)
    if(_SKIPME GREATER_EQUAL 0)
      continue()
    endif()
    add_test(NAME ${TNAME} COMMAND ${driver} ${TEST})
    set_tests_properties(${TNAME} PROPERTIES ENVIRONMENT "${GRANULAR_TEST_ENVIRONMENT}")
    set_tests_properties(${TNAME} PROPERTIES LABELS "${TEST_TAGS}")
  endforeach()
endfunction()

# DEM-01: freely falling particle bouncing off a wall
add_executable(test_dem_01 test_dem_01.cpp)
target_link_libraries(test_dem_01 PRIVATE lammps granular_tests)
register_dem_tests(test_dem_01 dem01 DEM01)

# DEM-02: repeated bouncing particle (convergence to the hard-sphere limit)
add_executable(test_dem_02 test_dem_02.cpp)
target_link_libraries(test_dem_02 PRIVATE lammps granular_tests)
register_dem_tests(test_dem_02 dem02 DEM02)

# DEM-03: two stacked particles in continuous compression between two walls
add_executable(test_dem_03 test_dem_03.cpp)
target_link_libraries(test_dem_03 PRIVATE lammps granular_tests)
register_dem_tests(test_dem_03 dem03 DEM03)

# DEM-04: sphere sliding then rolling without slipping on a rough surface
add_executable(test_dem_04 test_dem_04.cpp)
target_link_libraries(test_dem_04 PRIVATE lammps granular_tests)
register_dem_tests(test_dem_04 dem04 DEM04)

# DEM-05: oblique collision of a sphere with a flat wall
add_executable(test_dem_05 test_dem_05.cpp)
target_link_libraries(test_dem_05 PRIVATE lammps granular_tests)
register_dem_tests(test_dem_05 dem05 DEM05)

# DEM-06: single particle settling to terminal velocity under fluid drag
add_executable(test_dem_06 test_dem_06.cpp)
target_link_libraries(test_dem_06 PRIVATE lammps granular_tests)
register_dem_tests(test_dem_06 dem06 DEM06)

# DEM-07: spinning sphere damped by rolling resistance (rolling sds)
add_executable(test_dem_07 test_dem_07.cpp)
target_link_libraries(test_dem_07 PRIVATE lammps granular_tests)
register_dem_tests(test_dem_07 dem07 DEM07)

# DEM-08: cohesive/adhesive contact (DMT/JKR pull-off force)
add_executable(test_dem_08 test_dem_08.cpp)
target_link_libraries(test_dem_08 PRIVATE lammps granular_tests)
register_dem_tests(test_dem_08 dem08 DEM08)

# DEM-09: particle-particle head-on and oblique (shear) collisions
add_executable(test_dem_09 test_dem_09.cpp)
target_link_libraries(test_dem_09 PRIVATE lammps granular_tests)
register_dem_tests(test_dem_09 dem09 DEM09)

# DEM-10: bulk angle of repose and multi-sphere rigid clumps
add_executable(test_dem_10 test_dem_10.cpp)
target_link_libraries(test_dem_10 PRIVATE lammps granular_tests)
register_dem_tests(test_dem_10 dem10 DEM10)

# DEM-11: elastic Hertzian normal impact (peak contact mechanics)
add_executable(test_dem_11 test_dem_11.cpp)
target_link_libraries(test_dem_11 PRIVATE lammps granular_tests)
register_dem_tests(test_dem_11 dem11 DEM11)
