# -*- CMake -*- file for tests of the BPM (bonded particle model) pair+bond styles

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
add_library(bpm_tests STATIC
  ../error_stats.cpp
  ../yaml_writer.cpp
  test_config_reader.cpp
  test_analytic_models.cpp
  test_bpm_common.cpp
  test_main.cpp)
if(YAML_FOUND)
  target_compile_definitions(bpm_tests PRIVATE TEST_INPUT_FOLDER=${TEST_INPUT_FOLDER})
else()
  # we always use static linkage with local compiled libyaml
  target_compile_definitions(bpm_tests PRIVATE TEST_INPUT_FOLDER=${TEST_INPUT_FOLDER} YAML_DECLARE_STATIC)
endif()
target_include_directories(bpm_tests PRIVATE ${LAMMPS_SOURCE_DIR})
target_link_libraries(bpm_tests PUBLIC GTest::GMock Yaml::Yaml lammps)

# propagate sanitizer options to test tools
if(ENABLE_SANITIZER AND (NOT (ENABLE_SANITIZER STREQUAL "none")))
  target_compile_options(bpm_tests PUBLIC -fsanitize=${ENABLE_SANITIZER})
  target_link_options(bpm_tests PUBLIC -fsanitize=${ENABLE_SANITIZER})
endif()

# prepare environment for testers
if(WIN32)
  set(BPM_TEST_ENVIRONMENT PYTHONPATH=${TEST_INPUT_FOLDER})
else()
  set(BPM_TEST_ENVIRONMENT PYTHONPATH=${TEST_INPUT_FOLDER}:${LAMMPS_PYTHON_DIR}:$ENV{PYTHONPATH})
endif()
list(APPEND BPM_TEST_ENVIRONMENT "PYTHONUNBUFFERED=1")
list(APPEND BPM_TEST_ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1")
list(APPEND BPM_TEST_ENVIRONMENT "OMP_PROC_BIND=false")
list(APPEND BPM_TEST_ENVIRONMENT "OMP_NUM_THREADS=4")
list(APPEND BPM_TEST_ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR}")
get_property(BUILD_IS_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(BUILD_IS_MULTI_CONFIG)
  set(LAMMPS_LIB_PATH ${CMAKE_BINARY_DIR}/$<CONFIG>)
else()
  set(LAMMPS_LIB_PATH ${CMAKE_BINARY_DIR})
endif()
if(APPLE)
  list(APPEND BPM_TEST_ENVIRONMENT "DYLD_LIBRARY_PATH=${LAMMPS_LIB_PATH}:$ENV{DYLD_LIBRARY_PATH}")
elseif(WIN32)
  list(APPEND BPM_TEST_ENVIRONMENT "LAMMPSDLLPATH=${LAMMPS_LIB_PATH}")
else()
  list(APPEND BPM_TEST_ENVIRONMENT "LD_LIBRARY_PATH=${LAMMPS_LIB_PATH}:$ENV{LD_LIBRARY_PATH}")
endif()

# helper to register all YAML files matching a prefix with a given test driver
function(register_bpm_tests driver prefix label)
  file(GLOB BPM_TESTS LIST_DIRECTORIES false CONFIGURE_DEPENDS ${TEST_INPUT_FOLDER}/${prefix}-*.yaml)
  foreach(TEST ${BPM_TESTS})
    string(REGEX REPLACE "^.*${prefix}-(.*)\.yaml" "${label}:\\1" TNAME ${TEST})
    extract_tags(TEST_TAGS ${TEST})
    list(FIND TEST_TAGS "no${CMAKE_SYSTEM_NAME}" _SKIPME)
    if(_SKIPME GREATER_EQUAL 0)
      continue()
    endif()
    add_test(NAME ${TNAME} COMMAND ${driver} ${TEST})
    set_tests_properties(${TNAME} PROPERTIES ENVIRONMENT "${BPM_TEST_ENVIRONMENT}")
    set_tests_properties(${TNAME} PROPERTIES LABELS "${TEST_TAGS}")
  endforeach()
endfunction()

# BPM peridynamics: pair_style bpm/peri + bond_style bpm/peri (pmb/lps/ves/eps)
add_executable(test_bpm_peri test_bpm_peri.cpp)
target_link_libraries(test_bpm_peri PRIVATE lammps bpm_tests)
register_bpm_tests(test_bpm_peri bpmperi BpmPeri)

# BPM zero bond style: force-free topology/breaking style (same shared driver)
add_executable(test_bpm_zero test_bpm_zero.cpp)
target_link_libraries(test_bpm_zero PRIVATE lammps bpm_tests)
register_bpm_tests(test_bpm_zero bpmzero BpmZero)

# BPM central-force bond styles: bpm/spring and bpm/spring/plastic
add_executable(test_bpm_spring test_bpm_spring.cpp)
target_link_libraries(test_bpm_spring PRIVATE lammps bpm_tests)
register_bpm_tests(test_bpm_spring bpmspring BpmSpring)
