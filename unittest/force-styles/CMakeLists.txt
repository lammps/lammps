
find_package(YAML)
if(NOT YAML_FOUND)
  # download and build a local copy of libyaml
  include(YAML)
endif()

if(CMAKE_VERSION VERSION_LESS 3.12)
  # adjust so we find Python 3 versions before Python 2 on old systems with old CMake
  set(Python_ADDITIONAL_VERSIONS 3.8 3.7 3.6 3.5)
  find_package(PythonInterp) # Deprecated since version 3.12
  if(PYTHONINTERP_FOUND)
    set(Python_EXECUTABLE ${PYTHON_EXECUTABLE})
  endif()
else()
  find_package(Python COMPONENTS Interpreter)
endif()
if(Python_EXECUTABLE)
  add_custom_target(check-tests
    ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_tests.py
    -s ${LAMMPS_SOURCE_DIR} -t ${CMAKE_CURRENT_SOURCE_DIR}/tests
    COMMENT "Check completeness of force style tests")
endif()

set(TEST_INPUT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR}/tests)
add_library(style_tests STATIC yaml_writer.cpp error_stats.cpp test_config_reader.cpp test_main.cpp)
target_compile_definitions(style_tests PRIVATE -DTEST_INPUT_FOLDER=${TEST_INPUT_FOLDER})
target_include_directories(style_tests PRIVATE ${LAMMPS_SOURCE_DIR})
target_link_libraries(style_tests PUBLIC GTest::GTest GTest::GMock Yaml::Yaml lammps)
if(BUILD_MPI)
  target_link_libraries(style_tests PUBLIC MPI::MPI_CXX)
else()
  target_link_libraries(style_tests PUBLIC mpi_stubs)
endif()
# propagate sanitizer options to test tools
if(ENABLE_SANITIZER AND (NOT (ENABLE_SANITIZER STREQUAL "none")))
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_compile_options(style_tests PUBLIC -fsanitize=${ENABLE_SANITIZER})
    target_link_options(style_tests PUBLIC -fsanitize=${ENABLE_SANITIZER})
  endif()
endif()

# unit test for error stats class
add_executable(test_error_stats test_error_stats.cpp)
target_include_directories(test_error_stats PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${LAMMPS_SOURCE_DIR})
target_link_libraries(test_error_stats PRIVATE GTest::GTestMain GTest::GTest)
add_test(NAME ErrorStats COMMAND test_error_stats)

# pair style tester
add_executable(test_pair_style test_pair_style.cpp)
target_link_libraries(test_pair_style PRIVATE lammps style_tests)
if(FFT_SINGLE)
  target_compile_definitions(test_pair_style PRIVATE -DFFT_SINGLE)
endif()

# tests for molecular systems and related pair styles
file(GLOB MOL_PAIR_TESTS LIST_DIRECTORIES false ${TEST_INPUT_FOLDER}/mol-pair-*.yaml)
# cannot test MSM with single precision data
if(FFT_SINGLE)
  list(FILTER MOL_PAIR_TESTS EXCLUDE REGEX "msm")
endif()
foreach(TEST ${MOL_PAIR_TESTS})
  string(REGEX REPLACE "^.*mol-pair-(.*)\.yaml" "MolPairStyle:\\1" TNAME ${TEST})
  add_test(NAME ${TNAME} COMMAND test_pair_style ${TEST} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(${TNAME} PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR};PYTHONPATH=${TEST_INPUT_FOLDER}:$ENV{PYTHONPATH}")
endforeach()

# tests for metal-like atomic systems and related pair styles
file(GLOB ATOMIC_PAIR_TESTS LIST_DIRECTORIES false ${TEST_INPUT_FOLDER}/atomic-pair-*.yaml)
foreach(TEST ${ATOMIC_PAIR_TESTS})
  string(REGEX REPLACE "^.*atomic-pair-(.*)\.yaml" "AtomicPairStyle:\\1" TNAME ${TEST})
  add_test(NAME ${TNAME} COMMAND test_pair_style ${TEST} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(${TNAME} PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR};PYTHONPATH=${TEST_INPUT_FOLDER}:$ENV{PYTHONPATH}")
endforeach()

# tests for Si-like manybody systems and related pair styles
file(GLOB MANYBODY_PAIR_TESTS LIST_DIRECTORIES false ${TEST_INPUT_FOLDER}/manybody-pair-*.yaml)
foreach(TEST ${MANYBODY_PAIR_TESTS})
  string(REGEX REPLACE "^.*manybody-pair-(.*)\.yaml" "ManybodyPairStyle:\\1" TNAME ${TEST})
  add_test(NAME ${TNAME} COMMAND test_pair_style ${TEST} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(${TNAME} PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR};PYTHONPATH=${TEST_INPUT_FOLDER}:$ENV{PYTHONPATH}")
endforeach()

# bond style tester
add_executable(test_bond_style test_bond_style.cpp)
target_link_libraries(test_bond_style PRIVATE lammps style_tests)

file(GLOB BOND_TESTS LIST_DIRECTORIES false ${TEST_INPUT_FOLDER}/bond-*.yaml)
foreach(TEST ${BOND_TESTS})
  string(REGEX REPLACE "^.*bond-(.*)\.yaml" "BondStyle:\\1" TNAME ${TEST})
  add_test(NAME ${TNAME} COMMAND test_bond_style ${TEST} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(${TNAME} PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR};PYTHONPATH=${TEST_INPUT_FOLDER}:$ENV{PYTHONPATH}")
endforeach()

# angle style tester
add_executable(test_angle_style test_angle_style.cpp)
target_link_libraries(test_angle_style PRIVATE lammps style_tests)

file(GLOB ANGLE_TESTS LIST_DIRECTORIES false ${TEST_INPUT_FOLDER}/angle-*.yaml)
foreach(TEST ${ANGLE_TESTS})
  string(REGEX REPLACE "^.*angle-(.*)\.yaml" "AngleStyle:\\1" TNAME ${TEST})
  add_test(NAME ${TNAME} COMMAND test_angle_style ${TEST} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(${TNAME} PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR};PYTHONPATH=${TEST_INPUT_FOLDER}:$ENV{PYTHONPATH}")
endforeach()

# kspace style tester, currently uses the pair style tool
file(GLOB KSPACE_TESTS LIST_DIRECTORIES false ${TEST_INPUT_FOLDER}/kspace-*.yaml)
# cannot test MSM with single precision data
if(FFT_SINGLE)
  list(FILTER KSPACE_TESTS EXCLUDE REGEX "msm")
endif()
foreach(TEST ${KSPACE_TESTS})
  string(REGEX REPLACE "^.*kspace-(.*)\.yaml" "KSpaceStyle:\\1" TNAME ${TEST})
  add_test(NAME ${TNAME} COMMAND test_pair_style ${TEST} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(${TNAME} PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR};PYTHONPATH=${TEST_INPUT_FOLDER}:$ENV{PYTHONPATH}")
endforeach()

# tester for timestepping fixes
add_executable(test_fix_timestep test_fix_timestep.cpp)
if(NOT BUILD_SHARED_LIBS)
  target_compile_definitions(test_fix_timestep PRIVATE USING_STATIC_LIBS=1)
endif()
target_link_libraries(test_fix_timestep PRIVATE lammps style_tests)

# tests for timestep related fixes (time integration, thermostat, force manipulation, constraints/restraints)
file(GLOB FIX_TIMESTEP_TESTS LIST_DIRECTORIES false ${TEST_INPUT_FOLDER}/fix-timestep-*.yaml)
foreach(TEST ${FIX_TIMESTEP_TESTS})
  string(REGEX REPLACE "^.*fix-timestep-(.*)\.yaml" "FixTimestep:\\1" TNAME ${TEST})
  add_test(NAME ${TNAME} COMMAND test_fix_timestep ${TEST} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(${TNAME} PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR};PYTHONPATH=${TEST_INPUT_FOLDER}:${LAMMPS_PYTHON_DIR}:$ENV{PYTHONPATH}")
endforeach()

# dihedral style tester
add_executable(test_dihedral_style test_dihedral_style.cpp)
target_link_libraries(test_dihedral_style PRIVATE lammps style_tests)

file(GLOB DIHEDRAL_TESTS LIST_DIRECTORIES false ${TEST_INPUT_FOLDER}/dihedral-*.yaml)
foreach(TEST ${DIHEDRAL_TESTS})
  string(REGEX REPLACE "^.*dihedral-(.*)\.yaml" "DihedralStyle:\\1" TNAME ${TEST})
  add_test(NAME ${TNAME} COMMAND test_dihedral_style ${TEST} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(${TNAME} PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR};PYTHONPATH=${TEST_INPUT_FOLDER}:$ENV{PYTHONPATH}")
endforeach()

# improper style tester
add_executable(test_improper_style test_improper_style.cpp)
target_link_libraries(test_improper_style PRIVATE lammps style_tests)

file(GLOB IMPROPER_TESTS LIST_DIRECTORIES false ${TEST_INPUT_FOLDER}/improper-*.yaml)
foreach(TEST ${IMPROPER_TESTS})
  string(REGEX REPLACE "^.*improper-(.*)\.yaml" "ImproperStyle:\\1" TNAME ${TEST})
  add_test(NAME ${TNAME} COMMAND test_improper_style ${TEST} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(${TNAME} PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR};PYTHONPATH=${TEST_INPUT_FOLDER}:$ENV{PYTHONPATH}")
endforeach()
