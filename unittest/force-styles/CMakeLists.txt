
find_package(YAML)
if(NOT YAML_FOUND)
  message(STATUS "Skipping tests because libyaml is not found")
  return()
endif()

set(TEST_COMMON_FILES yaml_writer.cpp;error_stats.cpp;test_config_reader.cpp;test_main.cpp)
set(TEST_INPUT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR}/tests)

# pair style tester
add_executable(pair_style pair_style.cpp ${TEST_COMMON_FILES})
target_link_libraries(pair_style PRIVATE lammps GTest::GMock GTest::GTest Yaml::Yaml)
target_compile_definitions(pair_style PRIVATE TEST_INPUT_FOLDER=${TEST_INPUT_FOLDER})

# tests for a molecular systems and related pair styles
file(GLOB MOL_PAIR_TESTS LIST_DIRECTORIES false ${TEST_INPUT_FOLDER}/mol-pair-*.yaml)
foreach(TEST ${MOL_PAIR_TESTS})
  string(REGEX REPLACE "^.*mol-pair-(.*)\.yaml" "MolPairStyle:\\1" TNAME ${TEST})
  add_test(NAME ${TNAME} COMMAND pair_style ${TEST}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endforeach()

# tests for atomic systems and related pair styles
file(GLOB ATOMIC_PAIR_TESTS LIST_DIRECTORIES false ${TEST_INPUT_FOLDER}/atomic-pair-*.yaml)
foreach(TEST ${ATOMIC_PAIR_TESTS})
  string(REGEX REPLACE "^.*atomic-pair-(.*)\.yaml" "AtomicPairStyle:\\1" TNAME ${TEST})
  add_test(NAME ${TNAME} COMMAND pair_style ${TEST}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(${TNAME} PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR}")
endforeach()

# bond style tester
add_executable(bond_style bond_style.cpp ${TEST_COMMON_FILES})
target_link_libraries(bond_style PRIVATE lammps GTest::GMock GTest::GTest Yaml::Yaml)
target_compile_definitions(bond_style PRIVATE TEST_INPUT_FOLDER=${TEST_INPUT_FOLDER})

file(GLOB BOND_TESTS LIST_DIRECTORIES false ${TEST_INPUT_FOLDER}/bond-*.yaml)
foreach(TEST ${BOND_TESTS})
  string(REGEX REPLACE "^.*bond-(.*)\.yaml" "BondStyle:\\1" TNAME ${TEST})
  add_test(NAME ${TNAME} COMMAND bond_style ${TEST}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endforeach()

# angle style tester
add_executable(angle_style angle_style.cpp ${TEST_COMMON_FILES})
target_link_libraries(angle_style PRIVATE lammps GTest::GMock GTest::GTest Yaml::Yaml)
target_compile_definitions(angle_style PRIVATE TEST_INPUT_FOLDER=${TEST_INPUT_FOLDER})

file(GLOB ANGLE_TESTS LIST_DIRECTORIES false ${TEST_INPUT_FOLDER}/angle-*.yaml)
foreach(TEST ${ANGLE_TESTS})
  string(REGEX REPLACE "^.*angle-(.*)\.yaml" "AngleStyle:\\1" TNAME ${TEST})
  add_test(NAME ${TNAME} COMMAND angle_style ${TEST}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endforeach()

