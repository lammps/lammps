
find_package(YAML)
if(NOT YAML_FOUND)
  message(STATUS "Skipping tests because libyaml is not found")
  return()
endif()

set(TEST_COMMON_FILES yaml_writer.cpp;error_stats.cpp;test_config_reader.cpp;test_main.cpp)

# pair style tester
add_executable(pair_style pair_style.cpp ${TEST_COMMON_FILES})
target_link_libraries(pair_style PRIVATE lammps GTest::GMock GTest::GTest Yaml::Yaml)
target_compile_definitions(pair_style PRIVATE TEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB MOL_PAIR_TESTS LIST_DIRECTORIES false mol-pair-*.yaml)
foreach(TEST ${MOL_PAIR_TESTS})
  string(REGEX REPLACE "^.*mol-pair-(.*)\.yaml" "MolPairStyle:\\1" TNAME ${TEST})
  add_test(NAME ${TNAME} COMMAND pair_style ${TEST}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endforeach()

file(GLOB MOL_PAIR_TESTS LIST_DIRECTORIES false atomic-pair-*.yaml)
foreach(TEST ${MOL_PAIR_TESTS})
  string(REGEX REPLACE "^.*atomic-pair-(.*)\.yaml" "AtomicPairStyle:\\1" TNAME ${TEST})
  add_test(NAME ${TNAME} COMMAND pair_style ${TEST}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(${TNAME} PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR}")
endforeach()

# bond style tester
add_executable(bond_style bond_style.cpp ${TEST_COMMON_FILES})
target_link_libraries(bond_style PRIVATE lammps GTest::GMock GTest::GTest Yaml::Yaml)
target_compile_definitions(bond_style PRIVATE TEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB BOND_TESTS LIST_DIRECTORIES false bond-*.yaml)
foreach(TEST ${BOND_TESTS})
  string(REGEX REPLACE "^.*bond-(.*)\.yaml" "BondStyle:\\1" TNAME ${TEST})
  add_test(NAME ${TNAME} COMMAND bond_style ${TEST}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endforeach()

