
find_package(YAML)
if(NOT YAML_FOUND)
  message(STATUS "Skipping tests because libyaml is not found")
  return()
endif()

# pair style tester
add_executable(pair-style pair-style.cpp)
target_link_libraries(pair-style PRIVATE lammps GTest::GMock GTest::GTest Yaml::Yaml)

file(GLOB MOL_PAIR_TESTS LIST_DIRECTORIES false mol-pair-*.yaml)
foreach(TEST ${MOL_PAIR_TESTS})
  string(REGEX REPLACE "^.*mol-pair-(.*)\.yaml" "MolPairStyle:\\1" TNAME ${TEST})
  add_test(NAME ${TNAME} COMMAND pair-style ${TEST}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endforeach()

file(GLOB MOL_PAIR_TESTS LIST_DIRECTORIES false atomic-pair-*.yaml)
foreach(TEST ${MOL_PAIR_TESTS})
  string(REGEX REPLACE "^.*atomic-pair-(.*)\.yaml" "AtomicPairStyle:\\1" TNAME ${TEST})
  add_test(NAME ${TNAME} COMMAND pair-style ${TEST}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(${TNAME} PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR}")
endforeach()
# copy files read at runtime to the working directory
file(COPY in.fourmol data.fourmol in.metal data.metal DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
