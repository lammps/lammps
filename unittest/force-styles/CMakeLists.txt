
find_package(YAML)
if(NOT YAML_FOUND)
  message(STATUS "Skipping tests because libyaml is not found")
  return()
endif()

add_executable(mol-pair-style mol-pair-style.cpp)
target_link_libraries(mol-pair-style PRIVATE lammps "GTest::GTest;Yaml::Yaml")
file(GLOB MOL_PAIR_TESTS LIST_DIRECTORIES false mol-pair-*.yaml)

foreach(TEST ${MOL_PAIR_TESTS})
  string(REGEX REPLACE "^.*mol-pair-(.*)\.yaml" "MolPairStyle:\\1" TNAME ${TEST})
  add_test(NAME ${TNAME} COMMAND mol-pair-style ${TEST})
endforeach()

file(INSTALL in.fourmol data.fourmol DESTINATION .)
