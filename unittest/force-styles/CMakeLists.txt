
find_package(Yaml QUIET)
if(NOT Yaml_FOUND)
  message(STATUS "Skipping tests because libyaml is not found")
  return()
endif()

add_executable(mol-pair-style mol-pair-style.cpp)
target_link_libraries(mol-pair-style PRIVATE lammps "GTest::GTest;Yaml::Yaml")
add_test(MolPairStyle mol-pair-style)

file(INSTALL in.fourmol data.fourmol DESTINATION .)
