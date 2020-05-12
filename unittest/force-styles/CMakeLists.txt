
add_executable(mol-pair-style mol-pair-style.cpp)
target_link_libraries(mol-pair-style PRIVATE lammps ${GTEST_LIBRARY})
add_test(MolPairStyle mol-pair-style)

file(INSTALL in.fourmol data.fourmol DESTINATION .)
