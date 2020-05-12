
add_executable(lammps-class lammps-class.cpp)
target_link_libraries(lammps-class PRIVATE lammps ${GMOCK_MAIN_LIBRARY} ${GTEST_LIBRARY} ${GMOCK_LIBRARY})
add_test(LammpsClass lammps-class)

add_executable(input-class input-class.cpp)
target_link_libraries(input-class PRIVATE lammps ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})
add_test(InputClass input-class)
