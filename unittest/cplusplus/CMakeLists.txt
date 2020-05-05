
add_executable(lammps-class lammps-class.cpp)
target_link_libraries(lammps-class PRIVATE lammps GTest::GTest)
add_test(BaseClasses lammps-class)
