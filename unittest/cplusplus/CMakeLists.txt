
add_executable(lammps-class lammps-class.cpp)
target_link_libraries(lammps-class PRIVATE lammps GTest::GMockMain GTest::GTest GTest::GMock)
add_test(LammpsClass lammps-class)

add_executable(input-class input-class.cpp)
target_link_libraries(input-class PRIVATE lammps GTest::GTest GTest::GTestMain)
add_test(InputClass input-class)
