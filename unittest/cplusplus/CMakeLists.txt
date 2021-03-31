
add_executable(test_lammps_class test_lammps_class.cpp)
target_link_libraries(test_lammps_class PRIVATE lammps GTest::GMockMain GTest::GTest GTest::GMock)
add_test(LammpsClass test_lammps_class)

add_executable(test_input_class test_input_class.cpp)
target_link_libraries(test_input_class PRIVATE lammps GTest::GTest GTest::GTestMain)
add_test(InputClass test_input_class)
