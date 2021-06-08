
add_executable(test_lammps_class test_lammps_class.cpp)
target_link_libraries(test_lammps_class PRIVATE lammps GTest::GMockMain GTest::GTest GTest::GMock)
add_test(LammpsClass test_lammps_class)
set_tests_properties(LammpsClass PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")

add_executable(test_input_class test_input_class.cpp)
target_link_libraries(test_input_class PRIVATE lammps GTest::GTest GTest::GTestMain)
add_test(InputClass test_input_class)

add_executable(test_error_class test_error_class.cpp)
target_link_libraries(test_error_class PRIVATE lammps GTest::GMock GTest::GTest)
add_test(ErrorClass test_error_class)
