# -*- CMake -*- file for testing the direct calls to individual C++ classes

add_executable(test_lammps_class test_lammps_class.cpp)
target_link_libraries(test_lammps_class PRIVATE lammps GTest::GMockMain)
add_test(NAME LammpsClass COMMAND test_lammps_class)
set_tests_properties(LammpsClass PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")

add_executable(test_input_class test_input_class.cpp)
target_link_libraries(test_input_class PRIVATE lammps GTest::GTestMain)
add_test(NAME InputClass COMMAND test_input_class)

add_executable(test_advanced_utils test_advanced_utils.cpp)
target_link_libraries(test_advanced_utils PRIVATE lammps GTest::GMock)
add_test(NAME AdvancedUtils COMMAND test_advanced_utils)

add_executable(test_error_class test_error_class.cpp)
target_link_libraries(test_error_class PRIVATE lammps GTest::GMock)
add_test(NAME ErrorClass COMMAND test_error_class)
