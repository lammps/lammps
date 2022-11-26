
add_executable(test_lammps_class test_lammps_class.cpp)
target_link_libraries(test_lammps_class PRIVATE lammps GTest::gmock_main)
add_test(NAME LammpsClass COMMAND test_lammps_class)
set_tests_properties(LammpsClass PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")

add_executable(test_input_class test_input_class.cpp)
target_link_libraries(test_input_class PRIVATE lammps GTest::gtest_main)
add_test(NAME InputClass COMMAND test_input_class)

add_executable(test_advanced_utils test_advanced_utils.cpp)
target_link_libraries(test_advanced_utils PRIVATE lammps GTest::gmock)
add_test(NAME AdvancedUtils COMMAND test_advanced_utils)

add_executable(test_error_class test_error_class.cpp)
target_link_libraries(test_error_class PRIVATE lammps GTest::gmock)
add_test(NAME ErrorClass COMMAND test_error_class)
