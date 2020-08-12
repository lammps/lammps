
add_executable(test_atom_styles test_atom_styles.cpp)
target_link_libraries(test_atom_styles PRIVATE lammps GTest::GMock GTest::GTest)
add_test(NAME AtomStyles COMMAND test_atom_styles WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_pair_unit_convert test_pair_unit_convert.cpp)
target_link_libraries(test_pair_unit_convert PRIVATE lammps GTest::GMock GTest::GTest)
add_test(NAME PairUnitConvert COMMAND test_pair_unit_convert WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(PairUnitConvert PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR}")

add_executable(test_potential_file_reader test_potential_file_reader.cpp)
target_link_libraries(test_potential_file_reader PRIVATE lammps GTest::GMock GTest::GTest)
add_test(NAME PotentialFileReader COMMAND test_potential_file_reader WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(PotentialFileReader PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR}")

if (PKG_MANYBODY)
    add_executable(test_eim_potential_file_reader test_eim_potential_file_reader.cpp)
    target_link_libraries(test_eim_potential_file_reader PRIVATE lammps GTest::GMock GTest::GTest)
    add_test(NAME EIMPotentialFileReader COMMAND test_eim_potential_file_reader WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(EIMPotentialFileReader PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR}")
endif()

add_executable(test_file_operations test_file_operations.cpp)
target_link_libraries(test_file_operations PRIVATE lammps GTest::GMock GTest::GTest)
add_test(NAME FileOperations COMMAND test_file_operations WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
