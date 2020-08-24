
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

add_executable(test_dump_atom test_dump_atom.cpp)
target_link_libraries(test_dump_atom PRIVATE lammps GTest::GMock GTest::GTest)
add_test(NAME DumpAtom COMMAND test_dump_atom WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(DumpAtom PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR}")

add_executable(test_dump_custom test_dump_custom.cpp)
target_link_libraries(test_dump_custom PRIVATE lammps GTest::GMock GTest::GTest)
add_test(NAME DumpCustom COMMAND test_dump_custom WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(DumpCustom PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR}")

add_executable(test_dump_cfg test_dump_cfg.cpp)
target_link_libraries(test_dump_cfg PRIVATE lammps GTest::GMock GTest::GTest)
add_test(NAME DumpCfg COMMAND test_dump_cfg WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(DumpCfg PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR}")

if (BUILD_TOOLS)
    set_tests_properties(DumpAtom PROPERTIES ENVIRONMENT "BINARY2TXT_BINARY=$<TARGET_FILE:binary2txt>")
    set_tests_properties(DumpCustom PROPERTIES ENVIRONMENT "BINARY2TXT_BINARY=$<TARGET_FILE:binary2txt>")
endif()
