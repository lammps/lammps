
add_executable(test_library_open test_library_open.cpp)
target_link_libraries(test_library_open PRIVATE lammps GTest::GTest GTest::GTestMain)
add_test(LibraryOpen test_library_open)

add_executable(test_library_commands test_library_commands.cpp)
target_link_libraries(test_library_commands PRIVATE lammps GTest::GTest GTest::GTestMain)
add_test(LibraryCommands test_library_commands)

add_executable(test_library_properties test_library_properties.cpp)
target_link_libraries(test_library_properties PRIVATE lammps GTest::GTest GTest::GTestMain)
add_test(LibraryProperties test_library_properties)

