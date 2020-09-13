
add_executable(test_library_open test_library_open.cpp test_main.cpp)
target_link_libraries(test_library_open PRIVATE lammps GTest::GTest GTest::GMock)
add_test(LibraryOpen test_library_open)

add_executable(test_library_commands test_library_commands.cpp test_main.cpp)
target_link_libraries(test_library_commands PRIVATE lammps GTest::GTest GTest::GMock)
add_test(LibraryCommands test_library_commands)

add_executable(test_library_properties test_library_properties.cpp test_main.cpp)
target_link_libraries(test_library_properties PRIVATE lammps GTest::GTest GTest::GMock)
target_compile_definitions(test_library_properties PRIVATE -DTEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR})
add_test(LibraryProperties test_library_properties)
