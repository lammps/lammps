
add_executable(library-open library-open.cpp)
target_link_libraries(library-open PRIVATE lammps GTest::GTest GTest::GTestMain)
add_test(LibraryOpen library-open)

add_executable(library-commands library-commands.cpp)
target_link_libraries(library-commands PRIVATE lammps GTest::GTest GTest::GTestMain)
add_test(LibraryCommands library-commands)

add_executable(library-properties library-properties.cpp)
target_link_libraries(library-properties PRIVATE lammps GTest::GTest GTest::GTestMain)
add_test(LibraryProperties library-properties)

