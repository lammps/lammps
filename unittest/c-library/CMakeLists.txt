
add_executable(library-open library-open.cpp)
target_link_libraries(library-open PRIVATE lammps ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})
add_test(LibraryOpen library-open)

add_executable(library-commands library-commands.cpp)
target_link_libraries(library-commands PRIVATE lammps ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})
add_test(LibraryCommands library-commands)

