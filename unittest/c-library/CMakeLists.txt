
add_executable(library-create library-create.cpp)
target_link_libraries(library-create PRIVATE lammps GTest::GTest)
add_test(LibraryInterface library-create)
