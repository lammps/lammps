add_executable(test_tokenizer test_tokenizer.cpp)
target_link_libraries(test_tokenizer PRIVATE lammps GTest::GMockMain GTest::GMock GTest::GTest)
add_test(Tokenizer test_tokenizer)

add_executable(test_mempool test_mempool.cpp)
target_link_libraries(test_mempool PRIVATE lammps GTest::GMockMain GTest::GMock GTest::GTest)
add_test(MemPool test_mempool)

add_executable(test_argutils test_argutils.cpp)
target_link_libraries(test_argutils PRIVATE lammps GTest::GMockMain GTest::GMock GTest::GTest)
add_test(ArgUtils test_argutils)

add_executable(test_utils test_utils.cpp)
target_link_libraries(test_utils PRIVATE lammps GTest::GMockMain GTest::GMock GTest::GTest)
add_test(Utils test_utils)
set_tests_properties(Utils PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR}")

add_executable(test_fmtlib test_fmtlib.cpp)
target_link_libraries(test_fmtlib PRIVATE lammps GTest::GMockMain GTest::GMock GTest::GTest)
add_test(FmtLib test_fmtlib)

add_executable(test_math_eigen_impl test_math_eigen_impl.cpp)
target_include_directories(test_math_eigen_impl PRIVATE ${LAMMPS_SOURCE_DIR})
add_test(MathEigen test_math_eigen_impl 10 5)
