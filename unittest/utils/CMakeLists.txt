add_executable(test_tokenizer test_tokenizer.cpp)
target_link_libraries(test_tokenizer PRIVATE lammps GTest::GMockMain GTest::GMock GTest::GTest)
add_test(Tokenizer test_tokenizer)

add_executable(test_utils test_utils.cpp)
target_link_libraries(test_utils PRIVATE lammps GTest::GMockMain GTest::GMock GTest::GTest)
add_test(Utils test_utils)
