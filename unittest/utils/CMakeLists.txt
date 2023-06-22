
add_executable(test_tokenizer test_tokenizer.cpp)
target_link_libraries(test_tokenizer PRIVATE lammps GTest::GMockMain)
add_test(NAME Tokenizer COMMAND test_tokenizer)

add_executable(test_mempool test_mempool.cpp)
target_link_libraries(test_mempool PRIVATE lammps GTest::GMockMain)
add_test(NAME MemPool COMMAND test_mempool)

add_executable(test_lmptype test_lmptype.cpp)
target_link_libraries(test_lmptype PRIVATE lammps GTest::GMockMain)
add_test(NAME LmpType COMMAND test_lmptype)

add_executable(test_argutils test_argutils.cpp)
target_link_libraries(test_argutils PRIVATE lammps GTest::GMockMain)
add_test(NAME ArgUtils COMMAND test_argutils)

add_executable(test_utils test_utils.cpp)
target_link_libraries(test_utils PRIVATE lammps GTest::GMockMain)
add_test(NAME Utils COMMAND test_utils)

add_executable(test_platform test_platform.cpp)
target_link_libraries(test_platform PRIVATE lammps GTest::GMockMain)
add_test(NAME Platform COMMAND test_platform)

if(PKG_LEPTON)
  add_executable(test_lepton test_lepton.cpp)
  target_link_libraries(test_lepton PRIVATE lepton lammps GTest::GMockMain)
  add_test(NAME Lepton COMMAND test_lepton)
endif()

set_tests_properties(Utils Platform PROPERTIES
  ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR}")

if(BUILD_SHARED_LIBS)
  enable_language(C)
  add_library(testsharedobj MODULE testshared.c)
  set_target_properties(testsharedobj PROPERTIES PREFIX "" WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
  add_library(testsharedlib SHARED testshared.c)
  set_target_properties(testsharedlib PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
  add_dependencies(test_platform testsharedobj testsharedlib)
  target_compile_definitions(test_platform PRIVATE TEST_SHARED_LOAD=1
    TEST_SHARED_LIB=$<TARGET_FILE:testsharedlib>
    TEST_SHARED_OBJ=$<TARGET_FILE:testsharedobj>)
endif()

add_executable(test_fmtlib test_fmtlib.cpp)
target_link_libraries(test_fmtlib PRIVATE lammps GTest::GMockMain)
add_test(NAME FmtLib COMMAND test_fmtlib)

add_executable(test_math_eigen_impl test_math_eigen_impl.cpp)
target_include_directories(test_math_eigen_impl PRIVATE ${LAMMPS_SOURCE_DIR})
add_test(NAME MathEigen COMMAND test_math_eigen_impl 10 5)
