# workaround to allow passing extra arguments to test runs
# through ctest via a TEST_ARGS environment variable
# This can be used to, e.g. reset reference data for individual
# tests from the build folder with "env TEST_ARGS=-u ctest -R sometest"
execute_process(COMMAND ${TEST_EXECUTABLE} ${TEST_INPUT} $ENV{TEST_ARGS} RESULT_VARIABLE rv)
if(NOT "${rv}" STREQUAL "0")
  message(FATAL_ERROR "Test ${TEST_NAME} failed with status ${rv}")
endif()
