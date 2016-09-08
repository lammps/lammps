
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINRARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(SOURCES
  PerfTestMain.cpp
  PerfTestHost.cpp
  PerfTestCuda.cpp
  )

# Per #374, we always want to build this test, but we only want to run
# it as a PERFORMANCE test.  That's why we separate building the test
# from running the test.

TRIBITS_ADD_EXECUTABLE(
  PerfTestExec
  SOURCES ${SOURCES}
  COMM serial mpi
  TESTONLYLIBS kokkos_gtest
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  PerfTest
  NAME PerfTestExec
  COMM serial mpi
  NUM_MPI_PROCS 1
  CATEGORIES PERFORMANCE
  FAIL_REGULAR_EXPRESSION "  FAILED  "
  )
