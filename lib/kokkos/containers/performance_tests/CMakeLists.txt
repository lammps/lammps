
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src )

SET(SOURCES
  TestMain.cpp 
  TestCuda.cpp
  )

IF(Kokkos_ENABLE_Pthread)
  LIST( APPEND SOURCES TestThreads.cpp)
ENDIF()

IF(Kokkos_ENABLE_OpenMP)
  LIST( APPEND SOURCES TestOpenMP.cpp)
ENDIF()

# Per #374, we always want to build this test, but we only want to run
# it as a PERFORMANCE test.  That's why we separate building the test
# from running the test.

TRIBITS_ADD_EXECUTABLE(
  PerfTestExec
  SOURCES ${SOURCES}
  COMM serial mpi
  TESTONLYLIBS kokkos_gtest
  )

TRIBITS_ADD_TEST(
  PerformanceTest
  NAME PerfTestExec
  COMM serial mpi
  NUM_MPI_PROCS 1
  CATEGORIES PERFORMANCE
  FAIL_REGULAR_EXPRESSION "  FAILED  "
  )
