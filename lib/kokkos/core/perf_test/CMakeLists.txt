
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

IF(NOT KOKKOS_HAS_TRILINOS)
  IF(KOKKOS_SEPARATE_LIBS)
    set(TEST_LINK_TARGETS kokkoscore)
  ELSE()
    set(TEST_LINK_TARGETS kokkos)
  ENDIF()
ENDIF()

# warning: PerfTest_CustomReduction.cpp uses
# ../../algorithms/src/Kokkos_Random.hpp
# we'll just allow it to be included, but note
# that in TriBITS KokkosAlgorithms can be disabled...
INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/../../algorithms/src")

SET(SOURCES
  PerfTestMain.cpp
  PerfTestGramSchmidt.cpp
  PerfTestHexGrad.cpp
  PerfTest_CustomReduction.cpp
  )

# Per #374, we always want to build this test, but we only want to run
# it as a PERFORMANCE test.  That's why we separate building the test
# from running the test.

TRIBITS_ADD_EXECUTABLE(
  PerfTestExec
  SOURCES ${SOURCES}
  COMM serial mpi
  TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
  )

TRIBITS_ADD_TEST(
  PerfTest
  NAME PerfTestExec
  COMM serial mpi
  NUM_MPI_PROCS 1
  CATEGORIES PERFORMANCE
  FAIL_REGULAR_EXPRESSION "  FAILED  "
  )

TRIBITS_ADD_EXECUTABLE(
  PerformanceTest_TaskDAG
  SOURCES test_taskdag.cpp
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  PerformanceTest_TaskDAG
  NAME PerformanceTest_TaskDAG
  COMM serial mpi
  NUM_MPI_PROCS 1
  CATEGORIES PERFORMANCE
  )
