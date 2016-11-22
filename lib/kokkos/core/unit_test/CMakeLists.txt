#
# Add test-only library for gtest to be reused by all the subpackages
#

SET(GTEST_SOURCE_DIR ${${PARENT_PACKAGE_NAME}_SOURCE_DIR}/tpls/gtest)

INCLUDE_DIRECTORIES(${GTEST_SOURCE_DIR})
TRIBITS_ADD_LIBRARY(
  kokkos_gtest
  HEADERS ${GTEST_SOURCE_DIR}/gtest/gtest.h
  SOURCES ${GTEST_SOURCE_DIR}/gtest/gtest-all.cc
  TESTONLY
  )

#
# Define the tests
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(Kokkos_ENABLE_Serial)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Serial
    SOURCES UnitTestMain.cpp TestSerial.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest
  )
ENDIF()

IF(Kokkos_ENABLE_Pthread)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Threads
    SOURCES UnitTestMain.cpp TestThreads.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest
  )
ENDIF()

IF(Kokkos_ENABLE_OpenMP)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_OpenMP
    SOURCES UnitTestMain.cpp TestOpenMP.cpp TestOpenMP_a.cpp TestOpenMP_b.cpp TestOpenMP_c.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest
  )
ENDIF()

IF(Kokkos_ENABLE_QTHREAD)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Qthread
    SOURCES UnitTestMain.cpp TestQthread.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest
  )
ENDIF()

IF(Kokkos_ENABLE_Cuda)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Cuda
    SOURCES UnitTestMain.cpp TestCuda.cpp TestCuda_a.cpp TestCuda_b.cpp TestCuda_c.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest
  )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_Default
  SOURCES UnitTestMain.cpp TestDefaultDeviceType.cpp TestDefaultDeviceType_a.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest
)

foreach(INITTESTS_NUM RANGE 1 16)
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_DefaultInit_${INITTESTS_NUM}
  SOURCES UnitTestMain.cpp TestDefaultDeviceTypeInit_${INITTESTS_NUM}.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest
)
endforeach(INITTESTS_NUM)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_HWLOC
  SOURCES UnitTestMain.cpp  TestHWLOC.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest
)

