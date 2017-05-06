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
INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

IF(Kokkos_ENABLE_Serial)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Serial
    SOURCES
      UnitTestMain.cpp
      serial/TestSerial_Atomics.cpp
      serial/TestSerial_Other.cpp
      serial/TestSerial_Reductions.cpp
      serial/TestSerial_SubView_a.cpp
      serial/TestSerial_SubView_b.cpp
      serial/TestSerial_SubView_c01.cpp
      serial/TestSerial_SubView_c02.cpp
      serial/TestSerial_SubView_c03.cpp
      serial/TestSerial_SubView_c04.cpp
      serial/TestSerial_SubView_c05.cpp
      serial/TestSerial_SubView_c06.cpp
      serial/TestSerial_SubView_c07.cpp
      serial/TestSerial_SubView_c08.cpp
      serial/TestSerial_SubView_c09.cpp
      serial/TestSerial_SubView_c10.cpp
      serial/TestSerial_SubView_c11.cpp
      serial/TestSerial_SubView_c12.cpp
      serial/TestSerial_Team.cpp
      serial/TestSerial_ViewAPI_a.cpp
      serial/TestSerial_ViewAPI_b.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest
  )
ENDIF()

IF(Kokkos_ENABLE_Pthread)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Threads
    SOURCES
      UnitTestMain.cpp
      threads/TestThreads_Atomics.cpp
      threads/TestThreads_Other.cpp
      threads/TestThreads_Reductions.cpp
      threads/TestThreads_SubView_a.cpp
      threads/TestThreads_SubView_b.cpp
      threads/TestThreads_SubView_c01.cpp
      threads/TestThreads_SubView_c02.cpp
      threads/TestThreads_SubView_c03.cpp
      threads/TestThreads_SubView_c04.cpp
      threads/TestThreads_SubView_c05.cpp
      threads/TestThreads_SubView_c06.cpp
      threads/TestThreads_SubView_c07.cpp
      threads/TestThreads_SubView_c08.cpp
      threads/TestThreads_SubView_c09.cpp
      threads/TestThreads_SubView_c10.cpp
      threads/TestThreads_SubView_c11.cpp
      threads/TestThreads_SubView_c12.cpp
      threads/TestThreads_Team.cpp
      threads/TestThreads_ViewAPI_a.cpp
      threads/TestThreads_ViewAPI_b.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest
  )
ENDIF()

IF(Kokkos_ENABLE_OpenMP)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_OpenMP
    SOURCES
      UnitTestMain.cpp
      openmp/TestOpenMP_Atomics.cpp
      openmp/TestOpenMP_Other.cpp
      openmp/TestOpenMP_Reductions.cpp
      openmp/TestOpenMP_SubView_a.cpp
      openmp/TestOpenMP_SubView_b.cpp
      openmp/TestOpenMP_SubView_c01.cpp
      openmp/TestOpenMP_SubView_c02.cpp
      openmp/TestOpenMP_SubView_c03.cpp
      openmp/TestOpenMP_SubView_c04.cpp
      openmp/TestOpenMP_SubView_c05.cpp
      openmp/TestOpenMP_SubView_c06.cpp
      openmp/TestOpenMP_SubView_c07.cpp
      openmp/TestOpenMP_SubView_c08.cpp
      openmp/TestOpenMP_SubView_c09.cpp
      openmp/TestOpenMP_SubView_c10.cpp
      openmp/TestOpenMP_SubView_c11.cpp
      openmp/TestOpenMP_SubView_c12.cpp
      openmp/TestOpenMP_Team.cpp
      openmp/TestOpenMP_ViewAPI_a.cpp
      openmp/TestOpenMP_ViewAPI_b.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest
  )
ENDIF()

IF(Kokkos_ENABLE_Qthreads)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Qthreads
    SOURCES
      UnitTestMain.cpp
      qthreads/TestQthreads_Atomics.cpp
      qthreads/TestQthreads_Other.cpp
      qthreads/TestQthreads_Reductions.cpp
      qthreads/TestQthreads_SubView_a.cpp
      qthreads/TestQthreads_SubView_b.cpp
      qthreads/TestQthreads_SubView_c01.cpp
      qthreads/TestQthreads_SubView_c02.cpp
      qthreads/TestQthreads_SubView_c03.cpp
      qthreads/TestQthreads_SubView_c04.cpp
      qthreads/TestQthreads_SubView_c05.cpp
      qthreads/TestQthreads_SubView_c06.cpp
      qthreads/TestQthreads_SubView_c07.cpp
      qthreads/TestQthreads_SubView_c08.cpp
      qthreads/TestQthreads_SubView_c09.cpp
      qthreads/TestQthreads_SubView_c10.cpp
      qthreads/TestQthreads_SubView_c11.cpp
      qthreads/TestQthreads_SubView_c12.cpp
      qthreads/TestQthreads_Team.cpp
      qthreads/TestQthreads_ViewAPI_a.cpp
      qthreads/TestQthreads_ViewAPI_b.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest
  )
ENDIF()

IF(Kokkos_ENABLE_Cuda)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Cuda
    SOURCES
      UnitTestMain.cpp
      cuda/TestCuda_Atomics.cpp
      cuda/TestCuda_Other.cpp
      cuda/TestCuda_Reductions_a.cpp
      cuda/TestCuda_Reductions_b.cpp
      cuda/TestCuda_Spaces.cpp
      cuda/TestCuda_SubView_a.cpp
      cuda/TestCuda_SubView_b.cpp
      cuda/TestCuda_SubView_c01.cpp
      cuda/TestCuda_SubView_c02.cpp
      cuda/TestCuda_SubView_c03.cpp
      cuda/TestCuda_SubView_c04.cpp
      cuda/TestCuda_SubView_c05.cpp
      cuda/TestCuda_SubView_c06.cpp
      cuda/TestCuda_SubView_c07.cpp
      cuda/TestCuda_SubView_c08.cpp
      cuda/TestCuda_SubView_c09.cpp
      cuda/TestCuda_SubView_c10.cpp
      cuda/TestCuda_SubView_c11.cpp
      cuda/TestCuda_SubView_c12.cpp
      cuda/TestCuda_Team.cpp
      cuda/TestCuda_ViewAPI_a.cpp
      cuda/TestCuda_ViewAPI_b.cpp
      cuda/TestCuda_ViewAPI_c.cpp
      cuda/TestCuda_ViewAPI_d.cpp
      cuda/TestCuda_ViewAPI_e.cpp
      cuda/TestCuda_ViewAPI_f.cpp
      cuda/TestCuda_ViewAPI_g.cpp
      cuda/TestCuda_ViewAPI_h.cpp
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
