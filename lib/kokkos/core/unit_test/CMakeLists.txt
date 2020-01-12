#
# Add test-only library for gtest to be reused by all the subpackages
#

IF(NOT KOKKOS_HAS_TRILINOS)
  IF(KOKKOS_SEPARATE_LIBS)
    set(TEST_LINK_TARGETS kokkoscore)
  ELSE()
    set(TEST_LINK_TARGETS kokkos)
  ENDIF()
ENDIF()

SET(GTEST_SOURCE_DIR ${${PARENT_PACKAGE_NAME}_SOURCE_DIR}/tpls/gtest)
# TODO get the C++ standard flag from KOKKOS_CXX_STANDARD
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGTEST_HAS_PTHREAD=0")

INCLUDE_DIRECTORIES(${GTEST_SOURCE_DIR})
TRIBITS_ADD_LIBRARY(
  kokkos_gtest
  HEADERS ${GTEST_SOURCE_DIR}/gtest/gtest.h
  SOURCES ${GTEST_SOURCE_DIR}/gtest/gtest-all.cc
  TESTONLY
  )

IF(NOT KOKKOS_HAS_TRILINOS)
target_compile_options(
  kokkos_gtest
  PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${KOKKOS_CXX_FLAGS}>
)
target_link_libraries(
  kokkos_gtest
  PUBLIC ${KOKKOS_LD_FLAGS}
)
ENDIF()

#
# Define the tests
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

IF(Kokkos_ENABLE_Serial)
  IF(KOKKOS_SEPARATE_TESTS)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_Serial_Atomics
      SOURCES
        UnitTestMainInit.cpp
        serial/TestSerial_AtomicOperations_int.cpp
        serial/TestSerial_AtomicOperations_unsignedint.cpp
        serial/TestSerial_AtomicOperations_longint.cpp
        serial/TestSerial_AtomicOperations_unsignedlongint.cpp
        serial/TestSerial_AtomicOperations_longlongint.cpp
        serial/TestSerial_AtomicOperations_double.cpp
        serial/TestSerial_AtomicOperations_float.cpp
        serial/TestSerial_AtomicOperations_complexdouble.cpp
        serial/TestSerial_AtomicOperations_complexfloat.cpp
        serial/TestSerial_AtomicViews.cpp
        serial/TestSerial_Atomics.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_Serial_SubView
      SOURCES
        UnitTestMainInit.cpp
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
        serial/TestSerial_SubView_c13.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_Serial_ViewAPI
      SOURCES
        UnitTestMainInit.cpp
        serial/TestSerial_ViewAPI_a.cpp
        serial/TestSerial_ViewAPI_b.cpp
        serial/TestSerial_ViewAPI_c.cpp
        serial/TestSerial_ViewAPI_d.cpp
        serial/TestSerial_ViewAPI_e.cpp
        serial/TestSerial_ViewOfClass.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_Serial_ViewMapping
      SOURCES
        UnitTestMainInit.cpp
        serial/TestSerial_ViewMapping_a.cpp
        serial/TestSerial_ViewMapping_b.cpp
        serial/TestSerial_ViewMapping_subview.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_Serial_Reducers
      SOURCES
        UnitTestMainInit.cpp
        serial/TestSerial_Reductions.cpp
        serial/TestSerial_Reducers_a.cpp
        serial/TestSerial_Reducers_b.cpp
        serial/TestSerial_Reducers_c.cpp
        serial/TestSerial_Reducers_d.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_Serial_MDRange
      SOURCES
        UnitTestMainInit.cpp
        serial/TestSerial_MDRange_a.cpp
        serial/TestSerial_MDRange_b.cpp
        serial/TestSerial_MDRange_c.cpp
        serial/TestSerial_MDRange_d.cpp
        serial/TestSerial_MDRange_e.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_Serial_Team
      SOURCES
        UnitTestMainInit.cpp
        serial/TestSerial_Team.cpp
        serial/TestSerial_TeamReductionScan.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_Serial_Tasking
      SOURCES
        UnitTestMainInit.cpp
        serial/TestSerial_Task.cpp
        serial/TestSerial_WorkGraph.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_Serial_Misc
      SOURCES
        UnitTestMainInit.cpp
        serial/TestSerial_Complex.cpp
        serial/TestSerial_Init.cpp
        serial/TestSerial_Other.cpp
        serial/TestSerial_RangePolicy.cpp
        serial/TestSerial_Scan.cpp
        serial/TestSerial_SharedAlloc.cpp
        serial/TestSerial_Crs.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
  ELSE()
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_Serial
      SOURCES
        UnitTestMainInit.cpp
        serial/TestSerial_AtomicOperations_int.cpp
        serial/TestSerial_AtomicOperations_unsignedint.cpp
        serial/TestSerial_AtomicOperations_longint.cpp
        serial/TestSerial_AtomicOperations_unsignedlongint.cpp
        serial/TestSerial_AtomicOperations_longlongint.cpp
        serial/TestSerial_AtomicOperations_double.cpp
        serial/TestSerial_AtomicOperations_float.cpp
        serial/TestSerial_AtomicOperations_complexdouble.cpp
        serial/TestSerial_AtomicOperations_complexfloat.cpp
        serial/TestSerial_AtomicViews.cpp
        serial/TestSerial_Atomics.cpp
        serial/TestSerial_Complex.cpp
        serial/TestSerial_Init.cpp
        serial/TestSerial_MDRange_a.cpp
        serial/TestSerial_MDRange_b.cpp
        serial/TestSerial_MDRange_c.cpp
        serial/TestSerial_MDRange_d.cpp
        serial/TestSerial_MDRange_e.cpp
        serial/TestSerial_Other.cpp
        serial/TestSerial_RangePolicy.cpp
        serial/TestSerial_Reductions.cpp
        serial/TestSerial_Reducers_a.cpp
        serial/TestSerial_Reducers_b.cpp
        serial/TestSerial_Reducers_c.cpp
        serial/TestSerial_Reducers_d.cpp
        serial/TestSerial_Scan.cpp
        serial/TestSerial_SharedAlloc.cpp
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
        serial/TestSerial_SubView_c13.cpp
        serial/TestSerial_Task.cpp
        serial/TestSerial_Team.cpp
        serial/TestSerial_TeamReductionScan.cpp
        serial/TestSerial_TeamScratch.cpp
        serial/TestSerial_ViewAPI_a.cpp
        serial/TestSerial_ViewAPI_b.cpp
        serial/TestSerial_ViewAPI_c.cpp
        serial/TestSerial_ViewAPI_d.cpp
        serial/TestSerial_ViewAPI_e.cpp
        serial/TestSerial_ViewMapping_a.cpp
        serial/TestSerial_ViewMapping_b.cpp
        serial/TestSerial_ViewMapping_subview.cpp
        serial/TestSerial_ViewOfClass.cpp
        serial/TestSerial_Crs.cpp
        serial/TestSerial_WorkGraph.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
  ENDIF()
ENDIF()

IF(Kokkos_ENABLE_Pthread)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Threads
    SOURCES
      UnitTestMainInit.cpp
      threads/TestThreads_AtomicOperations_int.cpp
      threads/TestThreads_AtomicOperations_unsignedint.cpp
      threads/TestThreads_AtomicOperations_longint.cpp
      threads/TestThreads_AtomicOperations_unsignedlongint.cpp
      threads/TestThreads_AtomicOperations_longlongint.cpp
      threads/TestThreads_AtomicOperations_double.cpp
      threads/TestThreads_AtomicOperations_float.cpp
      threads/TestThreads_AtomicOperations_complexdouble.cpp
      threads/TestThreads_AtomicOperations_complexfloat.cpp
      threads/TestThreads_AtomicViews.cpp
      threads/TestThreads_Atomics.cpp
      threads/TestThreads_Complex.cpp
      threads/TestThreads_Init.cpp
      threads/TestThreads_MDRange_a.cpp
      threads/TestThreads_MDRange_b.cpp
      threads/TestThreads_MDRange_c.cpp
      threads/TestThreads_MDRange_d.cpp
      threads/TestThreads_MDRange_e.cpp
      threads/TestThreads_Other.cpp
      threads/TestThreads_RangePolicy.cpp
      threads/TestThreads_Reductions.cpp
      threads/TestThreads_Reducers_a.cpp
      threads/TestThreads_Reducers_b.cpp
      threads/TestThreads_Reducers_c.cpp
      threads/TestThreads_Reducers_d.cpp
      threads/TestThreads_Scan.cpp
      threads/TestThreads_SharedAlloc.cpp
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
      threads/TestThreads_SubView_c13.cpp
      threads/TestThreads_Team.cpp
      threads/TestThreads_TeamReductionScan.cpp
      threads/TestThreads_TeamScratch.cpp
      threads/TestThreads_ViewAPI_a.cpp
      threads/TestThreads_ViewAPI_b.cpp
      threads/TestThreads_ViewAPI_c.cpp
      threads/TestThreads_ViewAPI_d.cpp
      threads/TestThreads_ViewAPI_e.cpp
      threads/TestThreads_ViewMapping_a.cpp
      threads/TestThreads_ViewMapping_b.cpp
      threads/TestThreads_ViewMapping_subview.cpp
      threads/TestThreads_ViewOfClass.cpp
      threads/TestThreads_Crs.cpp
      threads/TestThreads_WorkGraph.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
  )
ENDIF()

IF(Kokkos_ENABLE_OpenMP)
  IF(KOKKOS_SEPARATE_TESTS)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_OpenMP_Atomics
      SOURCES
        UnitTestMainInit.cpp
        openmp/TestOpenMP_AtomicOperations_int.cpp
        openmp/TestOpenMP_AtomicOperations_unsignedint.cpp
        openmp/TestOpenMP_AtomicOperations_longint.cpp
        openmp/TestOpenMP_AtomicOperations_unsignedlongint.cpp
        openmp/TestOpenMP_AtomicOperations_longlongint.cpp
        openmp/TestOpenMP_AtomicOperations_double.cpp
        openmp/TestOpenMP_AtomicOperations_float.cpp
        openmp/TestOpenMP_AtomicOperations_complexdouble.cpp
        openmp/TestOpenMP_AtomicOperations_complexfloat.cpp
        openmp/TestOpenMP_AtomicViews.cpp
        openmp/TestOpenMP_Atomics.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_OpenMP_SubView
      SOURCES
        UnitTestMainInit.cpp
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
        openmp/TestOpenMP_SubView_c13.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_OpenMP_ViewAPI
      SOURCES
        UnitTestMainInit.cpp
        openmp/TestOpenMP_ViewAPI_a.cpp
        openmp/TestOpenMP_ViewAPI_b.cpp
        openmp/TestOpenMP_ViewAPI_c.cpp
        openmp/TestOpenMP_ViewAPI_d.cpp
        openmp/TestOpenMP_ViewAPI_e.cpp
        openmp/TestOpenMP_ViewOfClass.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_OpenMP_ViewMapping
      SOURCES
        UnitTestMainInit.cpp
        openmp/TestOpenMP_ViewMapping_a.cpp
        openmp/TestOpenMP_ViewMapping_b.cpp
        openmp/TestOpenMP_ViewMapping_subview.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_OpenMP_Reducers
      SOURCES
        UnitTestMainInit.cpp
        openmp/TestOpenMP_Reductions.cpp
        openmp/TestOpenMP_Reducers_a.cpp
        openmp/TestOpenMP_Reducers_b.cpp
        openmp/TestOpenMP_Reducers_c.cpp
        openmp/TestOpenMP_Reducers_d.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_OpenMP_MDRange
      SOURCES
        UnitTestMainInit.cpp
        openmp/TestOpenMP_MDRange_a.cpp
        openmp/TestOpenMP_MDRange_b.cpp
        openmp/TestOpenMP_MDRange_c.cpp
        openmp/TestOpenMP_MDRange_d.cpp
        openmp/TestOpenMP_MDRange_e.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_OpenMP_Team
      SOURCES
        UnitTestMainInit.cpp
        openmp/TestOpenMP_Team.cpp
        openmp/TestOpenMP_TeamReductionScan.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_OpenMP_Tasking
      SOURCES
        UnitTestMainInit.cpp
        openmp/TestOpenMP_Task.cpp
        openmp/TestOpenMP_WorkGraph.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_OpenMP_Misc
      SOURCES
        UnitTestMainInit.cpp
        openmp/TestOpenMP_Complex.cpp
        openmp/TestOpenMP_Init.cpp
        openmp/TestOpenMP_Other.cpp
        openmp/TestOpenMP_RangePolicy.cpp
        openmp/TestOpenMP_Scan.cpp
        openmp/TestOpenMP_SharedAlloc.cpp
        openmp/TestOpenMP_Crs.cpp
        openmp/TestOpenMP_UniqueToken.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_OpenMPInterOp
      SOURCES
        UnitTestMain.cpp
        openmp/TestOpenMP_InterOp.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
  ELSE()
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_OpenMP
      SOURCES
        UnitTestMainInit.cpp
        openmp/TestOpenMP_AtomicOperations_int.cpp
        openmp/TestOpenMP_AtomicOperations_unsignedint.cpp
        openmp/TestOpenMP_AtomicOperations_longint.cpp
        openmp/TestOpenMP_AtomicOperations_unsignedlongint.cpp
        openmp/TestOpenMP_AtomicOperations_longlongint.cpp
        openmp/TestOpenMP_AtomicOperations_double.cpp
        openmp/TestOpenMP_AtomicOperations_float.cpp
        openmp/TestOpenMP_AtomicOperations_complexdouble.cpp
        openmp/TestOpenMP_AtomicOperations_complexfloat.cpp
        openmp/TestOpenMP_AtomicViews.cpp
        openmp/TestOpenMP_Atomics.cpp
        openmp/TestOpenMP_Complex.cpp
        openmp/TestOpenMP_Init.cpp
        openmp/TestOpenMP_MDRange_a.cpp
        openmp/TestOpenMP_MDRange_b.cpp
        openmp/TestOpenMP_MDRange_c.cpp
        openmp/TestOpenMP_MDRange_d.cpp
        openmp/TestOpenMP_MDRange_e.cpp
        openmp/TestOpenMP_Other.cpp
        openmp/TestOpenMP_RangePolicy.cpp
        openmp/TestOpenMP_Reductions.cpp
        openmp/TestOpenMP_Reducers_a.cpp
        openmp/TestOpenMP_Reducers_b.cpp
        openmp/TestOpenMP_Reducers_c.cpp
        openmp/TestOpenMP_Reducers_d.cpp
        openmp/TestOpenMP_Scan.cpp
        openmp/TestOpenMP_SharedAlloc.cpp
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
        openmp/TestOpenMP_SubView_c13.cpp
        openmp/TestOpenMP_Task.cpp
        openmp/TestOpenMP_Team.cpp
        openmp/TestOpenMP_TeamReductionScan.cpp
        openmp/TestOpenMP_ViewAPI_a.cpp
        openmp/TestOpenMP_ViewAPI_b.cpp
        openmp/TestOpenMP_ViewAPI_c.cpp
        openmp/TestOpenMP_ViewAPI_d.cpp
        openmp/TestOpenMP_ViewAPI_e.cpp
        openmp/TestOpenMP_ViewMapping_a.cpp
        openmp/TestOpenMP_ViewMapping_b.cpp
        openmp/TestOpenMP_ViewMapping_subview.cpp
        openmp/TestOpenMP_ViewOfClass.cpp
        openmp/TestOpenMP_Crs.cpp
        openmp/TestOpenMP_WorkGraph.cpp
        openmp/TestOpenMP_UniqueToken.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_OpenMPInterOp
      SOURCES
        UnitTestMain.cpp
        openmp/TestOpenMP_InterOp.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      FAIL_REGULAR_EXPRESSION "  FAILED  "
      TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
    )
  ENDIF()
ENDIF()

IF(Kokkos_ENABLE_HPX)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_HPX
    SOURCES
      UnitTestMainInit.cpp
      hpx/TestHPX_AtomicOperations_int.cpp
      hpx/TestHPX_AtomicOperations_unsignedint.cpp
      hpx/TestHPX_AtomicOperations_longint.cpp
      hpx/TestHPX_AtomicOperations_unsignedlongint.cpp
      hpx/TestHPX_AtomicOperations_longlongint.cpp
      hpx/TestHPX_AtomicOperations_double.cpp
      hpx/TestHPX_AtomicOperations_float.cpp
      hpx/TestHPX_AtomicViews.cpp
      hpx/TestHPX_Atomics.cpp
      hpx/TestHPX_Complex.cpp
      hpx/TestHPX_Init.cpp
      hpx/TestHPX_MDRange_a.cpp
      hpx/TestHPX_MDRange_b.cpp
      hpx/TestHPX_MDRange_c.cpp
      hpx/TestHPX_MDRange_d.cpp
      hpx/TestHPX_MDRange_e.cpp
      hpx/TestHPX_Other.cpp
      hpx/TestHPX_RangePolicy.cpp
      hpx/TestHPX_Reductions.cpp
      hpx/TestHPX_Reducers_a.cpp
      hpx/TestHPX_Reducers_b.cpp
      hpx/TestHPX_Reducers_c.cpp
      hpx/TestHPX_Reducers_d.cpp
      hpx/TestHPX_Scan.cpp
      hpx/TestHPX_SharedAlloc.cpp
      hpx/TestHPX_SubView_a.cpp
      hpx/TestHPX_SubView_b.cpp
      hpx/TestHPX_SubView_c01.cpp
      hpx/TestHPX_SubView_c02.cpp
      hpx/TestHPX_SubView_c03.cpp
      hpx/TestHPX_SubView_c04.cpp
      hpx/TestHPX_SubView_c05.cpp
      hpx/TestHPX_SubView_c06.cpp
      hpx/TestHPX_SubView_c07.cpp
      hpx/TestHPX_SubView_c08.cpp
      hpx/TestHPX_SubView_c09.cpp
      hpx/TestHPX_SubView_c10.cpp
      hpx/TestHPX_SubView_c11.cpp
      hpx/TestHPX_SubView_c12.cpp
      hpx/TestHPX_SubView_c13.cpp
      hpx/TestHPX_Task.cpp
      hpx/TestHPX_Team.cpp
      hpx/TestHPX_TeamReductionScan.cpp
      hpx/TestHPX_ViewAPI_a.cpp
      hpx/TestHPX_ViewAPI_b.cpp
      hpx/TestHPX_ViewAPI_c.cpp
      hpx/TestHPX_ViewAPI_d.cpp
      hpx/TestHPX_ViewAPI_e.cpp
      hpx/TestHPX_ViewMapping_a.cpp
      hpx/TestHPX_ViewMapping_b.cpp
      hpx/TestHPX_ViewMapping_subview.cpp
      hpx/TestHPX_ViewOfClass.cpp
      hpx/TestHPX_Crs.cpp
      hpx/TestHPX_WorkGraph.cpp
      hpx/TestHPX_UniqueToken.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
  )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_HPXInterOp
    SOURCES
      UnitTestMain.cpp
      hpx/TestHPX_InterOp.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
  )
ENDIF()

IF(Kokkos_ENABLE_Qthreads)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Qthreads
    SOURCES
      UnitTestMainInit.cpp
      qthreads/TestQthreads_Atomics.cpp
      qthreads/TestQthreads_Complex.cpp
      qthreads/TestQthreads_Other.cpp
      qthreads/TestQthreads_Reductions.cpp
      qthreads/TestQthreads_Reducers_a.cpp
      qthreads/TestQthreads_Reducers_b.cpp
      qthreads/TestQthreads_Reducers_c.cpp
      qthreads/TestQthreads_Reducers_d.cpp
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
      qthreads/TestQthreads_SubView_c13.cpp
      qthreads/TestQthreads_Team.cpp
      qthreads/TestQthreads_ViewAPI_a.cpp
      qthreads/TestQthreads_ViewAPI_b.cpp
      qthreads/TestQthreads_ViewAPI_c.cpp
      qthreads/TestQthreads_ViewAPI_d.cpp
      qthreads/TestQthreads_ViewAPI_e.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
  )
ENDIF()

IF(Kokkos_ENABLE_Cuda)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Cuda
    SOURCES
      UnitTestMainInit.cpp
      cuda/TestCudaHostPinned_SharedAlloc.cpp
      cuda/TestCudaHostPinned_ViewCopy.cpp
      cuda/TestCudaHostPinned_ViewAPI_a.cpp
      cuda/TestCudaHostPinned_ViewAPI_b.cpp
      cuda/TestCudaHostPinned_ViewAPI_c.cpp
      cuda/TestCudaHostPinned_ViewAPI_d.cpp
      cuda/TestCudaHostPinned_ViewAPI_e.cpp
      cuda/TestCudaHostPinned_ViewMapping_a.cpp
      cuda/TestCudaHostPinned_ViewMapping_b.cpp
      cuda/TestCudaHostPinned_ViewMapping_subview.cpp
      cuda/TestCudaUVM_SharedAlloc.cpp
      cuda/TestCudaUVM_ViewCopy.cpp
      cuda/TestCudaUVM_ViewAPI_a.cpp
      cuda/TestCudaUVM_ViewAPI_b.cpp
      cuda/TestCudaUVM_ViewAPI_c.cpp
      cuda/TestCudaUVM_ViewAPI_d.cpp
      cuda/TestCudaUVM_ViewAPI_e.cpp
      cuda/TestCudaUVM_ViewMapping_a.cpp
      cuda/TestCudaUVM_ViewMapping_b.cpp
      cuda/TestCudaUVM_ViewMapping_subview.cpp
      cuda/TestCuda_AtomicOperations_int.cpp
      cuda/TestCuda_AtomicOperations_unsignedint.cpp
      cuda/TestCuda_AtomicOperations_longint.cpp
      cuda/TestCuda_AtomicOperations_unsignedlongint.cpp
      cuda/TestCuda_AtomicOperations_longlongint.cpp
      cuda/TestCuda_AtomicOperations_double.cpp
      cuda/TestCuda_AtomicOperations_float.cpp
      cuda/TestCuda_AtomicOperations_complexdouble.cpp
      cuda/TestCuda_AtomicOperations_complexfloat.cpp
      cuda/TestCuda_AtomicViews.cpp
      cuda/TestCuda_Atomics.cpp
      cuda/TestCuda_Complex.cpp
      cuda/TestCuda_Init.cpp
      cuda/TestCuda_MDRange_a.cpp
      cuda/TestCuda_MDRange_b.cpp
      cuda/TestCuda_MDRange_c.cpp
      cuda/TestCuda_MDRange_d.cpp
      cuda/TestCuda_MDRange_e.cpp
      cuda/TestCuda_Other.cpp
      cuda/TestCuda_RangePolicy.cpp
      cuda/TestCuda_Reductions.cpp
      cuda/TestCuda_Reducers_a.cpp
      cuda/TestCuda_Reducers_b.cpp
      cuda/TestCuda_Reducers_c.cpp
      cuda/TestCuda_Reducers_d.cpp
      cuda/TestCuda_Scan.cpp
      cuda/TestCuda_SharedAlloc.cpp
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
      cuda/TestCuda_SubView_c13.cpp
      cuda/TestCuda_Task.cpp
      cuda/TestCuda_Team.cpp
      cuda/TestCuda_TeamReductionScan.cpp
      cuda/TestCuda_TeamScratch.cpp
      cuda/TestCuda_ViewAPI_a.cpp
      cuda/TestCuda_ViewAPI_b.cpp
      cuda/TestCuda_ViewAPI_c.cpp
      cuda/TestCuda_ViewAPI_d.cpp
      cuda/TestCuda_ViewAPI_e.cpp
      cuda/TestCuda_ViewMapping_a.cpp
      cuda/TestCuda_ViewMapping_b.cpp
      cuda/TestCuda_ViewMapping_subview.cpp
      cuda/TestCuda_ViewOfClass.cpp
      cuda/TestCuda_Crs.cpp
      cuda/TestCuda_WorkGraph.cpp
      cuda/TestCuda_UniqueToken.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
  )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_CudaInterOpInit
    SOURCES
      UnitTestMain.cpp
      cuda/TestCuda_InterOp_Init.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
  )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_CudaInterOpStreams
    SOURCES
      UnitTestMain.cpp
      cuda/TestCuda_InterOp_Streams.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
  )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_Default
  SOURCES
    UnitTestMainInit.cpp
    default/TestDefaultDeviceType.cpp
    default/TestDefaultDeviceType_a1.cpp
    default/TestDefaultDeviceType_b1.cpp
    default/TestDefaultDeviceType_c1.cpp
    default/TestDefaultDeviceType_a2.cpp
    default/TestDefaultDeviceType_b2.cpp
    default/TestDefaultDeviceType_c2.cpp
    default/TestDefaultDeviceType_a3.cpp
    default/TestDefaultDeviceType_b3.cpp
    default/TestDefaultDeviceType_c3.cpp
    default/TestDefaultDeviceType_d.cpp
    default/TestDefaultDeviceTypeResize.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_PushFinalizeHook
  SOURCES
    UnitTest_PushFinalizeHook.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "FAILED"
    TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
)

# This test is special, because it passes exactly when it prints the
# message "PASSED: I am the custom std::terminate handler.", AND calls
# std::terminate.  This means that we can't use
# TRIBITS_ADD_EXECUTABLE_AND_TEST.  See GitHub issue #2147.

TRIBITS_ADD_EXECUTABLE( push_finalize_hook_terminate
  SOURCES UnitTest_PushFinalizeHook_terminate.cpp
  TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
)

TRIBITS_ADD_ADVANCED_TEST( UnitTest_PushFinalizeHook_terminate
  TEST_0
    EXEC push_finalize_hook_terminate
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION
      "PASSED: I am the custom std::terminate handler."
    ALWAYS_FAIL_ON_ZERO_RETURN
)

foreach(INITTESTS_NUM RANGE 1 16)
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_DefaultInit_${INITTESTS_NUM}
  SOURCES UnitTestMain.cpp default/TestDefaultDeviceTypeInit_${INITTESTS_NUM}.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
)
endforeach(INITTESTS_NUM)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_HWLOC
  SOURCES UnitTestMain.cpp  TestHWLOC.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_HostBarrier
  SOURCES UnitTestMain.cpp  TestHostBarrier.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  "
    TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
)

#
# Compile-only tests
#
FUNCTION(KOKKOS_ADD_COMPILE_TEST TEST_NAME)

  SET(options LINK_KOKKOS)
  SET(oneValueArgs)
  SET(multiValueArgs)

  CMAKE_PARSE_ARGUMENTS(PARSE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  IF(PARSE_LINK_KOKKOS)
    SET(libs ${TEST_LINK_TARGETS})
  ELSE()
    SET(libs)
  ENDIF()

  TRIBITS_ADD_EXECUTABLE(
    ${TEST_NAME}
    TESTONLY
    COMM serial
    TESTONLYLIBS ${libs}
    ${PARSE_UNPARSED_ARGUMENTS}
  )

  target_compile_options(
    ${PACKAGE_NAME}_${TEST_NAME}
    PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${KOKKOS_CXX_FLAGS}>
  )
  target_link_libraries(
    ${PACKAGE_NAME}_${TEST_NAME}
    PUBLIC ${KOKKOS_LD_FLAGS}
  )

ENDFUNCTION()

