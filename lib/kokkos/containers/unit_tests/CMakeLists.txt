
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src )

IF(NOT KOKKOS_HAS_TRILINOS)
  IF(KOKKOS_SEPARATE_LIBS)
    set(TEST_LINK_TARGETS kokkoscore)
  ELSE()
    set(TEST_LINK_TARGETS kokkos)
  ENDIF()
ENDIF()

IF(Kokkos_ENABLE_Pthread)
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_Threads
  SOURCES
    UnitTestMain.cpp
    threads/TestThreads_BitSet.cpp
    threads/TestThreads_DualView.cpp
    threads/TestThreads_DynamicView.cpp
    threads/TestThreads_DynRankViewAPI_generic.cpp
    threads/TestThreads_DynRankViewAPI_rank12345.cpp
    threads/TestThreads_DynRankViewAPI_rank67.cpp
    threads/TestThreads_ErrorReporter.cpp
    threads/TestThreads_OffsetView.cpp
    threads/TestThreads_ScatterView.cpp
    threads/TestThreads_StaticCrsGraph.cpp
    threads/TestThreads_UnorderedMap.cpp
    threads/TestThreads_Vector.cpp
    threads/TestThreads_ViewCtorPropEmbeddedDim.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  "
  TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
  )
ENDIF()

IF(Kokkos_ENABLE_Serial)
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_Serial
  SOURCES
    UnitTestMain.cpp
    serial/TestSerial_BitSet.cpp
    serial/TestSerial_DualView.cpp
    serial/TestSerial_DynamicView.cpp
    serial/TestSerial_DynRankViewAPI_generic.cpp
    serial/TestSerial_DynRankViewAPI_rank12345.cpp
    serial/TestSerial_DynRankViewAPI_rank67.cpp
    serial/TestSerial_ErrorReporter.cpp
    serial/TestSerial_OffsetView.cpp
    serial/TestSerial_ScatterView.cpp
    serial/TestSerial_StaticCrsGraph.cpp
    serial/TestSerial_UnorderedMap.cpp
    serial/TestSerial_Vector.cpp
    serial/TestSerial_ViewCtorPropEmbeddedDim.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  "
  TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
  )
ENDIF()

IF(Kokkos_ENABLE_OpenMP)
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_OpenMP
  SOURCES
    UnitTestMain.cpp
    openmp/TestOpenMP_BitSet.cpp
    openmp/TestOpenMP_DualView.cpp
    openmp/TestOpenMP_DynamicView.cpp
    openmp/TestOpenMP_DynRankViewAPI_generic.cpp
    openmp/TestOpenMP_DynRankViewAPI_rank12345.cpp
    openmp/TestOpenMP_DynRankViewAPI_rank67.cpp
    openmp/TestOpenMP_ErrorReporter.cpp
    openmp/TestOpenMP_OffsetView.cpp
    openmp/TestOpenMP_ScatterView.cpp
    openmp/TestOpenMP_StaticCrsGraph.cpp
    openmp/TestOpenMP_UnorderedMap.cpp
    openmp/TestOpenMP_Vector.cpp
    openmp/TestOpenMP_ViewCtorPropEmbeddedDim.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  "
  TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
  )
ENDIF()

IF(Kokkos_ENABLE_HPX)
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_HPX
  SOURCES
    UnitTestMain.cpp
    hpx/TestHPX_BitSet.cpp
    hpx/TestHPX_DualView.cpp
    hpx/TestHPX_DynamicView.cpp
    hpx/TestHPX_DynRankViewAPI_generic.cpp
    hpx/TestHPX_DynRankViewAPI_rank12345.cpp
    hpx/TestHPX_DynRankViewAPI_rank67.cpp
    hpx/TestHPX_ErrorReporter.cpp
    hpx/TestHPX_OffsetView.cpp
    hpx/TestHPX_ScatterView.cpp
    hpx/TestHPX_StaticCrsGraph.cpp
    hpx/TestHPX_UnorderedMap.cpp
    hpx/TestHPX_Vector.cpp
    hpx/TestHPX_ViewCtorPropEmbeddedDim.cpp
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
    UnitTestMain.cpp
    cuda/TestCuda_BitSet.cpp
    cuda/TestCuda_DualView.cpp
    cuda/TestCuda_DynamicView.cpp
    cuda/TestCuda_DynRankViewAPI_generic.cpp
    cuda/TestCuda_DynRankViewAPI_rank12345.cpp
    cuda/TestCuda_DynRankViewAPI_rank67.cpp
    cuda/TestCuda_ErrorReporter.cpp
    cuda/TestCuda_OffsetView.cpp
    cuda/TestCuda_ScatterView.cpp
    cuda/TestCuda_StaticCrsGraph.cpp
    cuda/TestCuda_UnorderedMap.cpp
    cuda/TestCuda_Vector.cpp
    cuda/TestCuda_ViewCtorPropEmbeddedDim.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  "
  TESTONLYLIBS kokkos_gtest ${TEST_LINK_TARGETS}
  )
ENDIF()

