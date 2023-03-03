
#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
#INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})


# warning: PerfTest_CustomReduction.cpp uses
# ../../algorithms/src/Kokkos_Random.hpp
# we'll just allow it to be included, but note
# that in TriBITS KokkosAlgorithms can be disabled...
#INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/../../algorithms/src")

# FIXME_OPENMPTARGET - the NVIDIA HPC compiler nvc++ in the OpenMPTarget backend does not pass the perf_tests.
# FIXME_OPENACC - temporarily disabled due to unimplemented features
IF ((KOKKOS_ENABLE_OPENMPTARGET OR KOKKOS_ENABLE_OPENACC) AND KOKKOS_CXX_COMPILER_ID STREQUAL NVHPC)
  RETURN()
ENDIF()


SET(SOURCES
  PerfTestMain.cpp
  PerfTestGramSchmidt.cpp
  PerfTestHexGrad.cpp
  PerfTest_CustomReduction.cpp
  PerfTest_ExecSpacePartitioning.cpp
  PerfTest_ViewAllocate.cpp
  PerfTest_ViewFill_123.cpp
  PerfTest_ViewFill_45.cpp
  PerfTest_ViewFill_6.cpp
  PerfTest_ViewFill_7.cpp
  PerfTest_ViewFill_8.cpp
  PerfTest_ViewResize_123.cpp
  PerfTest_ViewResize_45.cpp
  PerfTest_ViewResize_6.cpp
  PerfTest_ViewResize_7.cpp
  PerfTest_ViewResize_8.cpp
  )

IF(Kokkos_ENABLE_OPENMPTARGET)
# FIXME OPENMPTARGET requires TeamPolicy Reductions and Custom Reduction
  LIST(REMOVE_ITEM SOURCES
    PerfTestGramSchmidt.cpp
    PerfTest_CustomReduction.cpp
    PerfTest_ExecSpacePartitioning.cpp
  )
ENDIF()

IF(KOKKOS_ENABLE_CUDA OR KOKKOS_ENABLE_HIP OR KOKKOS_ENABLE_SYCL)
  KOKKOS_ADD_EXECUTABLE (
    PerformanceTest_SharedSpace
    SOURCES test_sharedSpace.cpp
  )
ENDIF()

# Per #374, we always want to build this test, but we only want to run
# it as a PERFORMANCE test.  That's why we separate building the test
# from running the test.

#leave these as basic includes for now
#I don't need anything transitive
KOKKOS_INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/../../algorithms/src")
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

# This test currently times out for MSVC
IF(NOT KOKKOS_CXX_COMPILER_ID STREQUAL "MSVC")
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    PerfTestExec
    SOURCES ${SOURCES}
    CATEGORIES PERFORMANCE
  )
ENDIF()

KOKKOS_ADD_EXECUTABLE_AND_TEST(
  PerformanceTest_Atomic
  SOURCES test_atomic.cpp
  CATEGORIES PERFORMANCE
)

IF(NOT KOKKOS_ENABLE_CUDA OR KOKKOS_ENABLE_CUDA_LAMBDA)
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    PerformanceTest_Atomic_MinMax
    SOURCES test_atomic_minmax_simple.cpp
    CATEGORIES PERFORMANCE
  )
ENDIF()

# FIXME_NVHPC
IF(NOT KOKKOS_CXX_COMPILER_ID STREQUAL NVHPC)
KOKKOS_ADD_EXECUTABLE_AND_TEST(
  PerformanceTest_Mempool
  SOURCES test_mempool.cpp
  CATEGORIES PERFORMANCE
)
ENDIF()

IF(NOT Kokkos_ENABLE_OPENMPTARGET)
# FIXME OPENMPTARGET needs tasking
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    PerformanceTest_TaskDag
    SOURCES test_taskdag.cpp
    CATEGORIES PERFORMANCE
  )
ENDIF()


IF(NOT Kokkos_ENABLE_BENCHMARKS)
  RETURN()
ENDIF()

IF (KOKKOS_HAS_TRILINOS)
  message(FATAL_ERROR "Benchmarks are not supported when building as part of Trilinos")
ENDIF()

find_package(benchmark QUIET)
IF(benchmark_FOUND)
  MESSAGE(STATUS "Using google benchmark found in ${benchmark_DIR}")
ELSE()
  message(STATUS "No installed google benchmark found, fetching from GitHub")
  include(FetchContent)
  SET(BENCHMARK_ENABLE_TESTING OFF)

  list(APPEND CMAKE_MESSAGE_INDENT "    ")
  FetchContent_Declare(
    googlebenchmark
    URL https://github.com/google/benchmark/archive/refs/tags/v1.6.2.tar.gz
    URL_HASH MD5=14d14849e075af116143a161bc3b927b
  )
  FetchContent_MakeAvailable(googlebenchmark)
  list(POP_BACK CMAKE_MESSAGE_INDENT)

  include_directories(${benchmark_SOURCE_DIR}/include)

  # Suppress clang-tidy diagnostics on code that we do not have control over
  IF(CMAKE_CXX_CLANG_TIDY)
    SET_TARGET_PROPERTIES(benchmark PROPERTIES CXX_CLANG_TIDY "")
  ENDIF()

  target_compile_options(benchmark PRIVATE -w)
  target_compile_options(benchmark_main PRIVATE -w)
ENDIF()


FUNCTION(KOKKOS_ADD_BENCHMARK NAME)
  CMAKE_PARSE_ARGUMENTS(
    BENCHMARK
    ""
    ""
    "SOURCES"
    ${ARGN}
  )
  IF(DEFINED BENCHMARK_UNPARSED_ARGUMENTS)
    MESSAGE(
      WARNING
      "Unexpected arguments when adding a benchmark: "
      ${BENCHMARK_UNPARSED_ARGUMENTS}
    )
  ENDIF()

  SET(BENCHMARK_NAME ${PACKAGE_NAME}_${NAME})

  ADD_EXECUTABLE(
    ${BENCHMARK_NAME}
    ${BENCHMARK_SOURCES}
  )
  TARGET_LINK_LIBRARIES(
    ${BENCHMARK_NAME}
    PRIVATE benchmark::benchmark Kokkos::kokkos impl_git_version
  )
  FOREACH(SOURCE_FILE ${BENCHMARK_SOURCES})
    SET_SOURCE_FILES_PROPERTIES(
      ${SOURCE_FILE}
      PROPERTIES LANGUAGE ${KOKKOS_COMPILE_LANGUAGE}
    )
  ENDFOREACH()

  STRING(TIMESTAMP BENCHMARK_TIME "%Y-%m-%d_T%H-%M-%S" UTC)
  SET(
    BENCHMARK_ARGS
    --benchmark_counters_tabular=true
    --benchmark_out=${BENCHMARK_NAME}_${BENCHMARK_TIME}.json
  )

  ADD_TEST(
    NAME ${BENCHMARK_NAME}
    COMMAND ${BENCHMARK_NAME} ${BENCHMARK_ARGS}
  )
ENDFUNCTION()

SET(
  BENCHMARK_SOURCES
  BenchmarkMain.cpp
  PerfTest_ViewCopy_a123.cpp
  PerfTest_ViewCopy_b123.cpp
  PerfTest_ViewCopy_c123.cpp
  PerfTest_ViewCopy_d123.cpp
  PerfTest_ViewCopy_a45.cpp
  PerfTest_ViewCopy_b45.cpp
  PerfTest_ViewCopy_c45.cpp
  PerfTest_ViewCopy_d45.cpp
  PerfTest_ViewCopy_a6.cpp
  PerfTest_ViewCopy_b6.cpp
  PerfTest_ViewCopy_c6.cpp
  PerfTest_ViewCopy_d6.cpp
  PerfTest_ViewCopy_a7.cpp
  PerfTest_ViewCopy_b7.cpp
  PerfTest_ViewCopy_c7.cpp
  PerfTest_ViewCopy_d7.cpp
  PerfTest_ViewCopy_a8.cpp
  PerfTest_ViewCopy_b8.cpp
  PerfTest_ViewCopy_c8.cpp
  PerfTest_ViewCopy_d8.cpp
  PerfTest_ViewCopy_Raw.cpp
)

KOKKOS_ADD_BENCHMARK(
  PerformanceTest_Benchmark
  SOURCES ${BENCHMARK_SOURCES}
)
