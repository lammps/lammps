# FIXME_OPENMPTARGET - the NVIDIA HPC compiler nvc++ in the OpenMPTarget backend does not pass the perf_tests.
# FIXME_OPENACC - temporarily disabled due to unimplemented features
IF ((KOKKOS_ENABLE_OPENMPTARGET OR KOKKOS_ENABLE_OPENACC) AND KOKKOS_CXX_COMPILER_ID STREQUAL NVHPC)
  RETURN()
ENDIF()
IF (KOKKOS_ENABLE_OPENACC AND KOKKOS_CXX_COMPILER_ID STREQUAL Clang)
  RETURN()
ENDIF()

# all PerformanceTest_* executables are part of regular tests
# TODO: finish converting these into benchmarks (in progress)
IF(KOKKOS_ENABLE_TESTS)
  IF(KOKKOS_ENABLE_CUDA OR KOKKOS_ENABLE_HIP OR KOKKOS_ENABLE_SYCL)
    KOKKOS_ADD_EXECUTABLE (
      PerformanceTest_SharedSpace
      SOURCES test_sharedSpace.cpp
    )
  ENDIF()

  KOKKOS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

  IF(NOT Kokkos_ENABLE_OPENMPTARGET)
  # FIXME OPENMPTARGET needs tasking
    KOKKOS_ADD_EXECUTABLE_AND_TEST(
      PerformanceTest_TaskDag
      SOURCES test_taskdag.cpp
      CATEGORIES PERFORMANCE
    )
  ENDIF()
ENDIF()

IF(NOT Kokkos_ENABLE_BENCHMARKS)
  RETURN()
ENDIF()

IF (KOKKOS_HAS_TRILINOS)
  message(FATAL_ERROR "Benchmarks are not supported when building as part of Trilinos")
ENDIF()

# Find or download google/benchmark library
find_package(benchmark QUIET 1.5.6)
IF(benchmark_FOUND)
  MESSAGE(STATUS "Using google benchmark found in ${benchmark_DIR}")
ELSE()
  message(STATUS "No installed google benchmark found, fetching from GitHub")
  include(FetchContent)
  SET(BENCHMARK_ENABLE_TESTING OFF)

  list(APPEND CMAKE_MESSAGE_INDENT "[benchmark] ")
  FetchContent_Declare(
    googlebenchmark
    DOWNLOAD_EXTRACT_TIMESTAMP FALSE
    URL https://github.com/google/benchmark/archive/refs/tags/v1.6.2.tar.gz
    URL_HASH MD5=14d14849e075af116143a161bc3b927b
  )
  FetchContent_MakeAvailable(googlebenchmark)
  list(POP_BACK CMAKE_MESSAGE_INDENT)

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
  LIST(APPEND BENCHMARK_SOURCES
    BenchmarkMain.cpp
    Benchmark_Context.cpp
  )

  ADD_EXECUTABLE(
    ${BENCHMARK_NAME}
    ${BENCHMARK_SOURCES}
  )
  TARGET_LINK_LIBRARIES(
    ${BENCHMARK_NAME}
    PRIVATE benchmark::benchmark Kokkos::kokkos impl_git_version
  )
  TARGET_INCLUDE_DIRECTORIES(
    ${BENCHMARK_NAME}
    SYSTEM PRIVATE ${benchmark_SOURCE_DIR}/include
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
  PerfTestGramSchmidt.cpp
  PerfTest_CustomReduction.cpp
  PerfTest_ExecSpacePartitioning.cpp
  PerfTestHexGrad.cpp
  PerfTest_MallocFree.cpp
  PerfTest_ViewAllocate.cpp
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
  PerfTest_ViewFill_123.cpp
  PerfTest_ViewFill_45.cpp
  PerfTest_ViewFill_6.cpp
  PerfTest_ViewFill_7.cpp
  PerfTest_ViewFill_8.cpp
  PerfTest_ViewFill_Raw.cpp
  PerfTest_ViewResize_123.cpp
  PerfTest_ViewResize_45.cpp
  PerfTest_ViewResize_6.cpp
  PerfTest_ViewResize_7.cpp
  PerfTest_ViewResize_8.cpp
  PerfTest_ViewResize_Raw.cpp
)

IF(Kokkos_ENABLE_OPENMPTARGET)
# FIXME OPENMPTARGET requires TeamPolicy Reductions and Custom Reduction
  LIST(REMOVE_ITEM BENCHMARK_SOURCES
    PerfTestGramSchmidt.cpp
    PerfTest_CustomReduction.cpp
    PerfTest_ExecSpacePartitioning.cpp
  )
ENDIF()

KOKKOS_ADD_BENCHMARK(
  PerformanceTest_Benchmark
  SOURCES ${BENCHMARK_SOURCES}
)

IF(NOT KOKKOS_ENABLE_CUDA OR KOKKOS_ENABLE_CUDA_LAMBDA)
  KOKKOS_ADD_BENCHMARK(
    Benchmark_Atomic_MinMax
    SOURCES test_atomic_minmax_simple.cpp
  )
ENDIF()

# FIXME_NVHPC
IF(NOT KOKKOS_CXX_COMPILER_ID STREQUAL NVHPC)
  KOKKOS_ADD_BENCHMARK(
    PerformanceTest_Mempool
    SOURCES test_mempool.cpp
  )
ENDIF()

KOKKOS_ADD_BENCHMARK(
  PerformanceTest_Atomic
  SOURCES test_atomic.cpp
)
