# FIXME_OPENMPTARGET - the NVIDIA HPC compiler nvc++ in the OpenMPTarget backend does not pass the perf_tests.
# FIXME_OPENACC - temporarily disabled due to unimplemented features
if((KOKKOS_ENABLE_OPENMPTARGET OR KOKKOS_ENABLE_OPENACC) AND KOKKOS_CXX_COMPILER_ID STREQUAL NVHPC)
  return()
endif()
if(KOKKOS_ENABLE_OPENACC AND KOKKOS_CXX_COMPILER_ID STREQUAL Clang)
  return()
endif()

# all PerformanceTest_* executables are part of regular tests
# TODO: finish converting these into benchmarks (in progress)
if(KOKKOS_ENABLE_TESTS)
  if(KOKKOS_ENABLE_CUDA OR KOKKOS_ENABLE_HIP OR KOKKOS_ENABLE_SYCL)
    kokkos_add_executable(PerformanceTest_SharedSpace SOURCES test_sharedSpace.cpp)
  endif()

  kokkos_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
endif()

if(NOT Kokkos_ENABLE_BENCHMARKS)
  return()
endif()

# Find or download google/benchmark library
find_package(benchmark QUIET 1.8.3)
if(benchmark_FOUND)
  message(STATUS "Using google benchmark found in ${benchmark_DIR}")
else()
  message(STATUS "No installed google benchmark found, fetching from GitHub")
  include(FetchContent)
  set(BENCHMARK_ENABLE_TESTING OFF)

  list(APPEND CMAKE_MESSAGE_INDENT "[benchmark] ")
  FetchContent_Declare(
    googlebenchmark
    DOWNLOAD_EXTRACT_TIMESTAMP FALSE
    URL https://github.com/google/benchmark/archive/refs/tags/v1.8.3.tar.gz
    URL_HASH MD5=7b93dd03670665684f1b2e9b70ad17fe
  )
  FetchContent_MakeAvailable(googlebenchmark)
  list(POP_BACK CMAKE_MESSAGE_INDENT)

  # Suppress clang-tidy diagnostics on code that we do not have control over
  if(CMAKE_CXX_CLANG_TIDY)
    set_target_properties(benchmark PROPERTIES CXX_CLANG_TIDY "")
  endif()

  target_compile_options(benchmark PRIVATE -w)
  target_compile_options(benchmark_main PRIVATE -w)
endif()

function(KOKKOS_ADD_BENCHMARK NAME)
  cmake_parse_arguments(BENCHMARK "" "" "SOURCES" ${ARGN})
  if(DEFINED BENCHMARK_UNPARSED_ARGUMENTS)
    message(WARNING "Unexpected arguments when adding a benchmark: " ${BENCHMARK_UNPARSED_ARGUMENTS})
  endif()

  set(BENCHMARK_NAME Kokkos_${NAME})
  list(APPEND BENCHMARK_SOURCES BenchmarkMain.cpp Benchmark_Context.cpp)

  add_executable(${BENCHMARK_NAME} ${BENCHMARK_SOURCES})
  target_link_libraries(${BENCHMARK_NAME} PRIVATE benchmark::benchmark Kokkos::kokkos impl_git_version)
  target_include_directories(${BENCHMARK_NAME} SYSTEM PRIVATE ${benchmark_SOURCE_DIR}/include)

  foreach(SOURCE_FILE ${BENCHMARK_SOURCES})
    set_source_files_properties(${SOURCE_FILE} PROPERTIES LANGUAGE ${KOKKOS_COMPILE_LANGUAGE})
  endforeach()

  string(TIMESTAMP BENCHMARK_TIME "%Y-%m-%d_T%H-%M-%S" UTC)
  set(BENCHMARK_ARGS --benchmark_counters_tabular=true --benchmark_out=${BENCHMARK_NAME}_${BENCHMARK_TIME}.json)

  add_test(NAME ${BENCHMARK_NAME} COMMAND ${BENCHMARK_NAME} ${BENCHMARK_ARGS})
endfunction()

set(BENCHMARK_SOURCES
    PerfTestGramSchmidt.cpp
    PerfTest_CustomReduction.cpp
    PerfTest_ExecSpacePartitioning.cpp
    PerfTest_Gemv.cpp
    PerfTestHexGrad.cpp
    PerfTest_MallocFree.cpp
    PerfTest_Stream.cpp
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

if(Kokkos_ENABLE_OPENMPTARGET)
  # FIXME OPENMPTARGET requires TeamPolicy Reductions and Custom Reduction
  list(REMOVE_ITEM BENCHMARK_SOURCES PerfTestGramSchmidt.cpp PerfTest_CustomReduction.cpp
       PerfTest_ExecSpacePartitioning.cpp
  )
endif()

kokkos_add_benchmark(PerformanceTest_Benchmark SOURCES ${BENCHMARK_SOURCES})

kokkos_add_benchmark(Benchmark_Atomic_MinMax SOURCES test_atomic_minmax_simple.cpp)

kokkos_add_benchmark(
  PerformanceTest_ViewFirstTouch SOURCES PerfTest_ViewFirstTouch_Initialize.cpp PerfTest_ViewFirstTouch_ParallelFor.cpp
  PerfTest_ViewFirstTouch_DeepCopy.cpp
)

# FIXME_NVHPC
if(NOT KOKKOS_CXX_COMPILER_ID STREQUAL NVHPC)
  kokkos_add_benchmark(PerformanceTest_Mempool SOURCES test_mempool.cpp)
endif()

kokkos_add_benchmark(PerformanceTest_Atomic SOURCES test_atomic.cpp)
kokkos_add_benchmark(PerformanceTest_Reduction SOURCES test_reduction.cpp)
