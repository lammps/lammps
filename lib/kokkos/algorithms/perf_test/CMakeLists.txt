# FIXME: The following logic should be moved from here and also from `core/perf_test/CMakeLists.txt` to
# the root `CMakeLists.txt` in the form of a macro
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

  # FIXME: Check whether the following target_compile_options are needed.
  # If so, clarify why.
  target_compile_options(benchmark PRIVATE -w)
  target_compile_options(benchmark_main PRIVATE -w)
endif()

# FIXME: This function should be moved from here and also from `core/perf_test/CMakeLists.txt` to
# the root `CMakeLists.txt`
# FIXME: Could NAME be a one_value_keyword specified in cmake_parse_arguments?
function(KOKKOS_ADD_BENCHMARK NAME)
  cmake_parse_arguments(BENCHMARK "" "" "SOURCES" ${ARGN})
  if(DEFINED BENCHMARK_UNPARSED_ARGUMENTS)
    message(WARNING "Unexpected arguments when adding a benchmark: " ${BENCHMARK_UNPARSED_ARGUMENTS})
  endif()

  set(BENCHMARK_NAME Kokkos_${NAME})
  # FIXME: BenchmarkMain.cpp and Benchmark_Context.cpp should be moved to a common location from which
  # they can be used by all performance tests.
  list(APPEND BENCHMARK_SOURCES ../../core/perf_test/BenchmarkMain.cpp ../../core/perf_test/Benchmark_Context.cpp)

  add_executable(${BENCHMARK_NAME} ${BENCHMARK_SOURCES})
  target_link_libraries(${BENCHMARK_NAME} PRIVATE benchmark::benchmark Kokkos::kokkos impl_git_version)
  target_include_directories(${BENCHMARK_NAME} SYSTEM PRIVATE ${benchmark_SOURCE_DIR}/include)

  # FIXME: This alone will not work. It might need an architecture and standard which need to be defined on target level.
  # It will potentially go away with #7582.
  foreach(SOURCE_FILE ${BENCHMARK_SOURCES})
    set_source_files_properties(${SOURCE_FILE} PROPERTIES LANGUAGE ${KOKKOS_COMPILE_LANGUAGE})
  endforeach()

  string(TIMESTAMP BENCHMARK_TIME "%Y-%m-%d_T%H-%M-%S" UTC)
  set(BENCHMARK_ARGS --benchmark_counters_tabular=true --benchmark_out=${BENCHMARK_NAME}_${BENCHMARK_TIME}.json)

  add_test(NAME ${BENCHMARK_NAME} COMMAND ${BENCHMARK_NAME} ${BENCHMARK_ARGS})
endfunction()

kokkos_add_benchmark(PerformanceTest_InclusiveScan SOURCES test_inclusive_scan.cpp)
kokkos_add_benchmark(PerformanceTest_Random SOURCES test_random.cpp)
