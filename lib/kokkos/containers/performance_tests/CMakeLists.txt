kokkos_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkos_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
kokkos_include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)

foreach(Tag Threads;OpenMP;Cuda;HPX;HIP)
  string(TOUPPER ${Tag} DEVICE)
  string(TOLOWER ${Tag} dir)

  if(Kokkos_ENABLE_${DEVICE})
    message(STATUS "Sources Test${Tag}.cpp")

    set(SOURCES TestMain.cpp Test${Tag}.cpp)

    kokkos_add_executable_and_test(ContainersPerformanceTest_${Tag} SOURCES ${SOURCES})
  endif()
endforeach()
