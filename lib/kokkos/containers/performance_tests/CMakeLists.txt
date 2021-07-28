
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src )

foreach(Tag Threads;OpenMP;Cuda;HPX;HIP)
  # Because there is always an exception to the rule
  if(Tag STREQUAL "Threads")
    set(DEVICE "PTHREAD")
  else()
    string(TOUPPER ${Tag} DEVICE)
  endif()
  string(TOLOWER ${Tag} dir)

  if(Kokkos_ENABLE_${DEVICE})
	    message(STATUS "Sources Test${Tag}.cpp")

    set(SOURCES
        TestMain.cpp
        Test${Tag}.cpp
    )

    KOKKOS_ADD_EXECUTABLE_AND_TEST(
      PerformanceTest_${Tag}
      SOURCES ${SOURCES}
    )
  endif()
endforeach()
