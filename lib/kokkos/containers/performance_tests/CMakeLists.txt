
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src )

IF(Kokkos_ENABLE_CUDA)
   SET(SOURCES
     TestMain.cpp
     TestCuda.cpp
     )

   KOKKOS_ADD_EXECUTABLE_AND_TEST( PerformanceTest_Cuda
                               SOURCES ${SOURCES}
                             )
ENDIF()

IF(Kokkos_ENABLE_PTHREAD)
   SET(SOURCES
     TestMain.cpp
     TestThreads.cpp
   )
   KOKKOS_ADD_EXECUTABLE_AND_TEST( PerformanceTest_Threads
                               SOURCES ${SOURCES}
                             )
ENDIF()

IF(Kokkos_ENABLE_OPENMP)
   SET(SOURCES
     TestMain.cpp
     TestOpenMP.cpp
   )
   KOKKOS_ADD_EXECUTABLE_AND_TEST( PerformanceTest_OpenMP
                               SOURCES ${SOURCES}
                             )
ENDIF()

IF(Kokkos_ENABLE_HPX)
   SET(SOURCES
     TestMain.cpp
     TestHPX.cpp
   )
   KOKKOS_ADD_EXECUTABLE_AND_TEST( PerformanceTest_HPX
                               SOURCES ${SOURCES}
                             )
ENDIF()

