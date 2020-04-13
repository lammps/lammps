
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src )

IF(Kokkos_ENABLE_CUDA)
   SET(SOURCES
     TestMain.cpp 
     TestCuda.cpp
     )

   KOKKOS_ADD_TEST_EXECUTABLE( PerfTestExec_Cuda
                               SOURCES ${SOURCES}
                             )

   KOKKOS_ADD_TEST( NAME PerformanceTest_Cuda
                    EXE  PerfTestExec_Cuda
                  )
ENDIF()

IF(Kokkos_ENABLE_PTHREAD)
   SET(SOURCES
     TestMain.cpp 
     TestThreads.cpp
   )
   KOKKOS_ADD_TEST_EXECUTABLE( PerfTestExec_Threads
                               SOURCES ${SOURCES}
                             )

   KOKKOS_ADD_TEST( NAME PerformanceTest_Threads
                    EXE  PerfTestExec_Threads
                  )
ENDIF()

IF(Kokkos_ENABLE_OPENMP)
   SET(SOURCES
     TestMain.cpp 
     TestOpenMP.cpp
   )
   KOKKOS_ADD_TEST_EXECUTABLE( PerfTestExec_OpenMP
                               SOURCES ${SOURCES}
                             )

   KOKKOS_ADD_TEST( NAME PerformanceTest_OpenMP
                    EXE  PerfTestExec_OpenMP
                  )
ENDIF()

IF(Kokkos_ENABLE_HPX)
   SET(SOURCES
     TestMain.cpp 
     TestHPX.cpp
   )
   KOKKOS_ADD_TEST_EXECUTABLE( PerfTestExec_HPX
                               SOURCES ${SOURCES}
                             )

   KOKKOS_ADD_TEST( NAME PerformanceTest_HPX
                    EXE  PerfTestExec_HPX
                  )
ENDIF()

