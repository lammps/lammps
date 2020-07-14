
#Leave these here for now - I don't need transitive deps anyway
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src )


SET(GTEST_SOURCE_DIR ${${PARENT_PACKAGE_NAME}_SOURCE_DIR}/tpls/gtest)
KOKKOS_INCLUDE_DIRECTORIES(${GTEST_SOURCE_DIR})

# mfh 03 Nov 2017: The gtest library used here must have a different
# name than that of the gtest library built in KokkosCore.  We can't
# just refer to the library in KokkosCore's tests, because it's
# possible to build only (e.g.,) KokkosAlgorithms tests, without
# building KokkosCore tests.


KOKKOS_ADD_TEST_LIBRARY(
  kokkosalgorithms_gtest
  HEADERS ${GTEST_SOURCE_DIR}/gtest/gtest.h
  SOURCES ${GTEST_SOURCE_DIR}/gtest/gtest-all.cc
)
# WORKAROUND FOR HIPCC
IF(Kokkos_ENABLE_HIP)
  TARGET_COMPILE_DEFINITIONS(kokkosalgorithms_gtest PUBLIC "-DGTEST_HAS_PTHREAD=0 --amdgpu-target=gfx906")
ELSE()
  TARGET_COMPILE_DEFINITIONS(kokkosalgorithms_gtest PUBLIC "-DGTEST_HAS_PTHREAD=0")
ENDIF()

TARGET_COMPILE_FEATURES(kokkosalgorithms_gtest PUBLIC cxx_std_11)

SET(SOURCES
  UnitTestMain.cpp
)

IF(Kokkos_ENABLE_OPENMP)
  LIST( APPEND SOURCES
    TestOpenMP.cpp
    TestOpenMP_Sort1D.cpp
    TestOpenMP_Sort3D.cpp
    TestOpenMP_SortDynamicView.cpp
    TestOpenMP_Random.cpp
  )
ENDIF()

IF(Kokkos_ENABLE_HIP)
  LIST( APPEND SOURCES
    TestHIP.cpp
  )
ENDIF()

IF(Kokkos_ENABLE_CUDA)
  LIST( APPEND SOURCES
    TestCuda.cpp
  )
ENDIF()

IF(Kokkos_ENABLE_HPX)
  LIST( APPEND SOURCES
    TestHPX.cpp
  )
ENDIF()

IF(Kokkos_ENABLE_SERIAL)
  LIST( APPEND SOURCES
    TestSerial.cpp
  )
ENDIF()

IF(Kokkos_ENABLE_PTHREAD)
  LIST( APPEND SOURCES
    TestThreads.cpp
  )
ENDIF()

KOKKOS_ADD_EXECUTABLE_AND_TEST(
  UnitTest
  SOURCES ${SOURCES}
)
