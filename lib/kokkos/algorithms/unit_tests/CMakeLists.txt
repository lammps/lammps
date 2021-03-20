
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

# avoid deprecation warnings from MSVC
TARGET_COMPILE_DEFINITIONS(kokkosalgorithms_gtest PUBLIC GTEST_HAS_TR1_TUPLE=0 GTEST_HAS_PTHREAD=0)

IF((NOT (Kokkos_ENABLE_CUDA AND WIN32)) AND (NOT ("${KOKKOS_CXX_COMPILER_ID}" STREQUAL "Fujitsu")))
TARGET_COMPILE_FEATURES(kokkosalgorithms_gtest PUBLIC cxx_std_11)
ENDIF()

# Suppress clang-tidy diagnostics on code that we do not have control over
IF(CMAKE_CXX_CLANG_TIDY)
  SET_TARGET_PROPERTIES(kokkosalgorithms_gtest PROPERTIES CXX_CLANG_TIDY "")
ENDIF()

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
