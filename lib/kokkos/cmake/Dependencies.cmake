IF (CMAKE_CXX_STANDARD GREATER_EQUAL 17)
  SET(KOKKOS_SIMD_TEST_CLASS PT)
ELSE()
  SET(KOKKOS_SIMD_TEST_CLASS EX)
  IF (${PROJECT_NAME}_ENABLE_KokkosSimd)
    MESSAGE(WARNING "KokkosSimd is explicitly enabled but C++17 is not available")
  ELSE()
    MESSAGE(STATUS "Disabling KokkosSimd by default because C++17 is not available")
  ENDIF()
ENDIF()

TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
    #SubPackageName       Directory         Class    Req/Opt
    #
    # New Kokkos subpackages:
    Core                  core              PS       REQUIRED
    Containers            containers        PS       OPTIONAL
    Algorithms            algorithms        PS       OPTIONAL
    Simd                  simd              ${KOKKOS_SIMD_TEST_CLASS}       OPTIONAL
  )
