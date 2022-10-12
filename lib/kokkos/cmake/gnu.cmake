
FUNCTION(kokkos_set_gnu_flags full_standard int_standard)
  STRING(TOLOWER ${full_standard} FULL_LC_STANDARD)
  STRING(TOLOWER ${int_standard} INT_LC_STANDARD)
  # The following three blocks of code were copied from
  # /Modules/Compiler/Intel-CXX.cmake from CMake 3.7.2 and then modified.
  IF(CMAKE_CXX_SIMULATE_ID STREQUAL MSVC)
    SET(_std -Qstd)
    SET(_ext c++)
  ELSE()
    SET(_std -std)
    SET(_ext gnu++)
  ENDIF()

  IF (CMAKE_CXX_EXTENSIONS)
    SET(KOKKOS_CXX_STANDARD_FLAG "-std=gnu++${FULL_LC_STANDARD}" PARENT_SCOPE)
    SET(KOKKOS_CXX_INTERMEDIATE_STANDARD_FLAG "-std=gnu++${INT_LC_STANDARD}" PARENT_SCOPE)
  ELSE()
    SET(KOKKOS_CXX_STANDARD_FLAG "-std=c++${FULL_LC_STANDARD}" PARENT_SCOPE)
    SET(KOKKOS_CXX_INTERMEDIATE_STANDARD_FLAG "-std=c++${INT_LC_STANDARD}" PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()

