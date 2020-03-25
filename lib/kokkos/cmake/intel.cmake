
FUNCTION(kokkos_set_intel_flags full_standard int_standard)
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

  IF(NOT KOKKOS_CXX_STANDARD STREQUAL 11 AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15.0.2)
    #There is no gnu++14 value supported; figure out what to do.
    SET(KOKKOS_CXX_STANDARD_FLAG "${_std}=c++${FULL_LC_STANDARD}" PARENT_SCOPE)
    SET(KOKKOS_CXX_INTERMEDIATE_STANDARD_FLAG "${_std}=c++${INT_LC_STANDARD}" PARENT_SCOPE)
  ELSEIF(KOKKOS_CXX_STANDARD STREQUAL 11 AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13.0)
    IF (CMAKE_CXX_EXTENSIONS)
      SET(KOKKOS_CXX_STANDARD_FLAG "${_std}=${_ext}c++11" PARENT_SCOPE)
    ELSE()
      SET(KOKKOS_CXX_STANDARD_FLAG "${_std}=c++11" PARENT_SCOPE)
    ENDIF()
  ELSE()
    MESSAGE(FATAL_ERROR "Intel compiler version too low - need 13.0 for C++11 and 15.0 for C++14")
  ENDIF()

ENDFUNCTION()

