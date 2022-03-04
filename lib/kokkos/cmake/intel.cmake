
FUNCTION(kokkos_set_intel_flags full_standard int_standard)
  STRING(TOLOWER ${full_standard} FULL_LC_STANDARD)
  STRING(TOLOWER ${int_standard} INT_LC_STANDARD)
  # The following three blocks of code were copied from
  # /Modules/Compiler/Intel-CXX.cmake from CMake 3.18.1 and then modified.
  IF(CMAKE_CXX_SIMULATE_ID STREQUAL MSVC)
    SET(_std -Qstd)
    SET(_ext c++)
  ELSE()
    SET(_std -std)
    SET(_ext gnu++)
  ENDIF()
  SET(KOKKOS_CXX_STANDARD_FLAG             "${_std}=c++${FULL_LC_STANDARD}" PARENT_SCOPE)
  SET(KOKKOS_CXX_INTERMDIATE_STANDARD_FLAG "${_std}=${_ext}${INT_LC_STANDARD}" PARENT_SCOPE)
ENDFUNCTION()


