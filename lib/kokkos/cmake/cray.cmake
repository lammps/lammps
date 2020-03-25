

function(kokkos_set_cray_flags full_standard int_standard)
  STRING(TOLOWER ${full_standard} FULL_LC_STANDARD)
  STRING(TOLOWER ${int_standard} INT_LC_STANDARD)
  SET(KOKKOS_CXX_STANDARD_FLAG "-hstd=c++${FULL_LC_STANDARD}", PARENT_SCOPE)
  SET(KOKKOS_CXX_INTERMDIATE_STANDARD_FLAG "-hstd=c++${INT_LC_STANDARD}" PARENT_SCOPE)
endfunction()

