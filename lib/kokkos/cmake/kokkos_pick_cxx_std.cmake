# From CMake 3.10 documentation

#This can run at any time
kokkos_option(
  CXX_STANDARD
  ""
  STRING
  "[[DEPRECATED - USE CMAKE_CXX_STANDARD INSTEAD]] The C++ standard for Kokkos to use: 17 or 20. If empty, this will default to CMAKE_CXX_STANDARD. If both CMAKE_CXX_STANDARD and Kokkos_CXX_STANDARD are empty, this will default to 17"
)

# Set CXX standard flags
set(KOKKOS_ENABLE_CXX17 OFF)
set(KOKKOS_ENABLE_CXX20 OFF)
set(KOKKOS_ENABLE_CXX23 OFF)
set(KOKKOS_ENABLE_CXX26 OFF)
if(KOKKOS_CXX_STANDARD)
  message(
    FATAL_ERROR
      "Setting the variable Kokkos_CXX_STANDARD in configuration is deprecated - set CMAKE_CXX_STANDARD directly instead"
  )
endif()

if(NOT CMAKE_CXX_STANDARD)
  set(KOKKOS_CXX_STANDARD "17")
else()
  set(KOKKOS_CXX_STANDARD ${CMAKE_CXX_STANDARD})
endif()
message(STATUS "Setting default Kokkos CXX standard to ${KOKKOS_CXX_STANDARD}")
