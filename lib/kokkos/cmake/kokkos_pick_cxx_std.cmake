# From CMake 3.10 documentation

#This can run at any time
KOKKOS_OPTION(CXX_STANDARD "" STRING "[[DEPRECATED - USE CMAKE_CXX_STANDARD INSTEAD]] The C++ standard for Kokkos to use: 17 or 20. If empty, this will default to CMAKE_CXX_STANDARD. If both CMAKE_CXX_STANDARD and Kokkos_CXX_STANDARD are empty, this will default to 17")

# Set CXX standard flags
SET(KOKKOS_ENABLE_CXX17 OFF)
SET(KOKKOS_ENABLE_CXX20 OFF)
SET(KOKKOS_ENABLE_CXX23 OFF)
IF (KOKKOS_CXX_STANDARD)
  MESSAGE(FATAL_ERROR "Setting the variable Kokkos_CXX_STANDARD in configuration is deprecated - set CMAKE_CXX_STANDARD directly instead")
ENDIF()

IF (NOT CMAKE_CXX_STANDARD)
  SET(KOKKOS_CXX_STANDARD "17")
ELSE()
  SET(KOKKOS_CXX_STANDARD ${CMAKE_CXX_STANDARD})
ENDIF()
MESSAGE(STATUS "Setting default Kokkos CXX standard to ${KOKKOS_CXX_STANDARD}")
