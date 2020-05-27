# From CMake 3.10 documentation

#This can run at any time
KOKKOS_OPTION(CXX_STANDARD "" STRING "The C++ standard for Kokkos to use: 11, 14, 17, or 20. If empty, this will default to CMAKE_CXX_STANDARD. If both CMAKE_CXX_STANDARD and Kokkos_CXX_STANDARD are empty, this will default to 11")

# Set CXX standard flags
SET(KOKKOS_ENABLE_CXX11 OFF)
SET(KOKKOS_ENABLE_CXX14 OFF)
SET(KOKKOS_ENABLE_CXX17 OFF)
SET(KOKKOS_ENABLE_CXX20 OFF)
IF (KOKKOS_CXX_STANDARD)
  IF (${KOKKOS_CXX_STANDARD} STREQUAL "c++98")
    MESSAGE(FATAL_ERROR "Kokkos no longer supports C++98 - minimum C++11")
  ELSEIF (${KOKKOS_CXX_STANDARD} STREQUAL "c++11")
    MESSAGE(WARNING "Deprecated Kokkos C++ standard set as 'c++11'. Use '11' instead.")
    SET(KOKKOS_CXX_STANDARD "11")
  ELSEIF(${KOKKOS_CXX_STANDARD} STREQUAL "c++14")
    MESSAGE(WARNING "Deprecated Kokkos C++ standard set as 'c++14'. Use '14' instead.")
    SET(KOKKOS_CXX_STANDARD "14")
  ELSEIF(${KOKKOS_CXX_STANDARD} STREQUAL "c++17")
    MESSAGE(WARNING "Deprecated Kokkos C++ standard set as 'c++17'. Use '17' instead.")
    SET(KOKKOS_CXX_STANDARD "17")
  ELSEIF(${KOKKOS_CXX_STANDARD} STREQUAL "c++1y")
    MESSAGE(WARNING "Deprecated Kokkos C++ standard set as 'c++1y'. Use '1Y' instead.")
    SET(KOKKOS_CXX_STANDARD "1Y")
  ELSEIF(${KOKKOS_CXX_STANDARD} STREQUAL "c++1z")
    MESSAGE(WARNING "Deprecated Kokkos C++ standard set as 'c++1z'. Use '1Z' instead.")
    SET(KOKKOS_CXX_STANDARD "1Z")
  ELSEIF(${KOKKOS_CXX_STANDARD} STREQUAL "c++2a")
    MESSAGE(WARNING "Deprecated Kokkos C++ standard set as 'c++2a'. Use '2A' instead.")
    SET(KOKKOS_CXX_STANDARD "2A")
  ENDIF()
ENDIF()

IF (NOT KOKKOS_CXX_STANDARD AND NOT CMAKE_CXX_STANDARD)
  MESSAGE(STATUS "Setting default Kokkos CXX standard to 11")
  SET(KOKKOS_CXX_STANDARD "11")
ELSEIF(NOT KOKKOS_CXX_STANDARD)
  MESSAGE(STATUS "Setting default Kokkos CXX standard to ${CMAKE_CXX_STANDARD}")
  SET(KOKKOS_CXX_STANDARD ${CMAKE_CXX_STANDARD})
ENDIF()





