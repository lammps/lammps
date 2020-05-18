
KOKKOS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#I have to leave these here for tribits
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------------

FILE(GLOB HEADERS *.hpp)
FILE(GLOB SOURCES *.cpp)
LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h)

#-----------------------------------------------------------------------------

# We have to pass the sources in here for Tribits
# These will get ignored for standalone CMake and a true interface library made
KOKKOS_ADD_INTERFACE_LIBRARY(
  kokkosalgorithms
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
)
KOKKOS_LIB_INCLUDE_DIRECTORIES(kokkosalgorithms
  ${KOKKOS_TOP_BUILD_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)



