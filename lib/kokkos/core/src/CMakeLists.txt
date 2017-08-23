
TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_Serial
  KOKKOS_HAVE_SERIAL
  "Whether to enable the Kokkos::Serial device.  This device executes \"parallel\" kernels sequentially on a single CPU thread.  It is enabled by default.  If you disable this device, please enable at least one other CPU device, such as Kokkos::OpenMP or Kokkos::Threads."
  ON
  )

ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_CXX11)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_CUDA)

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------------

SET(TRILINOS_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

#-----------------------------------------------------------------------------

SET(HEADERS_PUBLIC "")
SET(HEADERS_PRIVATE "")
SET(SOURCES "")

FILE(GLOB HEADERS_PUBLIC Kokkos*.hpp)
LIST( APPEND HEADERS_PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h )

#-----------------------------------------------------------------------------

FILE(GLOB HEADERS_IMPL impl/*.hpp)
FILE(GLOB SOURCES_IMPL impl/*.cpp)

LIST(APPEND HEADERS_PRIVATE ${HEADERS_IMPL} )
LIST(APPEND SOURCES         ${SOURCES_IMPL} )

INSTALL(FILES ${HEADERS_IMPL} DESTINATION ${TRILINOS_INCDIR}/impl/)

#-----------------------------------------------------------------------------

FILE(GLOB HEADERS_THREADS Threads/*.hpp)
FILE(GLOB SOURCES_THREADS Threads/*.cpp)

LIST(APPEND HEADERS_PRIVATE ${HEADERS_THREADS} )
LIST(APPEND SOURCES         ${SOURCES_THREADS} )

INSTALL(FILES ${HEADERS_THREADS} DESTINATION ${TRILINOS_INCDIR}/Threads/)

#-----------------------------------------------------------------------------

FILE(GLOB HEADERS_OPENMP OpenMP/*.hpp)
FILE(GLOB SOURCES_OPENMP OpenMP/*.cpp)

LIST(APPEND HEADERS_PRIVATE ${HEADERS_OPENMP} )
LIST(APPEND SOURCES         ${SOURCES_OPENMP} )

INSTALL(FILES ${HEADERS_OPENMP} DESTINATION ${TRILINOS_INCDIR}/OpenMP/)

#-----------------------------------------------------------------------------

FILE(GLOB HEADERS_CUDA Cuda/*.hpp)
FILE(GLOB SOURCES_CUDA Cuda/*.cpp)

LIST(APPEND HEADERS_PRIVATE ${HEADERS_CUDA} )
LIST(APPEND SOURCES         ${SOURCES_CUDA} )

INSTALL(FILES ${HEADERS_CUDA} DESTINATION ${TRILINOS_INCDIR}/Cuda/)

#-----------------------------------------------------------------------------
FILE(GLOB HEADERS_QTHREADS Qthreads/*.hpp)
FILE(GLOB SOURCES_QTHREADS Qthreads/*.cpp)

LIST(APPEND HEADERS_PRIVATE ${HEADERS_QTHREADS} )
LIST(APPEND SOURCES         ${SOURCES_QTHREADS} )

INSTALL(FILES ${HEADERS_QTHREADS} DESTINATION ${TRILINOS_INCDIR}/Qthreads/)

#-----------------------------------------------------------------------------

TRIBITS_ADD_LIBRARY(
    kokkoscore
    HEADERS ${HEADERS_PUBLIC}
    NOINSTALLHEADERS ${HEADERS_PRIVATE}
    SOURCES ${SOURCES}
    DEPLIBS
    )
