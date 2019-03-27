

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------------

SET(TRILINOS_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

#-----------------------------------------------------------------------------

IF(KOKKOS_LEGACY_TRIBITS)

  MESSAGE("LEGACY STUFF GETTING CALLED")

  IF(KOKKOS_ENABLE_EXPLICIT_INSTANTIATION)
    MESSAGE("GOING INTO ETI DIR")
    ADD_SUBDIRECTORY(eti)
    INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/eti")
  ENDIF()

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_CXX11)
  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_CUDA)

  SET(HEADERS_PUBLIC "")
  SET(HEADERS_PRIVATE "")
  SET(SOURCES "")

  FILE(GLOB HEADERS_PUBLIC Kokkos*.hpp)
  LIST( APPEND HEADERS_PUBLIC ${CMAKE_BINARY_DIR}/${PACKAGE_NAME}_config.h )

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

  IF(KOKKOS_ENABLE_EXPLICIT_INSTANTIATION)
    LIST(APPEND SOURCES       ${ETI_SOURCES} )
  ENDIF()

  INSTALL(FILES ${HEADERS_QTHREADS} DESTINATION ${TRILINOS_INCDIR}/Qthreads/)

  TRIBITS_ADD_LIBRARY(
      kokkoscore
      HEADERS ${HEADERS_PUBLIC}
      NOINSTALLHEADERS ${HEADERS_PRIVATE}
      SOURCES ${SOURCES}
      DEPLIBS
      )

#-----------------------------------------------------------------------------
#  In the new build system, sources are calculated by Makefile.kokkos
else()

  INSTALL (DIRECTORY
           "${CMAKE_CURRENT_SOURCE_DIR}/"
           DESTINATION ${TRILINOS_INCDIR}
           FILES_MATCHING PATTERN "*.hpp"
  )

  TRIBITS_ADD_LIBRARY(
      kokkoscore
      SOURCES ${KOKKOS_CORE_SRCS}
      DEPLIBS
      )

endif()
#-----------------------------------------------------------------------------

# build and install pkgconfig file
CONFIGURE_FILE(kokkos.pc.in kokkos.pc @ONLY)
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/kokkos.pc DESTINATION lib/pkgconfig)
