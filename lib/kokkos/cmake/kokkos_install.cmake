IF (NOT KOKKOS_HAS_TRILINOS)
  INCLUDE(GNUInstallDirs)

  #Set all the variables needed for KokkosConfig.cmake
  GET_PROPERTY(KOKKOS_PROP_LIBS GLOBAL PROPERTY KOKKOS_LIBRARIES_NAMES)
  SET(KOKKOS_LIBRARIES ${KOKKOS_PROP_LIBS})

  INCLUDE(CMakePackageConfigHelpers)
  CONFIGURE_PACKAGE_CONFIG_FILE(
    cmake/KokkosConfig.cmake.in
    "${Kokkos_BINARY_DIR}/KokkosConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR}/cmake)

  INCLUDE(CMakePackageConfigHelpers)
  CONFIGURE_PACKAGE_CONFIG_FILE(
	  cmake/KokkosConfigCommon.cmake.in
	  "${Kokkos_BINARY_DIR}/KokkosConfigCommon.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR}/cmake)

  WRITE_BASIC_PACKAGE_VERSION_FILE("${Kokkos_BINARY_DIR}/KokkosConfigVersion.cmake"
      VERSION "${Kokkos_VERSION}"
      COMPATIBILITY SameMajorVersion)

  # Install the KokkosConfig*.cmake files
  install(FILES
    "${Kokkos_BINARY_DIR}/KokkosConfig.cmake"
    "${Kokkos_BINARY_DIR}/KokkosConfigCommon.cmake"
    "${Kokkos_BINARY_DIR}/KokkosConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Kokkos)
  install(EXPORT KokkosTargets NAMESPACE Kokkos:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Kokkos)
ELSE()
  CONFIGURE_FILE(cmake/KokkosConfigCommon.cmake.in ${Kokkos_BINARY_DIR}/KokkosConfigCommon.cmake @ONLY)
  file(READ ${Kokkos_BINARY_DIR}/KokkosConfigCommon.cmake KOKKOS_CONFIG_COMMON)
  file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/KokkosConfig_install.cmake" "${KOKKOS_CONFIG_COMMON}")
  CONFIGURE_FILE(cmake/KokkosTrilinosConfig.cmake.in ${Kokkos_BINARY_DIR}/KokkosTrilinosConfig.cmake @ONLY)
  file(READ ${Kokkos_BINARY_DIR}/KokkosTrilinosConfig.cmake KOKKOS_TRILINOS_CONFIG)
  file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/KokkosConfig_install.cmake" "${KOKKOS_TRILINOS_CONFIG}")
ENDIF()

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/KokkosCore_config.h DESTINATION ${KOKKOS_HEADER_DIR})

