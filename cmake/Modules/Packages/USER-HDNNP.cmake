find_package(N2P2 QUIET)
if(N2P2_FOUND)
  set(DOWNLOAD_N2P2_DEFAULT OFF)
else()
  set(DOWNLOAD_N2P2_DEFAULT ON)
endif()
option(DOWNLOAD_N2P2 "Download n2p2 library instead of using an already installed one)" ${DOWNLOAD_N2P2_DEFAULT})
if(DOWNLOAD_N2P2)
  set(N2P2_URL "https://github.com/CompPhysVienna/n2p2/archive/v2.1.2.tar.gz" CACHE STRING "URL for n2p2 tarball")
  set(N2P2_MD5 "20cf194d14b1f1c72f38879bafda67e2" CACHE STRING "MD5 checksum of N2P2 tarball")
  mark_as_advanced(N2P2_URL)
  mark_as_advanced(N2P2_MD5)

  include(ExternalProject)
  ExternalProject_Add(n2p2_build
    DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}
    URL     ${N2P2_URL}
    URL_MD5 ${N2P2_MD5}
    UPDATE_COMMAND ""
    SOURCE_SUBDIR src/
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND "make libnnpif"
    INSTALL_COMMAND ""
    #BUILD_BYPRODUCTS <INSTALL_DIR>/lib/libnnp.a <INSTALL_DIR>/lib/libnnpif.a
    )
  ExternalProject_get_property(n2p2_build INSTALL_DIR)
  add_library(LAMMPS::N2P2 STATIC IMPORTED GLOBAL)
  set_target_properties(LAMMPS::N2P2 PROPERTIES
    IMPORTED_LOCATION "${INSTALL_DIR}/lib/libnnp.a"
    INTERFACE_INCLUDE_DIRECTORIES "${INSTALL_DIR}/include")
  target_link_libraries(lammps PRIVATE LAMMPS::N2P2)
  add_dependencies(LAMMPS::N2P2 n2p2_build)
else()
  find_package(N2P2)
  if(NOT N2P2_FOUND)
    message(FATAL_ERROR "n2p2 not found, help CMake to find it by setting N2P2_DIR, or set DOWNLOAD_N2P2=ON to download it")
  endif()
  target_link_libraries(lammps PRIVATE N2P2::N2P2)
  include(${N2P2_CMAKE_EXTRAS})
endif()
