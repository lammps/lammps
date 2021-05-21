find_package(N2P2 QUIET)
if(N2P2_FOUND)
  set(DOWNLOAD_N2P2_DEFAULT OFF)
else()
  set(DOWNLOAD_N2P2_DEFAULT ON)
endif()
option(DOWNLOAD_N2P2 "Download n2p2 library instead of using an already installed one)" ${DOWNLOAD_N2P2_DEFAULT})
if(DOWNLOAD_N2P2)
  set(N2P2_URL "https://github.com/CompPhysVienna/n2p2/archive/v2.1.4.tar.gz" CACHE STRING "URL for n2p2 tarball")
  set(N2P2_MD5 "9595b066636cd6b90b0fef93398297a5" CACHE STRING "MD5 checksum of N2P2 tarball")
  mark_as_advanced(N2P2_URL)
  mark_as_advanced(N2P2_MD5)

  if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(N2P2_COMP llvm)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(N2P2_COMP intel)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(N2P2_COMP gnu)
  else() # default
    set(N2P2_COMP gnu)
  endif()

  if(NOT BUILD_MPI)
    set(N2P2_PROJECT_OPTIONS "-DN2P2_NO_MPI")
  endif()

  string(TOUPPER "${CMAKE_BUILD_TYPE}" BTYPE)
  set(N2P2_BUILD_FLAGS "${CMAKE_SHARED_LIBRARY_CXX_FLAGS} ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${BTYPE}}")
  set(N2P2_BUILD_OPTIONS INTERFACES=LAMMPS COMP=${N2P2_COMP} "PROJECT_OPTIONS=${N2P2_PROJECT_OPTIONS}"
    "PROJECT_CC=${CMAKE_CXX_COMPILER}" "PROJECT_MPICC=${MPI_CXX_COMPILER}" "PROJECT_CFLAGS=${N2P2_BUILD_FLAGS}")
  message(STATUS "N2P2 BUILD OPTIONS: ${N2P2_BUILD_OPTIONS}")

  include(ExternalProject)
  ExternalProject_Add(n2p2_build
    URL     ${N2P2_URL}
    URL_MD5 ${N2P2_MD5}
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make -f makefile libnnpif ${N2P2_BUILD_OPTIONS}
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    LOG_BUILD ON
    SOURCE_SUBDIR src/
    #BUILD_BYPRODUCTS <INSTALL_DIR>/lib/libnnp.a <INSTALL_DIR>/lib/libnnpif.a
    )
  ExternalProject_get_property(n2p2_build SOURCE_DIR)
  # n2p2 core library "libnnp"
  add_library(LAMMPS::N2P2::LIBNNP UNKNOWN IMPORTED)
  set_target_properties(LAMMPS::N2P2::LIBNNP PROPERTIES
    IMPORTED_LOCATION "${SOURCE_DIR}/lib/libnnp.a"
    INTERFACE_INCLUDE_DIRECTORIES "${SOURCE_DIR}/include")
  # n2p2 interface library "libnnpif"
  add_library(LAMMPS::N2P2::LIBNNPIF UNKNOWN IMPORTED)
  set_target_properties(LAMMPS::N2P2::LIBNNPIF PROPERTIES
    IMPORTED_LOCATION "${SOURCE_DIR}/lib/libnnpif.a"
    INTERFACE_INCLUDE_DIRECTORIES "${SOURCE_DIR}/include")
  # Put libnnp, libnnpif and include directory together.
  add_library(LAMMPS::N2P2 INTERFACE IMPORTED)
  set_property(TARGET LAMMPS::N2P2 PROPERTY
    INTERFACE_LINK_LIBRARIES LAMMPS::N2P2::LIBNNPIF LAMMPS::N2P2::LIBNNP)
  target_link_libraries(lammps PRIVATE LAMMPS::N2P2)

  add_dependencies(LAMMPS::N2P2 n2p2_build)
  file(MAKE_DIRECTORY "${SOURCE_DIR}/include")
  file(MAKE_DIRECTORY "${SOURCE_DIR}/lib")
  file(TOUCH "${SOURCE_DIR}/lib/lammps-extra.cmake")
  include("${SOURCE_DIR}/lib/lammps-extra.cmake")
else()
  find_package(N2P2)
  if(NOT N2P2_FOUND)
    message(FATAL_ERROR "n2p2 not found, help CMake to find it by setting N2P2_DIR, or set DOWNLOAD_N2P2=ON to download it")
  endif()
  target_link_libraries(lammps PRIVATE N2P2::N2P2)
  include(${N2P2_CMAKE_EXTRAS})
endif()
