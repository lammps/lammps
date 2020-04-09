find_package(VORO)
if(VORO_FOUND)
  set(DOWNLOAD_VORO_DEFAULT OFF)
else()
  set(DOWNLOAD_VORO_DEFAULT ON)
endif()
option(DOWNLOAD_VORO "Download and compile the Voro++ library instead of using an already installed one" ${DOWNLOAD_VORO_DEFAULT})
if(DOWNLOAD_VORO)
  message(STATUS "Voro++ download requested - we will build our own")
  include(ExternalProject)

  if(BUILD_SHARED_LIBS)
    set(VORO_BUILD_CFLAGS "${CMAKE_SHARED_LIBRARY_CXX_FLAGS} ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${BTYPE}}")
  else()
    set(VORO_BUILD_CFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${BTYPE}}")
  endif()
  if(APPLE)
    get_filename_component(VORO_CXX ${CMAKE_CXX_COMPILER} NAME_WE)
    set(VORO_BUILD_OPTIONS CXX=${VORO_CXX} CFLAGS=${VORO_BUILD_CFLAGS})
  else()
    set(VORO_BUILD_OPTIONS CXX=${CMAKE_CXX_COMPILER} CFLAGS=${VORO_BUILD_CFLAGS})
  endif()

  ExternalProject_Add(voro_build
    URL https://download.lammps.org/thirdparty/voro++-0.4.6.tar.gz
    URL_MD5 2338b824c3b7b25590e18e8df5d68af9
    CONFIGURE_COMMAND "" BUILD_COMMAND make ${VORO_BUILD_OPTIONS} BUILD_IN_SOURCE 1 INSTALL_COMMAND ""
    BUILD_BYPRODUCTS <SOURCE_DIR>/src/libvoro++.a
    )
  ExternalProject_get_property(voro_build SOURCE_DIR)
  file(MAKE_DIRECTORY ${SOURCE_DIR}/src)
  add_library(LAMMPS::VORO UNKNOWN IMPORTED)
  set_target_properties(LAMMPS::VORO PROPERTIES
    IMPORTED_LOCATION "${SOURCE_DIR}/src/libvoro++.a"
    INTERFACE_INCLUDE_DIRECTORIES "${SOURCE_DIR}/src")
  target_link_libraries(lammps PRIVATE LAMMPS::VORO)
  add_dependencies(LAMMPS::VORO voro_build)
  if(BUILD_LIB)
    install(CODE "MESSAGE(FATAL_ERROR \"Installing liblammps with downloaded libraries is currently not supported.\")")
  endif()
else()
  find_package(VORO)
  if(NOT VORO_FOUND)
    message(FATAL_ERROR "Voro++ library not found. Help CMake to find it by setting VORO_LIBRARY and VORO_INCLUDE_DIR, or set DOWNLOAD_VORO=ON to download it")
  endif()
  target_link_libraries(lammps PRIVATE VORO::VORO)
endif()
