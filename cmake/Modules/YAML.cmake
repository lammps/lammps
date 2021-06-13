message(STATUS "Downloading and building YAML library")

include(ExternalProject)
set(YAML_URL "https://pyyaml.org/download/libyaml/yaml-0.2.5.tar.gz" CACHE STRING "URL for libyaml tarball")
set(YAML_MD5 "bb15429d8fb787e7d3f1c83ae129a999" CACHE STRING "MD5 checksum of libyaml tarball")
mark_as_advanced(YAML_URL)
mark_as_advanced(YAML_MD5)

# support cross-compilation to windows
if(CMAKE_CROSSCOMPILING AND (CMAKE_SYSTEM_NAME STREQUAL "Windows"))
  if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86")
    set(YAML_CROSS_HOST --host=i686-mingw64)
  elseif(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
    set(YAML_CROSS_HOST --host=x86_64-mingw64)
  else()
    message(FATAL_ERROR "Unsupported cross-compilation "
      " for ${CMAKE_SYSTEM_NAME}/${CMAKE_SYSTEM_PROCESSOR}"
      " on ${CMAKE_HOST_SYSTEM}/${CMAKE_HOST_SYSTEM_PROCESSOR}")
  endif()
endif()

ExternalProject_Add(libyaml
                    URL               ${YAML_URL}
                    URL_MD5           ${YAML_MD5}
                    SOURCE_DIR        "${CMAKE_BINARY_DIR}/yaml-src"
                    BINARY_DIR        "${CMAKE_BINARY_DIR}/yaml-build"
                    CONFIGURE_COMMAND <SOURCE_DIR>/configure ${CONFIGURE_REQUEST_PIC}
                                      CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER}
                                      --prefix=<INSTALL_DIR> --disable-shared ${YAML_CROSS_HOST}
                    BUILD_BYPRODUCTS  <INSTALL_DIR>/lib/libyaml${CMAKE_STATIC_LIBRARY_SUFFIX}
                    TEST_COMMAND      "")

ExternalProject_Get_Property(libyaml INSTALL_DIR)
set(YAML_INCLUDE_DIR ${INSTALL_DIR}/include)
set(YAML_LIBRARY_DIR ${INSTALL_DIR}/lib)

# workaround for CMake 3.10 on ubuntu 18.04
file(MAKE_DIRECTORY ${YAML_INCLUDE_DIR})
file(MAKE_DIRECTORY ${YAML_LIBRARY_DIR})

set(YAML_LIBRARY_PATH ${INSTALL_DIR}/lib/libyaml${CMAKE_STATIC_LIBRARY_SUFFIX})

add_library(Yaml::Yaml UNKNOWN IMPORTED)
set_target_properties(Yaml::Yaml PROPERTIES
        IMPORTED_LOCATION ${YAML_LIBRARY_PATH}
        INTERFACE_INCLUDE_DIRECTORIES ${YAML_INCLUDE_DIR})
add_dependencies(Yaml::Yaml libyaml)
