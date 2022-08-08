##########################################
# CMake build system for plugin examples.
# The is meant to be used as a template for plugins that are
# distributed independent from the LAMMPS package.
##########################################

cmake_minimum_required(VERSION 3.10)

project(kimplugin VERSION 1.0 LANGUAGES CXX)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
include(CheckIncludeFileCXX)
include(LAMMPSInterfacePlugin)

##########################
# building the plugins

add_library(kimplugin MODULE kimplugin.cpp ${LAMMPS_SOURCE_DIR}/KIM/pair_kim.cpp
  ${LAMMPS_SOURCE_DIR}/KIM/fix_store_kim.cpp ${LAMMPS_SOURCE_DIR}/KIM/kim_command.cpp
  ${LAMMPS_SOURCE_DIR}/KIM/kim_init.cpp ${LAMMPS_SOURCE_DIR}/KIM/kim_interactions.cpp
  ${LAMMPS_SOURCE_DIR}/KIM/kim_param.cpp ${LAMMPS_SOURCE_DIR}/KIM/kim_property.cpp
  ${LAMMPS_SOURCE_DIR}/KIM/kim_query.cpp ${LAMMPS_SOURCE_DIR}/KIM/kim_units.cpp)
target_link_libraries(kimplugin PRIVATE lammps)
target_include_directories(kimplugin PRIVATE ${LAMMPS_SOURCE_DIR}/KIM)
set_target_properties(kimplugin PROPERTIES PREFIX "" SUFFIX ".so")

find_package(KIM-API 2.2.0 CONFIG REQUIRED)
target_link_libraries(kimplugin PRIVATE KIM-API::kim-api)

##########################
# need libcurl
find_package(CURL)
if(CURL_FOUND)
  if(CMAKE_VERSION VERSION_LESS 3.12)
    target_include_directories(kimplugin PRIVATE ${CURL_INCLUDE_DIRS})
    target_link_libraries(kimplugin PRIVATE ${CURL_LIBRARIES})
  else()
    target_link_libraries(kimplugin PRIVATE CURL::libcurl)
  endif()
  target_compile_definitions(kimplugin PRIVATE -DLMP_KIM_CURL)
  set(LMP_DEBUG_CURL OFF CACHE STRING "Set libcurl verbose mode on/off. If on, it displays a lot of verbose information about its operations.")
  mark_as_advanced(LMP_DEBUG_CURL)
  if(LMP_DEBUG_CURL)
    target_compile_definitions(kimplugin PRIVATE -DLMP_DEBUG_CURL)
  endif()
  set(LMP_NO_SSL_CHECK OFF CACHE STRING "Tell libcurl to not verify the peer. If on, the connection succeeds regardless of the names in the certificate. Insecure - Use with caution!")
  mark_as_advanced(LMP_NO_SSL_CHECK)
  if(LMP_NO_SSL_CHECK)
    target_compile_definitions(kimplugin PRIVATE -DLMP_NO_SSL_CHECK)
  endif()
endif()

# MacOS seems to need this
if(CMAKE_SYSTEM_NAME STREQUAL Darwin)
  set_target_properties(kimplugin PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
# tell CMake to export all symbols to a .dll on Windows with special case for MinGW cross-compilers
  set_target_properties(kimplugin PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
  if(CMAKE_CROSSCOMPILING)
    set_target_properties(kimplugin PROPERTIES LINK_FLAGS "-Wl,--export-all-symbols")
  endif()
else()
  set_target_properties(kimplugin PROPERTIES LINK_FLAGS "-rdynamic")
endif()
