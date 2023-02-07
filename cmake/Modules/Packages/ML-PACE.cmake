set(PACELIB_URL "https://github.com/ICAMS/lammps-user-pace/archive/refs/tags/v.2021.10.25.fix2.tar.gz" CACHE STRING "URL for PACE evaluator library sources")

set(PACELIB_MD5 "32394d799bc282bb57696c78c456e64f" CACHE STRING "MD5 checksum of PACE evaluator library tarball")
mark_as_advanced(PACELIB_URL)
mark_as_advanced(PACELIB_MD5)
GetFallbackURL(PACELIB_URL PACELIB_FALLBACK)

# download library sources to build folder
if(EXISTS ${CMAKE_BINARY_DIR}/libpace.tar.gz)
  file(MD5 ${CMAKE_BINARY_DIR}/libpace.tar.gz DL_MD5)
endif()
if(NOT "${DL_MD5}" STREQUAL "${PACELIB_MD5}")
  message(STATUS "Downloading ${PACELIB_URL}")
  file(DOWNLOAD ${PACELIB_URL} ${CMAKE_BINARY_DIR}/libpace.tar.gz STATUS DL_STATUS SHOW_PROGRESS)
  file(MD5 ${CMAKE_BINARY_DIR}/libpace.tar.gz DL_MD5)
  if((NOT DL_STATUS EQUAL 0) OR (NOT "${DL_MD5}" STREQUAL "${PACELIB_MD5}"))
    message(WARNING "Download from primary URL ${PACELIB_URL} failed\nTrying fallback URL ${PACELIB_FALLBACK}")
    file(DOWNLOAD ${PACELIB_FALLBACK} ${CMAKE_BINARY_DIR}/libpace.tar.gz EXPECTED_HASH MD5=${PACELIB_MD5} SHOW_PROGRESS)
  endif()
else()
  message(STATUS "Using already downloaded archive ${CMAKE_BINARY_DIR}/libpace.tar.gz")
endif()

# uncompress downloaded sources
execute_process(
  COMMAND ${CMAKE_COMMAND} -E remove_directory lammps-user-pace*
  COMMAND ${CMAKE_COMMAND} -E tar xzf libpace.tar.gz
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
get_newest_file(${CMAKE_BINARY_DIR}/lammps-user-pace-* lib-pace)

# enforce building libyaml-cpp as static library and turn off optional features
set(YAML_BUILD_SHARED_LIBS OFF)
set(YAML_CPP_BUILD_CONTRIB OFF)
set(YAML_CPP_BUILD_TOOLS OFF)
add_subdirectory(${lib-pace}/yaml-cpp build-yaml-cpp)
set(YAML_CPP_INCLUDE_DIR ${lib-pace}/yaml-cpp/include)

file(GLOB PACE_EVALUATOR_INCLUDE_DIR ${CONFIGURE_DEPENDS} ${lib-pace}/ML-PACE)
file(GLOB PACE_EVALUATOR_SOURCES ${CONFIGURE_DEPENDS} ${lib-pace}/ML-PACE/*.cpp)
list(FILTER PACE_EVALUATOR_SOURCES EXCLUDE REGEX pair_pace.cpp)

add_library(pace STATIC ${PACE_EVALUATOR_SOURCES})
set_target_properties(pace PROPERTIES CXX_EXTENSIONS ON OUTPUT_NAME lammps_pace${LAMMPS_MACHINE})
target_include_directories(pace PUBLIC ${PACE_EVALUATOR_INCLUDE_DIR} ${YAML_CPP_INCLUDE_DIR})


target_link_libraries(pace PRIVATE yaml-cpp-pace)
if(CMAKE_PROJECT_NAME STREQUAL "lammps")
  target_link_libraries(lammps PRIVATE pace)
endif()
