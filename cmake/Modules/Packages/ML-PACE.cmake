set(PACELIB_URL "https://github.com/ICAMS/lammps-user-pace/archive/refs/tags/v.2021.10.25.fix2.tar.gz" CACHE STRING "URL for PACE evaluator library sources")

set(PACELIB_MD5 "d1984c9d5bf88715fc20b23dcf95aa12" CACHE STRING "MD5 checksum of PACE evaluator library tarball")
mark_as_advanced(PACELIB_URL)
mark_as_advanced(PACELIB_MD5)

# download library sources to build folder
file(DOWNLOAD ${PACELIB_URL} ${CMAKE_BINARY_DIR}/libpace.tar.gz EXPECTED_HASH MD5=${PACELIB_MD5}) #SHOW_PROGRESS

# uncompress downloaded sources
execute_process(
  COMMAND ${CMAKE_COMMAND} -E remove_directory lammps-user-pace*
  COMMAND ${CMAKE_COMMAND} -E tar xzf libpace.tar.gz
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

file(GLOB lib-pace ${CMAKE_BINARY_DIR}/lammps-user-pace-*)
add_subdirectory(${lib-pace}/yaml-cpp build-yaml-cpp)
set(YAML_CPP_INCLUDE_DIR ${lib-pace}/yaml-cpp/include)

file(GLOB PACE_EVALUATOR_INCLUDE_DIR ${lib-pace}/ML-PACE)
file(GLOB PACE_EVALUATOR_SOURCES ${lib-pace}/ML-PACE/*.cpp)
list(FILTER PACE_EVALUATOR_SOURCES EXCLUDE REGEX pair_pace.cpp)

add_library(pace STATIC ${PACE_EVALUATOR_SOURCES})
set_target_properties(pace PROPERTIES CXX_EXTENSIONS ON OUTPUT_NAME lammps_pace${LAMMPS_MACHINE})
target_include_directories(pace PUBLIC ${PACE_EVALUATOR_INCLUDE_DIR} ${YAML_CPP_INCLUDE_DIR})


target_link_libraries(pace PRIVATE yaml-cpp-pace)

target_link_libraries(lammps PRIVATE pace)
