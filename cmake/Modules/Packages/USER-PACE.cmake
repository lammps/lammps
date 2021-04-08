
set(PACELIB_URL "https://github.com/ICAMS/lammps-user-pace/archive/refs/tags/v.2021.2.3.upd2.tar.gz" CACHE STRING "URL for PACE evaluator library sources")
set(PACELIB_MD5 "8fd1162724d349b930e474927197f20d" CACHE STRING "MD5 checksum of PACE evaluator library tarball")
mark_as_advanced(PACELIB_URL)
mark_as_advanced(PACELIB_MD5)

# download library sources to build folder
file(DOWNLOAD ${PACELIB_URL} ./libpace.tar.gz SHOW_PROGRESS EXPECTED_HASH MD5=${PACELIB_MD5})

# uncompress downloaded sources
execute_process(
  COMMAND ${CMAKE_COMMAND} -E remove_directory lammps-user-pace*
  COMMAND ${CMAKE_COMMAND} -E tar xzf libpace.tar.gz
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)


file(GLOB PACE_EVALUATOR_INCLUDE_DIR ${CMAKE_BINARY_DIR}/lammps-user-pace-*/USER-PACE)
file(GLOB PACE_EVALUATOR_SOURCES ${CMAKE_BINARY_DIR}/lammps-user-pace-*/USER-PACE/*.cpp)
list(FILTER PACE_EVALUATOR_SOURCES EXCLUDE REGEX pair_pace.cpp)

add_library(pace STATIC ${PACE_EVALUATOR_SOURCES})
set_target_properties(pace PROPERTIES OUTPUT_NAME pace${LAMMPS_MACHINE})
target_include_directories(pace PUBLIC ${PACE_EVALUATOR_INCLUDE_DIR})
target_link_libraries(lammps PRIVATE pace)

