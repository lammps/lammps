set(PACELIB_URL "https://github.com/ICAMS/lammps-user-pace/archive/refs/tags/v.2023.01.3.tar.gz" CACHE STRING "URL for PACE evaluator library sources")

set(PACELIB_MD5 "f418d32b60e531063ac4285bf702b468" CACHE STRING "MD5 checksum of PACE evaluator library tarball")
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
get_newest_file(${CMAKE_BINARY_DIR}/lammps-user-pace-* lib-pace)

add_subdirectory(${lib-pace} build-pace)
set_target_properties(pace PROPERTIES CXX_EXTENSIONS ON OUTPUT_NAME lammps_pace${LAMMPS_MACHINE})

if(CMAKE_PROJECT_NAME STREQUAL "lammps")
  target_link_libraries(lammps PRIVATE pace)
endif()
