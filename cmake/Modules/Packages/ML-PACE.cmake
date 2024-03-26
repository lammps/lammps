set(PACELIB_URL "https://github.com/ICAMS/lammps-user-pace/archive/refs/tags/v.2023.11.25.fix.tar.gz" CACHE STRING "URL for PACE evaluator library sources")

set(PACELIB_MD5 "b45de9a633f42ed65422567e3ce56f9f" CACHE STRING "MD5 checksum of PACE evaluator library tarball")
mark_as_advanced(PACELIB_URL)
mark_as_advanced(PACELIB_MD5)
GetFallbackURL(PACELIB_URL PACELIB_FALLBACK)

# LOCAL_ML-PACE points to top-level dir with local lammps-user-pace repo,
# to make it easier to check local build without going through the public github releases
if(LOCAL_ML-PACE)
 set(lib-pace "${LOCAL_ML-PACE}")
else()
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
endif()

add_subdirectory(${lib-pace} build-pace)
set_target_properties(pace PROPERTIES CXX_EXTENSIONS ON OUTPUT_NAME lammps_pace${LAMMPS_MACHINE})

if(CMAKE_PROJECT_NAME STREQUAL "lammps")
  target_link_libraries(lammps PRIVATE pace)
endif()
