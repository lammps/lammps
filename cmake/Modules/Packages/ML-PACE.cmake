# PACE library support for ML-PACE package
find_package(pace QUIET)

if(pace_FOUND)
    find_package(pace)
    target_link_libraries(lammps PRIVATE pace::pace)
else()
    # set policy to silence warnings about timestamps of downloaded files. review occasionally if it may be set to NEW
    if(POLICY CMP0135)
      cmake_policy(SET CMP0135 OLD)
    endif()

    set(PACELIB_URL "https://github.com/ICAMS/lammps-user-pace/archive/refs/tags/v.2023.11.25.fix2.tar.gz" CACHE STRING "URL for PACE evaluator library sources")
    set(PACELIB_MD5 "a53bd87cfee8b07d9f44bc17aad69c3f" CACHE STRING "MD5 checksum of PACE evaluator library tarball")
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

    # some preinstalled yaml-cpp versions don't provide a namespaced target
    find_package(yaml-cpp QUIET)
    if(TARGET yaml-cpp AND NOT TARGET yaml-cpp::yaml-cpp)
      add_library(yaml-cpp::yaml-cpp ALIAS yaml-cpp)
    endif()

    add_subdirectory(${lib-pace} build-pace)
    set_target_properties(pace PROPERTIES CXX_EXTENSIONS ON OUTPUT_NAME lammps_pace${LAMMPS_MACHINE})

    if(CMAKE_PROJECT_NAME STREQUAL "lammps")
      target_link_libraries(lammps PRIVATE pace)
    endif()
endif()
