# PACE library support for ML-PACE package
find_package(pace QUIET)

if(pace_FOUND)
    find_package(pace)
    target_link_libraries(lammps PRIVATE pace::pace)
else()
    # set policy to use the time of extraction as timestamps of files unpacked from downloaded
    # archives, so that updating an archive version triggers rebuilding all dependent objects
    if(POLICY CMP0135)
      cmake_policy(SET CMP0135 NEW)
    endif()

    SetDownloadSettings(PACELIB "PACE evaluator library"
      "https://github.com/ICAMS/lammps-user-pace/archive/refs/tags/v.2025.12.4.p1.tar.gz"
      "21e9d7ad2094eef0f19958d154866fc725fc6ccfa82ec3681ef2b006545ced96")
    GetFallbackURL(PACELIB_URL PACELIB_FALLBACK)

    # LOCAL_ML-PACE points to top-level dir with local lammps-user-pace repo,
    # to make it easier to check local build without going through the public github releases
    if(LOCAL_ML-PACE)
     set(lib-pace "${LOCAL_ML-PACE}")
    else()
      # download library sources to build folder
      if(EXISTS ${CMAKE_BINARY_DIR}/libpace.tar.gz)
        file(SHA256 ${CMAKE_BINARY_DIR}/libpace.tar.gz DL_SHA256)
      endif()
      if(NOT "${DL_SHA256}" STREQUAL "${PACELIB_SHA256}")
        message(STATUS "Downloading ${PACELIB_URL}")
        file(DOWNLOAD ${PACELIB_URL} ${CMAKE_BINARY_DIR}/libpace.tar.gz STATUS DL_STATUS SHOW_PROGRESS)
        file(SHA256 ${CMAKE_BINARY_DIR}/libpace.tar.gz DL_SHA256)
        if((NOT DL_STATUS EQUAL 0) OR (NOT "${DL_SHA256}" STREQUAL "${PACELIB_SHA256}"))
          message(WARNING "Download from primary URL ${PACELIB_URL} failed\nTrying fallback URL ${PACELIB_FALLBACK}")
          file(DOWNLOAD ${PACELIB_FALLBACK} ${CMAKE_BINARY_DIR}/libpace.tar.gz EXPECTED_HASH SHA256=${PACELIB_SHA256} SHOW_PROGRESS)
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

    # fixup yaml-cpp/emitterutils.cpp for GCC 15+ until patch is applied
    file(READ ${lib-pace}/yaml-cpp/src/emitterutils.cpp yaml_emitterutils)
    string(REPLACE "#include <sstream>" "#include <sstream>\n#include <cinttypes>" yaml_tmp_emitterutils "${yaml_emitterutils}")
    string(REPLACE "#include <cinttypes>\n#include <cinttypes>" "#include <cinttypes>" yaml_emitterutils "${yaml_tmp_emitterutils}")
    file(WRITE ${lib-pace}/yaml-cpp/src/emitterutils.cpp "${yaml_emitterutils}")

    add_subdirectory(${lib-pace} build-pace EXCLUDE_FROM_ALL)
    set_target_properties(pace PROPERTIES CXX_EXTENSIONS ON OUTPUT_NAME lammps_pace${LAMMPS_MACHINE})

    if(CMAKE_PROJECT_NAME STREQUAL "lammps")
      target_link_libraries(lammps PRIVATE pace)
    endif()
endif()
