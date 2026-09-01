# Build a CMake based external library as subdirectory.
# The sources will be unpacked to ${CMAKE_BINARY_DIR}/_deps/${target}-src
# The binaries will be built in ${CMAKE_BINARY_DIR}/_deps/${target}-build
#
function(ExternalCMakeProject target url hash basedir cmakedir cmakefile)
  # change settings locally
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)

  get_filename_component(archive ${url} NAME)
  set(source_dir ${CMAKE_BINARY_DIR}/_deps/${target}-src)
  # download and unpack only if this exact archive version has not been unpacked before:
  # the stamp file embeds the SHA256 checksum of the unpacked archive.  unpacking on
  # every CMake run would reset the file timestamps and thus break incremental builds.
  if(NOT EXISTS ${source_dir}/.extracted-${hash})
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/_deps/src)
    if(EXISTS ${CMAKE_BINARY_DIR}/_deps/${archive})
      file(SHA256 ${CMAKE_BINARY_DIR}/_deps/${archive} DL_SHA256)
    endif()
    if(NOT "${DL_SHA256}" STREQUAL "${hash}")
      message(STATUS "Downloading ${url}")
      file(DOWNLOAD ${url} ${CMAKE_BINARY_DIR}/_deps/${archive} STATUS DL_STATUS SHOW_PROGRESS)
      file(SHA256 ${CMAKE_BINARY_DIR}/_deps/${archive} DL_SHA256)
      if((NOT DL_STATUS EQUAL 0) OR (NOT "${DL_SHA256}" STREQUAL "${hash}"))
        set(${target}_URL ${url})
        GetFallbackURL(${target}_URL fallback)
        message(WARNING "Download from primary URL ${url} failed\nTrying fallback URL ${fallback}")
        file(DOWNLOAD ${fallback} ${CMAKE_BINARY_DIR}/_deps/${archive} EXPECTED_HASH SHA256=${hash} SHOW_PROGRESS)
      endif()
    else()
      message(STATUS "Using already downloaded archive ${CMAKE_BINARY_DIR}/_deps/${archive}")
    endif()
    message(STATUS "Unpacking ${archive}")
    # TOUCH sets the timestamps of the unpacked files to the current time instead of
    # keeping the (older) timestamps stored in the archive.  otherwise, after updating
    # to a new archive version, compiled objects in an existing build directory appear
    # newer than the updated headers and are not recompiled, which leads to failed
    # links or subtly inconsistent binaries.
    file(ARCHIVE_EXTRACT INPUT ${CMAKE_BINARY_DIR}/_deps/${archive}
      DESTINATION ${CMAKE_BINARY_DIR}/_deps/src TOUCH)
    file(GLOB TARGET_SOURCE "${CMAKE_BINARY_DIR}/_deps/src/${basedir}*")
    list(LENGTH TARGET_SOURCE _num)
    if(_num GREATER 1)
      message(FATAL_ERROR "Inconsistent ${target} library sources. "
        "Please delete ${CMAKE_BINARY_DIR}/_deps/src and re-run cmake")
    endif()
    file(REMOVE_RECURSE ${source_dir})
    file(RENAME ${TARGET_SOURCE} ${source_dir})
    file(WRITE ${source_dir}/.extracted-${hash} "")
  endif()
  if(NOT (cmakefile STREQUAL ""))
    file(COPY ${cmakefile} DESTINATION ${source_dir}/${cmakedir}/)
    get_filename_component(_cmakefile ${cmakefile} NAME)
    file(RENAME "${source_dir}/${cmakedir}/${_cmakefile}"
      "${source_dir}/${cmakedir}/CMakeLists.txt")
  endif()
  add_subdirectory("${source_dir}/${cmakedir}"
    "${CMAKE_BINARY_DIR}/_deps/${target}-build" EXCLUDE_FROM_ALL)
endfunction(ExternalCMakeProject)
