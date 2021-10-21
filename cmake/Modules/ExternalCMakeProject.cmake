# Build a CMake based external library as subdirectory.
# The sources will be unpacked to ${CMAKE_BINARY_DIR}/_deps/${target}-src
# The binaries will be built in ${CMAKE_BINARY_DIR}/_deps/${target}-build
#
function(ExternalCMakeProject target url hash basedir cmakedir)
  # change settings locally
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)

  if(CMAKE_VERSION VERSION_LESS 3.14)
    get_filename_component(archive ${url} NAME)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/_deps/src)
    message(STATUS "Downloading ${url}")
    file(DOWNLOAD ${url} ${CMAKE_BINARY_DIR}/_deps/${archive} EXPECTED_HASH MD5=${hash} SHOW_PROGRESS)
    message(STATUS "Unpacking and configuring ${archive}")
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_BINARY_DIR}/_deps/${archive}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/_deps/src)
    file(GLOB TARGET_SOURCE "${CMAKE_BINARY_DIR}/_deps/src/${basedir}*")
    list(LENGTH TARGET_SOURCE _num)
    if(_num GREATER 1)
      message(FATAL_ERROR "Inconsistent ${target} library sources. "
        "Please delete ${CMAKE_BINARY_DIR}/_deps/src and re-run cmake")
    endif()
    file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/_deps/${target}-src)
    file(RENAME ${TARGET_SOURCE} ${CMAKE_BINARY_DIR}/_deps/${target}-src)
    add_subdirectory("${CMAKE_BINARY_DIR}/_deps/${target}-src/${cmakedir}"
      "${CMAKE_BINARY_DIR}/_deps/${target}-build")
  else()
    include(FetchContent)
    message(STATUS "Downloading ${url}")
    FetchContent_Declare(${target} URL ${url} URL_HASH MD5=${hash} SOURCE_SUBDIR ${cmakedir})
    message(STATUS "Unpacking and configuring ${archive}")
    FetchContent_MakeAvailable(${target})
  endif()
endfunction(ExternalCMakeProject)
