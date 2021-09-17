# - Prevent in-source builds.
# https://stackoverflow.com/questions/1208681/with-cmake-how-would-you-disable-in-source-builds/

function(prevent_in_source_builds)
  # make sure the user doesn't play dirty with symlinks
  get_filename_component(srcdir "${CMAKE_SOURCE_DIR}" REALPATH)
  get_filename_component(srcdir2 "${CMAKE_SOURCE_DIR}/.." REALPATH)
  get_filename_component(srcdir3 "${CMAKE_SOURCE_DIR}/../src" REALPATH)
  get_filename_component(bindir "${CMAKE_BINARY_DIR}" REALPATH)

  # disallow in-source builds
  if(("${srcdir}" STREQUAL "${bindir}") OR ("${srcdir2}" STREQUAL "${bindir}") OR ("${srcdir3}" STREQUAL "${bindir}"))
    message(FATAL_ERROR "\

CMake must not to be run in the source directory. \
Rather create a dedicated build directory and run CMake there. \
To clean up after this aborted in-place compilation:
  rm -r CMakeCache.txt CMakeFiles
")
  endif()
endfunction()

prevent_in_source_builds()
