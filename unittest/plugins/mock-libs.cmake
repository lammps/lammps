# Initial cache file (cmake -C) for testing the plugins of packages with external libraries
# (examples/PACKAGES/*/plugin) without spending time on compiling those libraries.  It packs the
# stand-ins mock-plumed and mock-mbx next to this file into archives in the build folder and sets
# the download URLs and checksums of PLUMED and MBX to them.  A stand-in builds a small library
# with the same interface (configure script, make, files and their locations) whose functions do
# nothing: the plugins can be compiled, linked and loaded, but do not work.

# the build folder does not exist yet when CMake reads this file
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR})
foreach(_mock plumed mbx)
  set(_archive ${CMAKE_BINARY_DIR}/mock-${_mock}.tar.gz)
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar czf ${_archive} mock-${_mock}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR} RESULT_VARIABLE _result)
  if(NOT _result EQUAL 0)
    message(FATAL_ERROR "Could not create ${_archive}")
  endif()
  file(SHA256 ${_archive} _sha256)
  string(TOUPPER ${_mock} _name)
  if(_mock STREQUAL "mbx")
    set(_name MBXLIB)
  endif()
  set(${_name}_URL "file://${_archive}" CACHE STRING "" FORCE)
  set(${_name}_SHA256 "${_sha256}" CACHE STRING "" FORCE)
endforeach()

set(DOWNLOAD_PLUMED ON CACHE BOOL "" FORCE)
set(DOWNLOAD_MBX ON CACHE BOOL "" FORCE)
