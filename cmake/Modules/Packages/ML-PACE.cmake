if(NOT DEFINED LOCAL_ML-PACE)
  # standard user scenario: get latest ML-PACE from github
  set(PACELIB_URL "https://github.com/ICAMS/lammps-user-pace/archive/refs/tags/v.2021.10.25.fix2.tar.gz" CACHE STRING "URL for PACE evaluator library sources")

  set(PACELIB_MD5 "32394d799bc282bb57696c78c456e64f" CACHE STRING "MD5 checksum of PACE evaluator library tarball")
  mark_as_advanced(PACELIB_URL)
  mark_as_advanced(PACELIB_MD5)

#  message("-- [ML-PACE]: Downloading from  ${PACELIB_URL}")
  # download library sources to build folder
  file(DOWNLOAD ${PACELIB_URL} ${CMAKE_BINARY_DIR}/libpace.tar.gz EXPECTED_HASH MD5=${PACELIB_MD5}) #SHOW_PROGRESS

  # uncompress downloaded sources
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E remove_directory lammps-user-pace*
    COMMAND ${CMAKE_COMMAND} -E tar xzf libpace.tar.gz
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
  get_newest_file(${CMAKE_BINARY_DIR}/lammps-user-pace-* lib-pace)
else()
  # developer scenario: use local ML-PACE

  if(EXISTS ${LOCAL_ML-PACE})
    set(lib-pace ${LOCAL_ML-PACE})
    message("-- [ML-PACE]: Using local folder ${lib-pace}")
  else()
    message(FATAL_ERROR "-- [ML-PACE]: ERROR! Folder ${LOCAL_ML-PACE} does not exist")
  endif()
endif()

##yaml
# enforce building libyaml-cpp as static library and turn off optional features
set(YAML_BUILD_SHARED_LIBS OFF)
set(YAML_CPP_BUILD_CONTRIB OFF)
set(YAML_CPP_BUILD_TOOLS OFF)
add_subdirectory(${lib-pace}/yaml-cpp build-yaml-cpp)
set(YAML_CPP_INCLUDE_DIR ${lib-pace}/yaml-cpp/include)

## cnpy
set(CNPY_PATH ${lib-pace}/cnpy)
set(CNPY_INCLUDE_PATH ${CNPY_PATH})
set(CNPY_SRC ${CNPY_PATH}/cnpy.cpp)
add_library(cnpy-static STATIC ${CNPY_SRC})
set_target_properties(cnpy-static PROPERTIES LINKER_LANGUAGE CXX)

## winger-cpp
# this is header-only library
set(WIGNER_PATH ${lib-pace}/wigner-cpp)
set(WIGNER_INCLUDE_PATH ${WIGNER_PATH}/include/wigner)

# ML-PACE sources
file(GLOB PACE_EVALUATOR_INCLUDE_DIR ${lib-pace}/ML-PACE)
file(GLOB PACE_EVALUATOR_SOURCES ${lib-pace}/ML-PACE/*.cpp)
list(FILTER PACE_EVALUATOR_SOURCES EXCLUDE REGEX pair_pace*.cpp)

add_library(pace STATIC ${PACE_EVALUATOR_SOURCES})
set_target_properties(pace PROPERTIES CXX_EXTENSIONS ON OUTPUT_NAME lammps_pace${LAMMPS_MACHINE})
target_include_directories(pace PUBLIC ${PACE_EVALUATOR_INCLUDE_DIR} ${YAML_CPP_INCLUDE_DIR} ${CNPY_INCLUDE_PATH} ${WIGNER_INCLUDE_PATH})
target_compile_definitions(pace PUBLIC EXTRA_C_PROJECTIONS)
target_link_libraries(pace PRIVATE yaml-cpp-pace cnpy-static)
if(CMAKE_PROJECT_NAME STREQUAL "lammps")
  target_link_libraries(lammps PRIVATE pace)
endif()
