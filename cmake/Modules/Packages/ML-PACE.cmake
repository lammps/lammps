# try to find system libyaml-cpp v.0.6.3 library
find_package(yaml-cpp 0.6.3 EXACT QUIET)


if(yaml-cpp_FOUND)
    message(STATUS "Library yaml-cpp(v0.6.3)  found")    
    #message(STATUS "YAML_CPP_LIBRARY_DIRS=${YAML_CPP_LIBRARY_DIRS}")
    
#    get_cmake_property(_variableNames VARIABLES)
#    list (SORT _variableNames)
#    foreach (_variableName ${_variableNames})
#        message(STATUS "${_variableName}=${${_variableName}}")
#    endforeach()

    find_path(YAML_CPP_INCLUDE_DIR
            NAMES yaml.h
            PATHS ${YAML_CPP_INCLUDE_DIRS})

    if(NOT YAML_CPP_LIBRARIES)
        find_library(YAML_CPP_LIBRARY
                NAMES yaml-cpp
                PATHS ${YAML_CPP_LIBRARY_DIRS})
    else()
        set(YAML_CPP_LIBRARY ${YAML_CPP_LIBRARIES})
    endif()
else()
    message(STATUS "Library yaml-cpp(v0.6.3) not found, downloading")
    
    set(YAML_063_URL "https://github.com/jbeder/yaml-cpp/archive/refs/tags/yaml-cpp-0.6.3.tar.gz" CACHE STRING "URL for yaml-cpp (v.0.6.3) library sources")    
    set(YAML_063_MD5 "b45bf1089a382e81f6b661062c10d0c2" CACHE STRING "MD5 checksum of yaml-cpp (v.0.6.3) library tarball")
    mark_as_advanced(YAML_063_URL)
    mark_as_advanced(YAML_063_MD5)

    # download library sources to build folder
    file(DOWNLOAD ${YAML_063_URL} ${CMAKE_BINARY_DIR}/yaml-cpp-0.6.3.tar.gz SHOW_PROGRESS EXPECTED_HASH MD5=${YAML_063_MD5})
    
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E remove_directory yaml-cpp-yaml-cpp-0.6.3*
        COMMAND ${CMAKE_COMMAND} -E tar xzf yaml-cpp-0.6.3.tar.gz
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    
    set(YAML_DIR ${CMAKE_BINARY_DIR}/yaml-cpp-yaml-cpp-0.6.3)
    
    add_subdirectory(${YAML_DIR} build-yaml-cpp EXCLUDE_FROM_ALL)
    set(YAML_CPP_INCLUDE_DIR ${YAML_DIR}/include)
    set(YAML_CPP_LIBRARY yaml-cpp)    
endif()



set(PACELIB_URL "https://github.com/ICAMS/lammps-user-pace/archive/refs/tags/v.2021.9.28.upd1.tar.gz" CACHE STRING "URL for PACE evaluator library sources")
set(PACELIB_MD5 "ec75bc491edd75e10560cdbf129d91a7" CACHE STRING "MD5 checksum of PACE evaluator library tarball")
mark_as_advanced(PACELIB_URL)
mark_as_advanced(PACELIB_MD5)

# download library sources to build folder
file(DOWNLOAD ${PACELIB_URL} ${CMAKE_BINARY_DIR}/libpace.tar.gz SHOW_PROGRESS EXPECTED_HASH MD5=${PACELIB_MD5})

# uncompress downloaded sources
execute_process(
  COMMAND ${CMAKE_COMMAND} -E remove_directory lammps-user-pace*
  COMMAND ${CMAKE_COMMAND} -E tar xzf libpace.tar.gz
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

file(GLOB lib-pace ${CMAKE_BINARY_DIR}/lammps-user-pace-*)

message(STATUS "YAML_CPP_INCLUDE_DIR=${YAML_CPP_INCLUDE_DIR}")
message(STATUS "YAML_CPP_LIBRARY=${YAML_CPP_LIBRARY}")

file(GLOB PACE_EVALUATOR_INCLUDE_DIR ${lib-pace}/ML-PACE)
file(GLOB PACE_EVALUATOR_SOURCES ${lib-pace}/ML-PACE/*.cpp)
list(FILTER PACE_EVALUATOR_SOURCES EXCLUDE REGEX pair_pace.cpp)

add_library(pace STATIC ${PACE_EVALUATOR_SOURCES})
set_target_properties(pace PROPERTIES CXX_EXTENSIONS ON OUTPUT_NAME lammps_pace${LAMMPS_MACHINE})
target_include_directories(pace PRIVATE ${YAML_CPP_INCLUDE_DIR})
target_include_directories(pace PUBLIC ${PACE_EVALUATOR_INCLUDE_DIR})

target_link_libraries(lammps PRIVATE pace)
target_link_libraries(lammps PRIVATE ${YAML_CPP_LIBRARY})
