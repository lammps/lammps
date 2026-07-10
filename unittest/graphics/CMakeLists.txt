# -*- CMake -*- file for tests of dump image and GRAPHICS package rendering

function(extract_tags out yaml_file)
  file(STRINGS ${yaml_file} TAGS_LINE REGEX "^tags:")
  string(REPLACE "tags:" "" TAGS_LINE "${TAGS_LINE}")
  string(REGEX MATCHALL "[^, \t\r\n]+" TAGS "${TAGS_LINE}")
  set(${out} "${TAGS}" PARENT_SCOPE)
endfunction()

function(extract_nprocs out yaml_file)
  file(STRINGS ${yaml_file} NPROCS_LINE REGEX "^nprocs:")
  string(REPLACE "nprocs:" "" NPROCS_LINE "${NPROCS_LINE}")
  string(STRIP "${NPROCS_LINE}" NPROCS_LINE)
  if(NOT NPROCS_LINE)
    set(NPROCS_LINE 1)
  endif()
  set(${out} "${NPROCS_LINE}" PARENT_SCOPE)
endfunction()

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  string(REPLACE "/" "\\\\" TEST_INPUT_FOLDER "${CMAKE_CURRENT_SOURCE_DIR}/tests")
else()
  set(TEST_INPUT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR}/tests)
endif()

add_executable(test_dump_image
  ../yaml_writer.cpp
  test_config_reader.cpp
  test_dump_image.cpp)
if(YAML_FOUND)
  target_compile_definitions(test_dump_image PRIVATE TEST_INPUT_FOLDER=${TEST_INPUT_FOLDER})
else()
  # we always use static linkage with local compiled libyaml
  target_compile_definitions(test_dump_image PRIVATE TEST_INPUT_FOLDER=${TEST_INPUT_FOLDER} YAML_DECLARE_STATIC)
endif()
target_include_directories(test_dump_image PRIVATE "${LAMMPS_SOURCE_DIR};${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(test_dump_image PRIVATE GTest::GMock Yaml::Yaml lammps)

# propagate sanitizer options to test tools
if(ENABLE_SANITIZER AND (NOT (ENABLE_SANITIZER STREQUAL "none")))
  target_compile_options(test_dump_image PUBLIC -fsanitize=${ENABLE_SANITIZER})
  target_link_options(test_dump_image PUBLIC -fsanitize=${ENABLE_SANITIZER})
endif()

# prepare environment for testers
set(GRAPHICS_TEST_ENVIRONMENT "OMP_PROC_BIND=false")
list(APPEND GRAPHICS_TEST_ENVIRONMENT "OMP_NUM_THREADS=4")
list(APPEND GRAPHICS_TEST_ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR}")
get_property(BUILD_IS_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(BUILD_IS_MULTI_CONFIG)
  set(LAMMPS_LIB_PATH ${CMAKE_BINARY_DIR}/$<CONFIG>)
else()
  set(LAMMPS_LIB_PATH ${CMAKE_BINARY_DIR})
endif()
if(APPLE)
  list(APPEND GRAPHICS_TEST_ENVIRONMENT "DYLD_LIBRARY_PATH=${LAMMPS_LIB_PATH}:$ENV{DYLD_LIBRARY_PATH}")
elseif(WIN32)
  list(APPEND GRAPHICS_TEST_ENVIRONMENT "LAMMPSDLLPATH=${LAMMPS_LIB_PATH}")
else()
  list(APPEND GRAPHICS_TEST_ENVIRONMENT "LD_LIBRARY_PATH=${LAMMPS_LIB_PATH}:$ENV{LD_LIBRARY_PATH}")
endif()

# register all YAML files; tests with nprocs > 1 need an MPI launcher
file(GLOB GRAPHICS_TESTS LIST_DIRECTORIES false CONFIGURE_DEPENDS ${TEST_INPUT_FOLDER}/*.yaml)
foreach(TEST ${GRAPHICS_TESTS})
  string(REGEX REPLACE "^.*tests.(.*)\.yaml" "GRAPHICS:\\1" TNAME ${TEST})
  extract_tags(TEST_TAGS ${TEST})
  list(FIND TEST_TAGS "no${CMAKE_SYSTEM_NAME}" _SKIPME)
  if(_SKIPME GREATER_EQUAL 0)
    continue()
  endif()
  extract_nprocs(TEST_NPROCS ${TEST})
  if(TEST_NPROCS GREATER 1)
    if(BUILD_MPI AND MPIEXEC_EXECUTABLE)
      add_test(NAME ${TNAME}
               COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${TEST_NPROCS}
                       $<TARGET_FILE:test_dump_image> ${TEST})
    else()
      continue()
    endif()
  else()
    add_test(NAME ${TNAME} COMMAND test_dump_image ${TEST})
  endif()
  set_tests_properties(${TNAME} PROPERTIES ENVIRONMENT "${GRAPHICS_TEST_ENVIRONMENT}")
  set_tests_properties(${TNAME} PROPERTIES LABELS "${TEST_TAGS}")
endforeach()
