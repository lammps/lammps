###############################################################################
# Testing
###############################################################################
option(ENABLE_TESTING "Enable testing" OFF)
if(ENABLE_TESTING)
  enable_testing()
  option(LAMMPS_TESTING_SOURCE_DIR "Location of lammps-testing source directory" "")
  option(LAMMPS_TESTING_GIT_TAG    "Git tag of lammps-testing" "master")
  mark_as_advanced(LAMMPS_TESTING_SOURCE_DIR LAMMPS_TESTING_GIT_TAG)

  if (CMAKE_VERSION VERSION_GREATER "3.10.3" AND NOT LAMMPS_TESTING_SOURCE_DIR)
    include(FetchContent)

    FetchContent_Declare(lammps-testing
      GIT_REPOSITORY https://github.com/lammps/lammps-testing.git
      GIT_TAG ${LAMMPS_TESTING_GIT_TAG}
    )

    FetchContent_GetProperties(lammps-testing)
    if(NOT lammps-testing_POPULATED)
      message(STATUS "Downloading tests...")
      FetchContent_Populate(lammps-testing)
    endif()

    set(LAMMPS_TESTING_SOURCE_DIR ${lammps-testing_SOURCE_DIR})
  elseif(NOT LAMMPS_TESTING_SOURCE_DIR)
    message(WARNING "Full test-suite requires CMake >= 3.11 or copy of\n"
                    "https://github.com/lammps/lammps-testing in LAMMPS_TESTING_SOURCE_DIR")
  endif()

  add_test(NAME ShowHelp COMMAND $<TARGET_FILE:lmp> -help)

  if(EXISTS ${LAMMPS_TESTING_SOURCE_DIR})
    message(STATUS "Running test discovery...")

    file(GLOB_RECURSE TEST_SCRIPTS ${LAMMPS_TESTING_SOURCE_DIR}/tests/core/*/in.*)
    foreach(script_path ${TEST_SCRIPTS})
      get_filename_component(TEST_NAME ${script_path} EXT)
      get_filename_component(SCRIPT_NAME ${script_path} NAME)
      get_filename_component(PARENT_DIR ${script_path} DIRECTORY)
      string(SUBSTRING ${TEST_NAME} 1 -1 TEST_NAME)
      string(REPLACE "-" "_" TEST_NAME ${TEST_NAME})
      string(REPLACE "+" "_" TEST_NAME ${TEST_NAME})
      set(TEST_NAME "test_core_${TEST_NAME}_serial")
      add_test(NAME ${TEST_NAME} COMMAND $<TARGET_FILE:lmp> -in ${SCRIPT_NAME} WORKING_DIRECTORY ${PARENT_DIR})
    endforeach()
    list(LENGTH TEST_SCRIPTS NUM_TESTS)

    message(STATUS "Found ${NUM_TESTS} tests.")
  endif()
endif()
