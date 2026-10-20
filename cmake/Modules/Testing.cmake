###############################################################################
# Testing
###############################################################################
option(ENABLE_TESTING "Enable testing" OFF)
if(ENABLE_TESTING)
  find_program(VALGRIND_BINARY NAMES valgrind)
  # generate custom suppression file
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lammps.supp "\n")
  file(GLOB VALGRIND_SUPPRESSION_FILES CONFIGURE_DEPENDS ${LAMMPS_TOOLS_DIR}/valgrind/[^.]*.supp)
  foreach(SUPP ${VALGRIND_SUPPRESSION_FILES})
    file(READ ${SUPP} SUPPRESSIONS)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/lammps.supp "${SUPPRESSIONS}")
  endforeach()
  # blocks that are still reachable when the process exits are not leaks. reporting
  # them buries the actual defects under thousands of one-time initializations in
  # OpenSSL, libcurl, CPython, and GoogleTest, plus everything that is still in use
  # when a command like "quit" or an error exit terminates the process on purpose.
  # --max-stackframe is raised because Fortran libraries (e.g. QUIP) place large
  # automatic arrays on the stack; with the default limit valgrind misreads the
  # resulting stack pointer jump as a stack switch and then floods the log with
  # millions of bogus invalid access and uninitialized value errors.
  set(VALGRIND_DEFAULT_OPTIONS "--leak-check=full --show-leak-kinds=definite,indirect,possible --track-origins=yes --max-stackframe=16777216 --suppressions=${CMAKE_BINARY_DIR}/lammps.supp")

  set(MEMORYCHECK_COMMAND "${VALGRIND_BINARY}" CACHE FILEPATH "Memory Check Command")
  set(MEMORYCHECK_COMMAND_OPTIONS "${VALGRIND_DEFAULT_OPTIONS}" CACHE STRING "Memory Check Command Options")

  # we need to build and link a LOT of tester executables, so it is worth checking if
  # a faster linker is available. requires GNU or Clang compiler, newer CMake.
  # also only verified with Fedora Linux > 30 and Ubuntu 18.04 or 22.04+(Ubuntu 20.04 fails)
  if((CMAKE_SYSTEM_NAME STREQUAL "Linux")
      AND ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")))
    if(((CMAKE_LINUX_DISTRO STREQUAL "Ubuntu") AND (CMAKE_DISTRO_VERSION VERSION_GREATER_EQUAL 22.04))
        OR ((CMAKE_LINUX_DISTRO STREQUAL "Fedora") AND (CMAKE_DISTRO_VERSION VERSION_GREATER 30)))
      include(CheckCXXCompilerFlag)
      set(LAMMPS_CUSTOM_LINKER_DEFAULT default)
      check_cxx_compiler_flag(--ld-path=${CMAKE_LINKER} HAVE_LD_PATH_FLAG)
      check_cxx_compiler_flag(-fuse-ld=mold HAVE_MOLD_LINKER_FLAG)
      check_cxx_compiler_flag(-fuse-ld=lld HAVE_LLD_LINKER_FLAG)
      check_cxx_compiler_flag(-fuse-ld=gold HAVE_GOLD_LINKER_FLAG)
      check_cxx_compiler_flag(-fuse-ld=bfd HAVE_BFD_LINKER_FLAG)
      find_program(HAVE_MOLD_LINKER_BIN ld.mold)
      find_program(HAVE_LLD_LINKER_BIN lld ld.lld)
      find_program(HAVE_GOLD_LINKER_BIN ld.gold)
      find_program(HAVE_BFD_LINKER_BIN ld.bfd)
      if(HAVE_MOLD_LINKER_FLAG AND HAVE_MOLD_LINKER_BIN)
        set(LAMMPS_CUSTOM_LINKER_DEFAULT mold)
      elseif(HAVE_LLD_LINKER_FLAG AND HAVE_LLD_LINKER_BIN)
        set(LAMMPS_CUSTOM_LINKER_DEFAULT lld)
      elseif(HAVE_GOLD_LINKER_FLAG AND HAVE_GOLD_LINKER_BIN)
        set(LAMMPS_CUSTOM_LINKER_DEFAULT gold)
      elseif(HAVE_BFD_LINKER_FLAG AND HAVE_BFD_LINKER_BIN)
        set(LAMMPS_CUSTOM_LINKER_DEFAULT bfd)
      endif()
      set(LAMMPS_CUSTOM_LINKER_VALUES mold lld gold bfd default)
      set(LAMMPS_CUSTOM_LINKER ${LAMMPS_CUSTOM_LINKER_DEFAULT} CACHE STRING "Choose a custom linker for faster linking (mold, lld, gold, bfd, default)")
      validate_option(LAMMPS_CUSTOM_LINKER LAMMPS_CUSTOM_LINKER_VALUES)
      mark_as_advanced(LAMMPS_CUSTOM_LINKER)
      if(NOT "${LAMMPS_CUSTOM_LINKER}" STREQUAL "default")
        target_link_options(lammps PUBLIC -fuse-ld=${LAMMPS_CUSTOM_LINKER})
      endif()
      if(HAVE_LD_PATH_FLAG)
        if("${LAMMPS_CUSTOM_LINKER}" STREQUAL "mold")
          target_link_options(lammps PUBLIC --ld-path=${HAVE_MOLD_LINKER_BIN})
        elseif("${LAMMPS_CUSTOM_LINKER}" STREQUAL "lld")
          target_link_options(lammps PUBLIC --ld-path=${HAVE_LLD_LINKER_BIN})
        elseif("${LAMMPS_CUSTOM_LINKER}" STREQUAL "gold")
          target_link_options(lammps PUBLIC --ld-path=${HAVE_GOLD_LINKER_BIN})
        elseif("${LAMMPS_CUSTOM_LINKER}" STREQUAL "bfd")
          target_link_options(lammps PUBLIC --ld-path=${HAVE_BFD_LINKER_BIN})
        endif()
      endif()
    endif()
  endif()

  get_filename_component(LAMMPS_UNITTEST_DIR ${LAMMPS_SOURCE_DIR}/../unittest ABSOLUTE)
  get_filename_component(LAMMPS_UNITTEST_BIN ${CMAKE_BINARY_DIR}/unittest ABSOLUTE)

  # With enabled GPU support, parallel tests may be dispatched to different NVIDIA GPUs through
  # the resource allocation feature of CTest and the CUDA_VISIBLE_DEVICES environment variable:
  # Every test requests one slot of a "gpus" resource, and a launcher script sets the
  # environment variable accordingly.  Available GPUs are listed in a resource spec file, which
  # is generated from the output of nvidia-smi if it does not exist. Without a resource spec file,
  # the tests run the same way as without this feature.
  set(LAMMPS_TEST_RESOURCE_SPEC_FILE ${CMAKE_BINARY_DIR}/ctest-resources.json CACHE FILEPATH
    "CTest resource spec file running GPU unit tests concurrently on multiple NVIDIA GPUS (empty to disable)")
  mark_as_advanced(LAMMPS_TEST_RESOURCE_SPEC_FILE)
  set(LAMMPS_TEST_LAUNCHER "")
  if(LAMMPS_TEST_RESOURCE_SPEC_FILE AND (PKG_GPU OR (PKG_KOKKOS AND Kokkos_ENABLE_CUDA)))
    # the TEST_LAUNCHER target property is required to apply the launcher to the tests
    if(CMAKE_VERSION VERSION_LESS 3.29 OR NOT UNIX)
      message(STATUS "Running unit tests across multiple GPUs requires CMake 3.29 or later on a Unix-like OS")
    else()
      if(NOT EXISTS ${LAMMPS_TEST_RESOURCE_SPEC_FILE})
        find_program(NVIDIA_SMI_EXECUTABLE nvidia-smi)
        mark_as_advanced(NVIDIA_SMI_EXECUTABLE)
        set(GPU_RESOURCES "")
        if(NVIDIA_SMI_EXECUTABLE)
          execute_process(COMMAND ${NVIDIA_SMI_EXECUTABLE} --query-gpu=index,name,pci.bus_id --format=csv,noheader
                          OUTPUT_VARIABLE GPU_INFO RESULT_VARIABLE GPU_STATUS ERROR_QUIET)
          if(GPU_STATUS EQUAL 0)
            string(REGEX MATCHALL "[^\r\n]+" GPU_INFO "${GPU_INFO}")
            foreach(GPU ${GPU_INFO})
              if(GPU MATCHES "^ *([0-9]+), *([^,\"]+), *([^,\"]+)$")
                list(APPEND GPU_RESOURCES "        { \"id\": \"${CMAKE_MATCH_1}\", \"slots\": 1, \"name\": \"${CMAKE_MATCH_2}\", \"pci_bus_id\": \"${CMAKE_MATCH_3}\" }")
              endif()
            endforeach()
          endif()
        endif()
        if(GPU_RESOURCES)
          list(LENGTH GPU_RESOURCES NUM_GPUS)
          list(JOIN GPU_RESOURCES ",\n" GPU_RESOURCES)
          message(STATUS "Writing CTest resource spec file for ${NUM_GPUS} GPU(s) to ${LAMMPS_TEST_RESOURCE_SPEC_FILE}")
          file(WRITE ${LAMMPS_TEST_RESOURCE_SPEC_FILE} "{
  \"comment\": [
    \"CTest resource spec file for running unit tests on multiple NVIDIA GPUs concurrently.\",
    \"Generated from the output of nvidia-smi.  Every test is assigned one slot of one GPU, and\",
    \"the GPU of each entry is selected by its 'id', which is the nvidia-smi index of the GPU.\",
  ],
  \"version\": { \"major\": 1, \"minor\": 0 },
  \"local\": [
    {
      \"gpus\": [
${GPU_RESOURCES}
      ]
    }
  ]
}
")
        else()
          message(STATUS "Found no NVIDIA GPUs with nvidia-smi for running unit tests concurrently")
        endif()
      endif()
      if(EXISTS ${LAMMPS_TEST_RESOURCE_SPEC_FILE})
        message(STATUS "Running unit tests concurrently on the GPUs listed in ${LAMMPS_TEST_RESOURCE_SPEC_FILE}")
        set(CTEST_RESOURCE_SPEC_FILE ${LAMMPS_TEST_RESOURCE_SPEC_FILE})
        set(LAMMPS_TEST_LAUNCHER /bin/sh ${LAMMPS_UNITTEST_DIR}/ctest-gpu-launcher.sh)
      endif()
    endif()
  endif()

  include(CTest)

  enable_testing()
  # tests with an executable target as command get the launcher through the TEST_LAUNCHER
  # target property, all others must have ${LAMMPS_TEST_LAUNCHER} prepended explicitly.
  # the launcher runs a test launcher set by the user in turn.
  if(LAMMPS_TEST_LAUNCHER)
    set(CMAKE_TEST_LAUNCHER ${LAMMPS_TEST_LAUNCHER} ${CMAKE_TEST_LAUNCHER})
  endif()
  add_subdirectory(${LAMMPS_UNITTEST_DIR} ${LAMMPS_UNITTEST_BIN})
  unset(CMAKE_TEST_LAUNCHER)

  # request one GPU for every test in the unittest folder and all its subfolders
  if(LAMMPS_TEST_LAUNCHER)
    function(request_test_gpus dir)
      get_property(TESTS DIRECTORY ${dir} PROPERTY TESTS)
      if(TESTS)
        set_tests_properties(${TESTS} DIRECTORY ${dir} PROPERTIES RESOURCE_GROUPS "gpus:1")
      endif()
      get_property(SUBDIRS DIRECTORY ${dir} PROPERTY SUBDIRECTORIES)
      foreach(SUBDIR ${SUBDIRS})
        request_test_gpus(${SUBDIR})
      endforeach()
    endfunction()
    request_test_gpus(${LAMMPS_UNITTEST_DIR})
  endif()
endif()

# Compiler specific features for testing
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  option(ENABLE_COVERAGE "Enable collecting code coverage data" OFF)
  mark_as_advanced(ENABLE_COVERAGE)
  if(ENABLE_COVERAGE)
    target_compile_options(lammps PUBLIC --coverage)
    target_link_options(lammps PUBLIC --coverage)
  endif()
endif()

#######################################
# add custom target for IWYU analysis
#######################################
set(ENABLE_IWYU OFF CACHE BOOL "Add 'iwyu' build target to call the include-what-you-use tool")
mark_as_advanced(ENABLE_IWYU)
if(ENABLE_IWYU)
  # enforce these settings
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Enable reporting compilation commands to compile_commands.json" FORCE)
  if(NOT ((CMAKE_CXX_COMPILER_ID STREQUAL "Clang") OR (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")))
    message(FATAL_ERROR "IWYU is only supported with Clang or GNU compilers")
  endif()
  # detect the "native" header folder so we can include them first
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} --print-search-dirs OUTPUT_VARIABLE IWYU_SEARCH_PATHS)
  string(REGEX REPLACE ".*libraries: *=([^:]+):.*" "\\1/include" IWYU_EXTRA_INCLUDE_DIR ${IWYU_SEARCH_PATHS})
  find_program(IWYU_EXE NAMES include-what-you-use iwyu)
  find_program(IWYU_TOOL NAMES iwyu_tool iwyu-tool iwyu_tool.py)
  if(IWYU_EXE AND IWYU_TOOL)
    add_custom_target(
      iwyu
      ${IWYU_TOOL} -o clang -p ${CMAKE_CURRENT_BINARY_DIR} -- -I${IWYU_EXTRA_INCLUDE_DIR} -Xiwyu --mapping_file=${CMAKE_CURRENT_SOURCE_DIR}/iwyu/iwyu-extra-map.imp
      COMMENT "Running IWYU")
    add_dependencies(iwyu lammps)
  else()
    message(FATAL_ERROR "To use IWYU you need the include-what-you-use/iwyu executable"
      "and the iwyu-tool/iwyu_tool script installed in your PATH")
  endif()
endif()

#######################################
# select code sanitizer options
#######################################
set(ENABLE_SANITIZER "none" CACHE STRING "Select a code sanitizer option (none (default), address, hwaddress, leak, thread, undefined)")
mark_as_advanced(ENABLE_SANITIZER)
set(ENABLE_SANITIZER_VALUES none address hwaddress leak thread undefined)
set_property(CACHE ENABLE_SANITIZER PROPERTY STRINGS ${ENABLE_SANITIZER_VALUES})
validate_option(ENABLE_SANITIZER ENABLE_SANITIZER_VALUES)
string(TOLOWER ${ENABLE_SANITIZER} ENABLE_SANITIZER)
if(NOT ENABLE_SANITIZER STREQUAL "none")
  if((${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU") OR (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
    target_compile_options(lammps PUBLIC -fsanitize=${ENABLE_SANITIZER})
    target_link_options(lammps PUBLIC -fsanitize=${ENABLE_SANITIZER})
  else()
    message(WARNING "ENABLE_SANITIZER option not supported by ${CMAKE_CXX_COMPILER_ID} compilers. Ignoring.")
    set(ENABLE_SANITIZER "none")
  endif()
endif()
