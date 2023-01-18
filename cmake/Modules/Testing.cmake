###############################################################################
# Testing
###############################################################################
option(ENABLE_TESTING "Enable testing" OFF)
if(ENABLE_TESTING)
  find_program(VALGRIND_BINARY NAMES valgrind)
  # generate custom suppression file
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lammps.supp "\n")
  file(GLOB VALGRIND_SUPPRESSION_FILES ${CONFIGURE_DEPENDS} ${LAMMPS_TOOLS_DIR}/valgrind/[^.]*.supp)
  foreach(SUPP ${VALGRIND_SUPPRESSION_FILES})
    file(READ ${SUPP} SUPPRESSIONS)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/lammps.supp "${SUPPRESSIONS}")
  endforeach()
  set(VALGRIND_DEFAULT_OPTIONS "--leak-check=full --show-leak-kinds=all --track-origins=yes --suppressions=${CMAKE_BINARY_DIR}/lammps.supp")

  set(MEMORYCHECK_COMMAND "${VALGRIND_BINARY}" CACHE FILEPATH "Memory Check Command")
  set(MEMORYCHECK_COMMAND_OPTIONS "${VALGRIND_DEFAULT_OPTIONS}" CACHE STRING "Memory Check Command Options")

  # we need to build and link a LOT of tester executables, so it is worth checking if
  # a faster linker is available. requires GNU or Clang compiler, newer CMake.
  # also only verified with Fedora Linux > 30 and Ubuntu <= 18.04 (Ubuntu 20.04 fails)
  if((CMAKE_SYSTEM_NAME STREQUAL "Linux") AND (CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
      AND ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")))
    if(((CMAKE_LINUX_DISTRO STREQUAL "Ubuntu") AND (CMAKE_DISTRO_VERSION VERSION_LESS_EQUAL 18.04))
        OR ((CMAKE_LINUX_DISTRO STREQUAL "Fedora") AND (CMAKE_DISTRO_VERSION VERSION_GREATER 30)))
      include(CheckCXXCompilerFlag)
      set(CMAKE_CUSTOM_LINKER_DEFAULT default)
      check_cxx_compiler_flag(-fuse-ld=lld HAVE_LLD_LINKER_FLAG)
      check_cxx_compiler_flag(-fuse-ld=gold HAVE_GOLD_LINKER_FLAG)
      check_cxx_compiler_flag(-fuse-ld=bfd HAVE_BFD_LINKER_FLAG)
      find_program(HAVE_LLD_LINKER_BIN lld ld.lld)
      find_program(HAVE_GOLD_LINKER_BIN ld.gold)
      find_program(HAVE_BFD_LINKER_BIN ld.bfd)
      if(HAVE_LLD_LINKER_FLAG AND HAVE_LLD_LINKER_BIN)
        set(CMAKE_CUSTOM_LINKER_DEFAULT lld)
      elseif(HAVE_GOLD_LINKER_FLAG AND HAVE_GOLD_LINKER_BIN)
        set(CMAKE_CUSTOM_LINKER_DEFAULT gold)
      elseif(HAVE_BFD_LINKER_FLAG AND HAVE_BFD_LINKER_BIN)
        set(CMAKE_CUSTOM_LINKER_DEFAULT bfd)
      endif()
      set(CMAKE_CUSTOM_LINKER_VALUES lld gold bfd default)
      set(CMAKE_CUSTOM_LINKER ${CMAKE_CUSTOM_LINKER_DEFAULT} CACHE STRING "Choose a custom linker for faster linking (lld, gold, bfd, default)")
      validate_option(CMAKE_CUSTOM_LINKER CMAKE_CUSTOM_LINKER_VALUES)
      mark_as_advanced(CMAKE_CUSTOM_LINKER)
      if(NOT "${CMAKE_CUSTOM_LINKER}" STREQUAL "default")
        target_link_options(lammps PUBLIC -fuse-ld=${CMAKE_CUSTOM_LINKER})
      endif()
    endif()
  endif()

  include(CTest)

  enable_testing()
  get_filename_component(LAMMPS_UNITTEST_DIR ${LAMMPS_SOURCE_DIR}/../unittest ABSOLUTE)
  get_filename_component(LAMMPS_UNITTEST_BIN ${CMAKE_BINARY_DIR}/unittest ABSOLUTE)
  add_subdirectory(${LAMMPS_UNITTEST_DIR} ${LAMMPS_UNITTEST_BIN})
endif()

# Compiler specific features for testing
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  option(ENABLE_COVERAGE "Enable collecting code coverage data" OFF)
  mark_as_advanced(ENABLE_COVERAGE)
  if(ENABLE_COVERAGE)
    if(CMAKE_VERSION VERSION_LESS 3.13)
      if(CMAKE_CXX_FLAGS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
      else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_${CMAKE_BUILD_TYPE}_FLAGS} --coverage")
      endif()
    else()
      target_compile_options(lammps PUBLIC --coverage)
      target_link_options(lammps PUBLIC --coverage)
    endif()
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
set(ENABLE_SANITIZER "none" CACHE STRING "Select a code sanitizer option (none (default), address, leak, thread, undefined)")
mark_as_advanced(ENABLE_SANITIZER)
set(ENABLE_SANITIZER_VALUES none address leak thread undefined)
set_property(CACHE ENABLE_SANITIZER PROPERTY STRINGS ${ENABLE_SANITIZER_VALUES})
validate_option(ENABLE_SANITIZER ENABLE_SANITIZER_VALUES)
string(TOLOWER ${ENABLE_SANITIZER} ENABLE_SANITIZER)
if(NOT ENABLE_SANITIZER STREQUAL "none")
  if((${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU") OR (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
    if(CMAKE_VERSION VERSION_LESS 3.13)
      if(CMAKE_CXX_FLAGS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=${ENABLE_SANITIZER}")
      else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_${CMAKE_BUILD_TYPE}_FLAGS} -fsanitize=${ENABLE_SANITIZER}")
      endif()
    else()
      target_compile_options(lammps PUBLIC -fsanitize=${ENABLE_SANITIZER})
      target_link_options(lammps PUBLIC -fsanitize=${ENABLE_SANITIZER})
    endif()
  else()
    message(WARNING "ENABLE_SANITIZER option not supported by ${CMAKE_CXX_COMPILER_ID} compilers. Ignoring.")
    set(ENABLE_SANITIZER "none")
  endif()
endif()
