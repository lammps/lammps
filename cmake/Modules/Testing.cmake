###############################################################################
# Testing
###############################################################################
option(ENABLE_TESTING "Enable testing" OFF)
if(ENABLE_TESTING)
  find_program(VALGRIND_BINARY NAMES valgrind)
  # generate custom suppression file
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lammps.supp "\n")
  file(GLOB VALGRIND_SUPPRESSION_FILES ${LAMMPS_TOOLS_DIR}/valgrind/[^.]*.supp)
  foreach(SUPP ${VALGRIND_SUPPRESSION_FILES})
    file(READ ${SUPP} SUPPRESSIONS)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/lammps.supp "${SUPPRESSIONS}")
  endforeach()
  set(VALGRIND_DEFAULT_OPTIONS "--leak-check=full --show-leak-kinds=all --track-origins=yes --suppressions=${CMAKE_BINARY_DIR}/lammps.supp")

  set(MEMORYCHECK_COMMAND "${VALGRIND_BINARY}" CACHE FILEPATH "Memory Check Command")
  set(MEMORYCHECK_COMMAND_OPTIONS "${VALGRIND_DEFAULT_OPTIONS}" CACHE STRING "Memory Check Command Options")

  include(CTest)

  enable_testing()
  get_filename_component(LAMMPS_UNITTEST_DIR ${LAMMPS_SOURCE_DIR}/../unittest ABSOLUTE)
  get_filename_component(LAMMPS_UNITTEST_BIN ${CMAKE_BINARY_DIR}/unittest ABSOLUTE)
  add_subdirectory(${LAMMPS_UNITTEST_DIR} ${LAMMPS_UNITTEST_BIN})
endif()
