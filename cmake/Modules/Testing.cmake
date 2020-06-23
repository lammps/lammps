###############################################################################
# Testing
###############################################################################
option(ENABLE_TESTING "Enable testing" OFF)
if(ENABLE_TESTING)
  find_program(VALGRIND_BINARY NAMES valgrind)
  set(VALGRIND_DEFAULT_OPTIONS "--leak-check=full --show-leak-kinds=all --track-origins=yes")

  if(BUILD_MPI)
    set(VALGRIND_DEFAULT_OPTIONS "${VALGRIND_DEFAULT_OPTIONS} --suppressions=${LAMMPS_TOOLS_DIR}/valgrind/OpenMP.supp")
  endif()

  if(BUILD_OMP)
    set(VALGRIND_DEFAULT_OPTIONS "${VALGRIND_DEFAULT_OPTIONS} --suppressions=${LAMMPS_TOOLS_DIR}/valgrind/OpenMP.supp")
  endif()

  if(PKG_PYTHON)
    set(VALGRIND_DEFAULT_OPTIONS "${VALGRIND_DEFAULT_OPTIONS} --suppressions=${LAMMPS_TOOLS_DIR}/valgrind/Python3.supp")
  endif()

  set(MEMORYCHECK_COMMAND "${VALGRIND_BINARY}" CACHE FILEPATH "Memory Check Command")
  set(MEMORYCHECK_COMMAND_OPTIONS "${VALGRIND_DEFAULT_OPTIONS}" CACHE STRING "Memory Check Command Options")

  include(CTest)

  enable_testing()
  get_filename_component(LAMMPS_UNITTEST_DIR ${LAMMPS_SOURCE_DIR}/../unittest ABSOLUTE)
  get_filename_component(LAMMPS_UNITTEST_BIN ${CMAKE_BINARY_DIR}/unittest ABSOLUTE)
  add_subdirectory(${LAMMPS_UNITTEST_DIR} ${LAMMPS_UNITTEST_BIN})
endif()
