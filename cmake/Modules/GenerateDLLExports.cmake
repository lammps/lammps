# GenerateDLLExports.cmake -- write the module definition file with the exports of liblammps.dll
#
# Usage: ${CMAKE_COMMAND} -D NM=<nm executable> -D OBJECTS=<file listing the object files>
#                         -D DEF_FILE=<.def file to write> -P GenerateDLLExports.cmake
#
# The export table of a Windows DLL can hold at most 65535 symbols.  When cross-compiling
# LAMMPS with MinGW and many packages, including KOKKOS, exporting all global symbols with
# "-Wl,--export-all-symbols" exceeds this limit and linking fails with "export ordinal too
# large".  This script is run as a PRE_LINK command of the LAMMPS library instead and writes
# a .def file that exports the global symbols defined in the object files of the library
# with the following exceptions:
#  - symbols from linked static libraries (bundled or system libraries) are not included
#  - symbols starting with '.' (the linker does not export those, either)
#  - symbols containing "Kokkos", i.e. the KOKKOS package and the Kokkos templates
# None of those are needed by code using the DLL: the Python module and LAMMPS-GUI use only
# the C library interface, the LAMMPS executable and plugins use the C library interface and
# the C++ classes from the core and packages other than KOKKOS.
# Data symbols (global variables, vtables, typeinfo) are marked as DATA, so that the import
# library provides only the __imp_ entry for them, same as for automatically exported symbols.

cmake_minimum_required(VERSION 3.27)

file(STRINGS "${OBJECTS}" objects)
list(LENGTH objects num_objects)
set(symbols)
# process the object files in batches to keep the command line short
set(batch_size 100)
foreach(first RANGE 0 ${num_objects} ${batch_size})
  list(SUBLIST objects ${first} ${batch_size} batch)
  if(NOT batch)
    continue()
  endif()
  execute_process(COMMAND ${NM} -g --defined-only -P ${batch}
    OUTPUT_VARIABLE output
    RESULT_VARIABLE result)
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "Running ${NM} on the LAMMPS library objects failed")
  endif()
  # nm -P prints "name type value size": keep global code (T) and data (D, R, B) symbols
  string(REPLACE "\n" ";" lines "${output}")
  list(FILTER lines INCLUDE REGEX "^[^. ][^ ]* [TDRB] ")
  list(FILTER lines EXCLUDE REGEX "Kokkos")
  list(APPEND symbols ${lines})
endforeach()

list(TRANSFORM symbols REPLACE "^([^ ]+) T .*$" "  \\1")
list(TRANSFORM symbols REPLACE "^([^ ]+) [DRB] .*$" "  \\1 DATA")
list(REMOVE_DUPLICATES symbols)
list(SORT symbols)
list(LENGTH symbols num_symbols)
if(num_symbols GREATER 65535)
  message(FATAL_ERROR "Too many symbols to export from the LAMMPS DLL: ${num_symbols} > 65535")
endif()
message(STATUS "Exporting ${num_symbols} symbols from the LAMMPS DLL")
list(JOIN symbols "\n" exports)
file(WRITE "${DEF_FILE}" "EXPORTS\n${exports}\n")
