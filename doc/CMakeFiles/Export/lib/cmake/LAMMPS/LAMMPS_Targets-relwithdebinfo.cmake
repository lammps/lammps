#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "LAMMPS::lmp" for configuration "RelWithDebInfo"
set_property(TARGET LAMMPS::lmp APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(LAMMPS::lmp PROPERTIES
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/bin/lmp"
  )

list(APPEND _IMPORT_CHECK_TARGETS LAMMPS::lmp )
list(APPEND _IMPORT_CHECK_FILES_FOR_LAMMPS::lmp "${_IMPORT_PREFIX}/bin/lmp" )

# Import target "LAMMPS::lammps" for configuration "RelWithDebInfo"
set_property(TARGET LAMMPS::lammps APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(LAMMPS::lammps PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELWITHDEBINFO "Python::Python"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/liblammps.0.dylib"
  IMPORTED_SONAME_RELWITHDEBINFO "@rpath/liblammps.0.dylib"
  )

list(APPEND _IMPORT_CHECK_TARGETS LAMMPS::lammps )
list(APPEND _IMPORT_CHECK_FILES_FOR_LAMMPS::lammps "${_IMPORT_PREFIX}/lib/liblammps.0.dylib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
