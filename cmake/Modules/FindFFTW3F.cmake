# Find the native single precision FFTW3 headers and libraries.
#
#  FFTW3F_INCLUDE_DIRS - where to find fftw3f.h, etc.
#  FFTW3F_LIBRARIES    - List of libraries when using fftw3f.
#  FFTW3F_OMP_LIBRARIES - List of libraries when using fftw3.
#  FFTW3F_FOUND        - True if fftw3f found.
#

find_package(PkgConfig)

pkg_check_modules(PC_FFTW3F fftw3f)
find_path(FFTW3F_INCLUDE_DIR fftw3.h HINTS ${PC_FFTW3F_INCLUDE_DIRS})
find_library(FFTW3F_LIBRARY NAMES fftw3f HINTS ${PC_FFTW3F_LIBRARY_DIRS})
find_library(FFTW3F_OMP_LIBRARY NAMES fftw3f_omp HINTS ${PC_FFTW3F_LIBRARY_DIRS})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set FFTW3F_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(FFTW3F DEFAULT_MSG FFTW3F_LIBRARY FFTW3F_INCLUDE_DIR)

# Copy the results to the output variables and target.
if(FFTW3F_FOUND)
  set(FFTW3F_LIBRARIES ${FFTW3F_LIBRARY} )
  set(FFTW3F_INCLUDE_DIRS ${FFTW3F_INCLUDE_DIR} )

  if(NOT TARGET FFTW3F::FFTW3F)
    add_library(FFTW3F::FFTW3F UNKNOWN IMPORTED)
    set_target_properties(FFTW3F::FFTW3F PROPERTIES
      IMPORTED_LINK_INTERFACE_LANGUAGES "C"
      IMPORTED_LOCATION "${FFTW3F_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${FFTW3F_INCLUDE_DIRS}")
  endif()
  if(FFTW3F_OMP_LIBRARY)
    set(FFTW3F_OMP_LIBRARIES ${FFTW3F_OMP_LIBRARY})
    if(NOT TARGET FFTW3F::FFTW3F_OMP)
      add_library(FFTW3F::FFTW3F_OMP UNKNOWN IMPORTED)
      set_target_properties(FFTW3F::FFTW3F_OMP PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "C"
	IMPORTED_LOCATION "${FFTW3F_OMP_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${FFTW3F_INCLUDE_DIRS}")
    endif()
  endif()
endif()

mark_as_advanced(FFTW3F_INCLUDE_DIR FFTW3F_LIBRARY FFTW3F_OMP_LIBRARY)
