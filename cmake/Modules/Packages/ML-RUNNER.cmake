# Enable Fortran language support
enable_language(Fortran)

# Build arguments.
option(DOWNLOAD_RUNNER "Force download and build of RuNNer. If this is OFF\
  the arguments RUNNER_LIB_DIR, RUNNER_LIB_NAME and RUNNER_SHARED_LIB must\
  be set to meaningful values." ON
)
option(RUNNER_SHARED_LIB "Use pre-compiled shared/dynamic RuNNer library. Only \
  considered if DOWNLOAD_RUNNER is OFF." OFF
)
set(RUNNER_LIB_DIR "" CACHE PATH "Directory containing \
  the RuNNer library. Only considered if DOWNLOAD_RUNNER is OFF."
)
set(RUNNER_LIB_NAME "libRuNNer_mpi" CACHE STRING "Name of the RuNNer library \
  (excluding file extension, this is controlled by RUNNER_SHARED_LIB). Only \
  considered if DOWNLOAD_RUNNER is OFF."
)

# Initialize RuNNer CMake arguments
set(USE_MKL OFF)
set(RUNNER_USE_FFTW_THREADED OFF)
set(RUNNER_DISABLE_FFTW OFF)

if (FFT STREQUAL "MKL")
  set(USE_MKL ON)
elseif(FFT STREQUAL "FFTW3")
  if (FFT_FFTW_THREADS)
    set(RUNNER_USE_FFTW_THREADED ON)
  else()
    set(RUNNER_USE_FFTW_THREADED OFF)
  endif()
else()
  set(RUNNER_DISABLE_FFTW ON)
  message(WARNING "No compatible FFT library found. The ML-RUNNER package only supports -DFFT=FFTW3 or -DFFT=MKL. RuNNer will fall back to being compiled WITHOUT FFT support, and certain features will not be available.")
endif()

if(BUILD_MPI)
  # Ensure the Fortran MPI components are found
  find_package(MPI REQUIRED COMPONENTS CXX Fortran)

  # Fail if the required mpi_f08 module is missing
  if(NOT MPI_Fortran_HAVE_F08_MODULE)
    message(FATAL_ERROR "RuNNer requires the Fortran 2008 MPI bindings (mpi_f08), but your MPI installation does not support it.")
  endif()
endif()

if(DOWNLOAD_RUNNER)
  # Include ExternalProject module
  include(ExternalProject)

  if(RUNNER_SHARED_LIB)
    set(RUNNER_LIB_EXT ".so")
    set(RUNNER_LIB_TYPE SHARED)
  else()
    set(RUNNER_LIB_EXT ".a")
    set(RUNNER_LIB_TYPE STATIC)
  endif()

  message(STATUS "DOWNLOAD_RUNNER is ON. Building RuNNer from source as a ${RUNNER_LIB_TYPE} library.")

  # Force using the shared library. RuNNer's cmake build always produces libRuNNer_mpi.so.
  if(BUILD_MPI)
    set(RUNNER_LIB_NAME "libRuNNer_mpi")
  else()
    set(RUNNER_LIB_NAME "libRuNNer")
  endif()

  set(RUNNER_LIB_FULL_NAME ${RUNNER_LIB_NAME}${RUNNER_LIB_EXT})

  ExternalProject_Add(runner_build
    GIT_REPOSITORY "https://gitlab.com/runner-suite/runner2.git"
    GIT_TAG "2.0.3_20260528"
    GIT_SHALLOW YES
    GIT_PROGRESS YES
    PATCH_COMMAND sed -i -e "s/--not --show-signature/--no-show-signature/" build_system/check_git.sh

    # Pass CMake arguments to RuNNer's build system
    CMAKE_ARGS
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -DCMAKE_POSITION_INDEPENDENT_CODE=yes
      -DBUILD_SHARED_LIBS=OFF
      -DENABLE_TESTS=no
      -DUSE_MKL=${USE_MKL}
      -DUSE_OPENMP=${BUILD_OMP}
      -DUSE_FFTW_THREADED=${RUNNER_USE_FFTW_THREADED}
      -DCMAKE_DISABLE_FIND_PACKAGE_FFTW=${RUNNER_DISABLE_FFTW}
      -DUSE_MPI=${BUILD_MPI}
      -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/runner_install
      -DBUILD_SHARED_LIB=${RUNNER_SHARED_LIB}

    # Define the build and install steps
    BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR>
    INSTALL_COMMAND ${CMAKE_COMMAND} --install <BINARY_DIR>

    # Specify the location of the built library
    BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/runner_install/RuNNer/lib/${RUNNER_LIB_FULL_NAME}
  )

  # Create an IMPORTED library target for RuNNer
  add_library(RuNNer::RuNNer ${RUNNER_LIB_TYPE} IMPORTED)
  set_target_properties(RuNNer::RuNNer PROPERTIES
    IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/runner_install/RuNNer/lib/${RUNNER_LIB_FULL_NAME}"
  )

  # Add a dependency to ensure RuNNer is built before the main target
  add_dependencies(lammps runner_build)

else()
  message(STATUS "DOWNLOAD_RUNNER is OFF. Looking for a pre-compiled RuNNer.")

  # Check if the directory specified by RUNNER_LIB_DIR exists, if the
  # user defined it.
  if(DEFINED RUNNER_LIB_DIR AND NOT RUNNER_LIB_DIR STREQUAL "")
    message(STATUS "RUNNER_LIB_DIR is set to: ${RUNNER_LIB_DIR}")

    if(NOT IS_DIRECTORY "${RUNNER_LIB_DIR}")
      message(FATAL_ERROR
        "The directory specified by RUNNER_LIB_DIR does not exist: ${RUNNER_LIB_DIR}"
      )
    endif()
  else()
    message(STATUS "RUNNER_LIB_DIR not set or empty; using default behaviour")
  endif()

  # Determine the correct file extension based on whether the library is shared or not.
  if(RUNNER_SHARED_LIB)
    set(RUNNER_LIB_EXT ".so")
    set(RUNNER_LIB_TYPE SHARED)
  else()
    set(RUNNER_LIB_EXT ".a")
    set(RUNNER_LIB_TYPE STATIC)
  endif()

  set(RUNNER_LIB_FULL_NAME ${RUNNER_LIB_NAME}${RUNNER_LIB_EXT})
  message(STATUS "Looking for the ${RUNNER_LIB_TYPE} library: ${RUNNER_LIB_FULL_NAME}")

  # Find the RuNNer library in the specified path
  find_library(RuNNer_LIBRARY
    NAMES ${RUNNER_LIB_FULL_NAME}
    HINTS
      "${RUNNER_LIB_DIR}"       # 1st priority: User's manual cache variable
      ENV LD_LIBRARY_PATH       # 2nd priority: Search the shell's LD_LIBRARY_PATH
      ENV LIBRARY_PATH          # 3rd priority: Search the compiler's link path
    PATH_SUFFIXES lib lib64     # Automatically look in <hint>/lib if not found in <hint>
  )

  if(RuNNer_LIBRARY)
    message(STATUS "Found RuNNer library: ${RuNNer_LIBRARY}")

    # Create an IMPORTED library target for the found RuNNer library
    add_library(RuNNer::RuNNer ${RUNNER_LIB_TYPE} IMPORTED)
    set_target_properties(RuNNer::RuNNer PROPERTIES
      IMPORTED_LOCATION "${RuNNer_LIBRARY}"
    )
  else()
    message(FATAL_ERROR "Could not find the RuNNer library in the specified \
      RUNNER_LIB_DIR. Searched for ${RUNNER_LIB_FULL_NAME} at ${RUNNER_LIB_DIR}"
    )
  endif()
endif()

# Link lammps to the found RuNNer library
target_link_libraries(lammps PRIVATE RuNNer::RuNNer ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

if(BUILD_MPI)
  # explicitly link Fortran MPI libraries so the C++ linker can resolve mpi_f08 symbols
  target_link_libraries(lammps PRIVATE MPI::MPI_Fortran)
endif()
