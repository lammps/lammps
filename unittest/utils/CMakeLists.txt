# -*- CMake -*- file for tests of utily functions and classes in LAMMPS

add_executable(test_tokenizer test_tokenizer.cpp)
target_link_libraries(test_tokenizer PRIVATE lammps GTest::GMockMain)
add_test(NAME Tokenizer COMMAND test_tokenizer)

add_executable(test_stl_reader test_stl_reader.cpp)
target_link_libraries(test_stl_reader PRIVATE lammps GTest::GMockMain)
add_test(NAME STLReader COMMAND test_stl_reader)

add_executable(test_mempool test_mempool.cpp)
target_link_libraries(test_mempool PRIVATE lammps GTest::GMockMain)
add_test(NAME MemPool COMMAND test_mempool)

add_executable(test_lmptype test_lmptype.cpp)
target_link_libraries(test_lmptype PRIVATE lammps GTest::GMockMain)
add_test(NAME LmpType COMMAND test_lmptype)

# Test nbody bigint type capacity (>INT_MAX values)
add_executable(test_nbody_bigint test_nbody_bigint.cpp)
target_link_libraries(test_nbody_bigint PRIVATE lammps GTest::GMockMain)
add_test(NAME NbodyBigint COMMAND test_nbody_bigint)

# Test nbody MPI communication with large values (requires MPI)
if(BUILD_MPI)
  add_executable(test_nbody_mpi test_nbody_mpi.cpp)
  target_link_libraries(test_nbody_mpi PRIVATE lammps GTest::GMockMain MPI::MPI_CXX)
  add_mpi_test(NAME NbodyMPI NUM_PROCS 4 COMMAND $<TARGET_FILE:test_nbody_mpi>)
  set_tests_properties(NbodyMPI PROPERTIES LABELS "mpi")
endif()

add_executable(test_argutils test_argutils.cpp)
target_link_libraries(test_argutils PRIVATE lammps GTest::GMockMain)
add_test(NAME ArgUtils COMMAND test_argutils)

add_executable(test_utils test_utils.cpp)
target_link_libraries(test_utils PRIVATE lammps GTest::GMockMain)
add_test(NAME Utils COMMAND test_utils)

add_executable(test_platform test_platform.cpp)
target_link_libraries(test_platform PRIVATE lammps GTest::GMockMain)
add_test(NAME Platform COMMAND test_platform)

if(PKG_LEPTON)
  add_executable(test_lepton test_lepton.cpp)
  target_link_libraries(test_lepton PRIVATE lepton lammps GTest::GMockMain)
  add_test(NAME Lepton COMMAND test_lepton)
endif()

set_tests_properties(Utils Platform PROPERTIES
  ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR}")

if(BUILD_SHARED_LIBS)
  enable_language(C)
  add_library(testsharedobj MODULE testshared.c)
  set_target_properties(testsharedobj PROPERTIES PREFIX "" WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
  add_library(testsharedlib SHARED testshared.c)
  set_target_properties(testsharedlib PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
  add_dependencies(test_platform testsharedobj testsharedlib)
  target_compile_definitions(test_platform PRIVATE TEST_SHARED_LOAD=1
    TEST_SHARED_LIB=$<TARGET_FILE:testsharedlib>
    TEST_SHARED_OBJ=$<TARGET_FILE:testsharedobj>)
endif()

add_executable(test_fmtlib test_fmtlib.cpp)
target_link_libraries(test_fmtlib PRIVATE lammps GTest::GMockMain)
add_test(NAME FmtLib COMMAND test_fmtlib)

add_executable(test_json test_json.cpp)
target_link_libraries(test_json PRIVATE lammps GTest::GMockMain)
add_test(NAME JSON COMMAND test_json)

add_executable(test_math_eigen_impl test_math_eigen_impl.cpp)
target_include_directories(test_math_eigen_impl PRIVATE ${LAMMPS_SOURCE_DIR})
add_test(NAME MathEigen COMMAND test_math_eigen_impl 10 5)

# =============================================================================
# FFT Testing Infrastructure
# =============================================================================

# -----------------------------------------------------------------------------
# FFT3D Tests (Standard KSPACE FFT)
# -----------------------------------------------------------------------------
if(PKG_KSPACE)
  # Build test_fft3d executable
  add_executable(test_fft3d test_fft3d.cpp ../error_stats.cpp)
  target_link_libraries(test_fft3d PRIVATE lammps GTest::GMockMain)
  target_include_directories(test_fft3d PRIVATE ${LAMMPS_SOURCE_DIR})

  # Add Kokkos if enabled
  if(PKG_KOKKOS)
    target_link_libraries(test_fft3d PRIVATE Kokkos::kokkos)
  endif()

  # Propagate FFT_SINGLE preprocessor definition
  if(FFT_SINGLE)
    target_compile_definitions(test_fft3d PRIVATE -DFFT_SINGLE)
  endif()

  # Register serial test
  add_test(NAME FFT3D COMMAND test_fft3d)
  set_tests_properties(FFT3D PROPERTIES
    LABELS "fft;kspace;utils"
    ENVIRONMENT "OMP_PROC_BIND=false")

  # Register MPI parallel tests (if MPI is enabled)
  if(BUILD_MPI)
    # Test with 2 processes
    add_mpi_test(NAME FFT3D_MPI_2proc NUM_PROCS 2 COMMAND $<TARGET_FILE:test_fft3d>)
    set_tests_properties(FFT3D_MPI_2proc PROPERTIES
      LABELS "fft;kspace;mpi;utils"
      ENVIRONMENT "OMP_PROC_BIND=false")

    # Test with 4 processes
    add_mpi_test(NAME FFT3D_MPI_4proc NUM_PROCS 4 COMMAND $<TARGET_FILE:test_fft3d>)
    set_tests_properties(FFT3D_MPI_4proc PROPERTIES
      LABELS "fft;kspace;mpi;utils"
      ENVIRONMENT "OMP_PROC_BIND=false")
  endif()
endif()

# -----------------------------------------------------------------------------
# FFT3D KOKKOS Tests (KOKKOS FFT - CPU and GPU)
# -----------------------------------------------------------------------------
if(PKG_KSPACE AND PKG_KOKKOS)
  # Build test_fft3d_kokkos executable
  add_executable(test_fft3d_kokkos test_fft3d_kokkos.cpp ../error_stats.cpp)
  target_link_libraries(test_fft3d_kokkos PRIVATE lammps GTest::GMockMain Kokkos::kokkos)
  target_include_directories(test_fft3d_kokkos PRIVATE ${LAMMPS_SOURCE_DIR} ${LAMMPS_SOURCE_DIR}/KSPACE)

  # Propagate FFT_SINGLE preprocessor definition
  if(FFT_SINGLE)
    target_compile_definitions(test_fft3d_kokkos PRIVATE -DFFT_SINGLE)
  endif()

  # Compile the test against the same Kokkos FFT backend as the library.  The
  # -DFFT_KOKKOS_<backend> macro is added PRIVATE to the lammps target (see
  # cmake/Modules/Packages/KOKKOS.cmake), so without re-adding it here the
  # KOKKOS FFT headers in this translation unit would default to KISS and
  # disagree with the actual (FFTW3/MKL/NVPL/cuFFT/...) build of liblammps --
  # an ODR mismatch in the FFT3dKokkos class layout.  Also link the matching
  # FFT library so its headers (fftw3.h, mkl_dfti.h, ...) are on the include
  # path; the imported targets carry the include directories.
  if(FFT_KOKKOS)
    target_compile_definitions(test_fft3d_kokkos PRIVATE -DFFT_KOKKOS_${FFT_KOKKOS})
    if(FFT_KOKKOS STREQUAL "FFTW3")
      if(FFT_SINGLE)
        set(_kokkos_fft_tgt FFTW3F::FFTW3F)
      else()
        set(_kokkos_fft_tgt FFTW3::FFTW3)
      endif()
    elseif(FFT_KOKKOS STREQUAL "NVPL")
      set(_kokkos_fft_tgt nvpl::fftw)
    elseif(FFT_KOKKOS STREQUAL "MKL" OR FFT_KOKKOS STREQUAL "MKL_GPU")
      set(_kokkos_fft_tgt MKL::MKL)
    elseif(FFT_KOKKOS STREQUAL "CUFFT")
      set(_kokkos_fft_tgt CUDA::cufft)
    elseif(FFT_KOKKOS STREQUAL "HIPFFT")
      set(_kokkos_fft_tgt hip::hipfft)
    else()
      set(_kokkos_fft_tgt "")
    endif()
    if(_kokkos_fft_tgt AND TARGET ${_kokkos_fft_tgt})
      target_link_libraries(test_fft3d_kokkos PRIVATE ${_kokkos_fft_tgt})
    endif()
    if((FFT_KOKKOS STREQUAL "FFTW3" OR FFT_KOKKOS STREQUAL "NVPL") AND FFT_FFTW_THREADS)
      target_compile_definitions(test_fft3d_kokkos PRIVATE -DFFT_KOKKOS_FFTW_THREADS)
    endif()
  endif()

  # Register serial test
  add_test(NAME FFT3DKokkos COMMAND test_fft3d_kokkos)
  set_tests_properties(FFT3DKokkos PROPERTIES
    LABELS "fft;kokkos;kspace;utils"
    ENVIRONMENT "OMP_PROC_BIND=false")

  # Register MPI parallel tests (if MPI is enabled)
  if(BUILD_MPI)
    # Test with 2 processes
    add_mpi_test(NAME FFT3DKokkos_MPI_2proc NUM_PROCS 2 COMMAND $<TARGET_FILE:test_fft3d_kokkos>)
    set_tests_properties(FFT3DKokkos_MPI_2proc PROPERTIES
      LABELS "fft;kokkos;kspace;mpi;utils"
      ENVIRONMENT "OMP_PROC_BIND=false")

    # Test with 4 processes
    add_mpi_test(NAME FFT3DKokkos_MPI_4proc NUM_PROCS 4 COMMAND $<TARGET_FILE:test_fft3d_kokkos>)
    set_tests_properties(FFT3DKokkos_MPI_4proc PROPERTIES
      LABELS "fft;kokkos;kspace;mpi;utils"
      ENVIRONMENT "OMP_PROC_BIND=false")
  endif()
endif()

# =============================================================================
# End of FFT Testing Infrastructure
# =============================================================================

# Extra math tests for superellipsoids
if(PKG_ASPHERE)
  add_executable(test_math_extra_superellipsoids test_math_extra_superellipsoids.cpp)
  target_link_libraries(test_math_extra_superellipsoids PRIVATE lammps GTest::GMockMain)
  add_test(NAME MathExtraSuperellipsoids COMMAND test_math_extra_superellipsoids)
endif()

add_executable(test_math_extra test_math_extra.cpp)
target_link_libraries(test_math_extra PRIVATE lammps GTest::GMockMain)
add_test(NAME MathExtra COMMAND test_math_extra)

add_executable(test_math_special test_math_special.cpp)
target_link_libraries(test_math_special PRIVATE lammps GTest::GMockMain)
add_test(NAME MathSpecial COMMAND test_math_special)

add_executable(test_random test_random.cpp)
target_link_libraries(test_random PRIVATE lammps GTest::GMockMain)
add_test(NAME Random COMMAND test_random)

add_executable(test_tabular_function test_tabular_function.cpp)
target_link_libraries(test_tabular_function PRIVATE lammps GTest::GMockMain)
add_test(NAME TabularFunction COMMAND test_tabular_function)
