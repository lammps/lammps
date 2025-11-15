# -*- CMake -*- file for tests of utily functions and classes in LAMMPS

add_executable(test_tokenizer test_tokenizer.cpp)
target_link_libraries(test_tokenizer PRIVATE lammps GTest::GMockMain)
add_test(NAME Tokenizer COMMAND test_tokenizer)

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
  add_executable(test_fft3d test_fft3d.cpp ../force-styles/error_stats.cpp)
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
  add_executable(test_fft3d_kokkos test_fft3d_kokkos.cpp ../force-styles/error_stats.cpp)
  target_link_libraries(test_fft3d_kokkos PRIVATE lammps GTest::GMockMain Kokkos::kokkos)
  target_include_directories(test_fft3d_kokkos PRIVATE ${LAMMPS_SOURCE_DIR} ${LAMMPS_SOURCE_DIR}/KSPACE)

  # Propagate FFT_SINGLE preprocessor definition
  if(FFT_SINGLE)
    target_compile_definitions(test_fft3d_kokkos PRIVATE -DFFT_SINGLE)
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
