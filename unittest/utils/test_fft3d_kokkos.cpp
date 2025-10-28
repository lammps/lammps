/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// Test suite for KOKKOS FFT3d wrapper
//
// Supported Backends (CPU only):
//   - Kokkos::Serial backend (always available)
//   - Kokkos::OpenMP backend (if enabled)
//   - Kokkos::Threads backend (if enabled)
//
// GPU Backend Limitations:
//   Currently, KOKKOS GPU backends (CUDA, HIP, SYCL) are NOT tested because:
//   1. When a GPU backend is configured, it MUST be used (cannot test threading only)
//   2. No reliable, non-crashing way exists to detect viable GPU hardware
//   Therefore, tests are skipped for KOKKOS/CUDA, KOKKOS/HIP, and KOKKOS/SYCL
//
// Test Categories:
//   1. Backend detection and configuration
//   2. Round-trip tests (forward + backward FFT)
//   3. Known answer tests (delta function, sine wave)
//   4. Multiple grid sizes

#include "lmpfftsettings.h"

#ifdef LMP_KOKKOS

#include "KOKKOS/fft3d_kokkos.h"
#include "info.h"
#include "lammps.h"

#include "../testing/core.h"
#include "../utils/fft_test_helpers.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <Kokkos_Core.hpp>
#include <cmath>
#include <cstring>
#include <mpi.h>
#include <random>
#include <vector>

using namespace LAMMPS_NS;
using namespace FFTTestHelpers;
using namespace FFTTestData;
using namespace FFTValidation;

// Verbose output control
bool verbose = false;

// Helper function to check if KOKKOS is using a GPU backend
static bool is_kokkos_gpu_backend()
{
    // Check for GPU backends in priority order
    if (Info::has_accelerator_feature("KOKKOS", "api", "cuda")) return true;
    if (Info::has_accelerator_feature("KOKKOS", "api", "hip")) return true;
    if (Info::has_accelerator_feature("KOKKOS", "api", "sycl")) return true;
    return false;
}

// =============================================================================
// Test Fixture for KOKKOS FFT Tests
// =============================================================================

class FFT3DKokkosTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        // Check if KOKKOS package is available
        if (!Info::has_package("KOKKOS")) {
            GTEST_SKIP() << "Test requires KOKKOS package";
        }

        // Skip GPU backends (no safe way to detect GPU hardware)
        if (is_kokkos_gpu_backend()) {
            GTEST_SKIP() << "KOKKOS GPU backend not testable (no safe GPU detection)";
        }

        // Add -k on to enable KOKKOS in LAMMPS (this creates lmp->kokkos)
        args = {"-log", "none", "-echo", "screen", "-nocite", "-k", "on"};

        // Initialize LAMMPS test
        LAMMPSTest::SetUp();

        // Initialize Kokkos if not already initialized
        if (!Kokkos::is_initialized()) {
            int argc = 0;
            char **argv = nullptr;
            Kokkos::initialize(argc, argv);
            kokkos_initialized_here = true;
        } else {
            kokkos_initialized_here = false;
        }

        // Initialize FFT-related members
        fft = nullptr;
        nfast = 0;
        nmid = 0;
        nslow = 0;
    }

    void TearDown() override
    {
        // Clean up FFT object (note: fft is void* and must not be deleted without cast)
        // The specific test methods will handle deletion with proper casting
        fft = nullptr;

        // Clean up LAMMPS instance
        LAMMPSTest::TearDown();

        // NOTE: Do NOT finalize Kokkos between tests, as FFT3dKokkos constructor
        // may call initialize() again. Kokkos will be finalized automatically at program exit.
    }

    // Helper: Create serial KOKKOS FFT3d object (no MPI decomposition)
    template<typename DeviceType>
    void create_serial_fft(int nfast_in, int nmid_in, int nslow_in)
    {
        nfast = nfast_in;
        nmid = nmid_in;
        nslow = nslow_in;

        // Serial FFT: entire grid on one processor
        int in_ilo = 0, in_ihi = nfast - 1;
        int in_jlo = 0, in_jhi = nmid - 1;
        int in_klo = 0, in_khi = nslow - 1;

        int out_ilo = 0, out_ihi = nfast - 1;
        int out_jlo = 0, out_jhi = nmid - 1;
        int out_klo = 0, out_khi = nslow - 1;

        // FFT parameters
        int scaled = 0;         // No scaling
        int permute = 0;        // No permutation
        int nbuf = 0;           // Buffer size (output)
        int usecollective = 0;  // Use point-to-point communication
        int usegpu = 0;         // Let KOKKOS decide based on backend

        // Create FFT3dKokkos object
        BEGIN_HIDE_OUTPUT();
        fft = new FFT3dKokkos<DeviceType>(
            lmp, MPI_COMM_WORLD, nfast, nmid, nslow, in_ilo, in_ihi, in_jlo, in_jhi, in_klo,
            in_khi, out_ilo, out_ihi, out_jlo, out_jhi, out_klo, out_khi, scaled, permute, &nbuf,
            usecollective, usegpu);
        END_HIDE_OUTPUT();

        ASSERT_NE(fft, nullptr);
    }

    // Helper: Perform round-trip test (forward + backward FFT)
    template<typename DeviceType>
    void run_roundtrip_test(int nfast_in, int nmid_in, int nslow_in)
    {
        create_serial_fft<DeviceType>(nfast_in, nmid_in, nslow_in);

        int nsize = nfast * nmid * nslow;

        // Create Kokkos views for input and output data on device
        typedef FFTArrayTypes<DeviceType> FFT_AT;
        typename FFT_AT::t_FFT_SCALAR_1d d_input("fft_input", 2 * nsize);
        typename FFT_AT::t_FFT_SCALAR_1d d_output("fft_output", 2 * nsize);

        // Create host mirror for initialization
        auto h_input = Kokkos::create_mirror_view(d_input);

        // Generate random complex data on host
        RandomComplexGenerator gen(42424);
        std::vector<FFT_SCALAR> temp_data(2 * nsize);
        gen.generate(temp_data.data(), nfast, nmid, nslow);
        for (int i = 0; i < 2 * nsize; ++i) {
            h_input(i) = temp_data[i];
        }

        // Copy input data to device
        Kokkos::deep_copy(d_input, h_input);

        // Perform forward FFT
        auto fft_ptr = static_cast<FFT3dKokkos<DeviceType>*>(fft);
        fft_ptr->compute(d_input, d_output, FFT3dKokkos<DeviceType>::FORWARD);

        // Perform backward FFT (in-place: output becomes input)
        fft_ptr->compute(d_output, d_output, FFT3dKokkos<DeviceType>::BACKWARD);

        // Copy result back to host for validation
        auto h_output = Kokkos::create_mirror_view(d_output);
        Kokkos::deep_copy(h_output, d_output);

        // Convert Kokkos views to raw pointers for validator
        std::vector<FFT_SCALAR> input_vec(2 * nsize);
        std::vector<FFT_SCALAR> output_vec(2 * nsize);

        for (int i = 0; i < 2 * nsize; ++i) {
            input_vec[i] = h_input(i);
            output_vec[i] = h_output(i) / static_cast<FFT_SCALAR>(nsize);  // Apply normalization
        }

        // Validate round-trip: output should equal input after normalization
        RoundTripValidator validator(input_vec.data(), output_vec.data(), nfast, nmid, nslow,
                                      ROUNDTRIP_TOLERANCE, verbose);
        bool valid = validator.validate();

        if (verbose || !valid) {
            std::cout << "Round-trip test:" << std::endl;
            std::cout << "  Max error: " << validator.get_error_stats().max() << std::endl;
            std::cout << "  Avg error: " << validator.get_error_stats().avg() << std::endl;
            std::cout << "  Status: " << (valid ? "PASSED" : "FAILED") << std::endl;
        }

        // Clean up FFT object with proper cast
        delete static_cast<FFT3dKokkos<DeviceType>*>(fft);
        fft = nullptr;

        EXPECT_TRUE(valid) << "Round-trip test failed";
    }

    // Helper: Perform known-answer test (delta function)
    template<typename DeviceType>
    void run_delta_test(int nfast_in, int nmid_in, int nslow_in)
    {
        create_serial_fft<DeviceType>(nfast_in, nmid_in, nslow_in);

        int nsize = nfast * nmid * nslow;

        // Create Kokkos views
        typedef FFTArrayTypes<DeviceType> FFT_AT;
        typename FFT_AT::t_FFT_SCALAR_1d d_input("fft_input", 2 * nsize);
        typename FFT_AT::t_FFT_SCALAR_1d d_output("fft_output", 2 * nsize);

        // Create host mirror and generate delta function
        auto h_input = Kokkos::create_mirror_view(d_input);

        DeltaFunctionGenerator gen;
        gen.generate(h_input.data(), nfast, nmid, nslow);

        // Copy to device
        Kokkos::deep_copy(d_input, h_input);

        // Perform forward FFT
        auto fft_ptr = static_cast<FFT3dKokkos<DeviceType>*>(fft);
        fft_ptr->compute(d_input, d_output, FFT3dKokkos<DeviceType>::FORWARD);

        // Copy result back to host
        auto h_output = Kokkos::create_mirror_view(d_output);
        Kokkos::deep_copy(h_output, d_output);

        // Validate: FFT(delta) should be constant (all values = 1+0i)
        // Create expected output: all bins = (1.0, 0.0)
        std::vector<FFT_SCALAR> output_vec(2 * nsize);
        std::vector<FFT_SCALAR> expected_vec(2 * nsize);

        for (int i = 0; i < 2 * nsize; ++i) {
            output_vec[i] = h_output(i);
        }

        for (int i = 0; i < nsize; ++i) {
            set_complex_linear(expected_vec.data(), i, std::complex<FFT_SCALAR>(1.0, 0.0));
        }

        KnownAnswerValidator validator(output_vec.data(), expected_vec.data(), nfast, nmid, nslow,
                                        1e-10, verbose);
        bool valid = validator.validate();

        if (verbose || !valid) {
            std::cout << "Delta function test:" << std::endl;
            std::cout << "  Max error: " << validator.get_error_stats().max() << std::endl;
            std::cout << "  Avg error: " << validator.get_error_stats().avg() << std::endl;
            std::cout << "  Status: " << (valid ? "PASSED" : "FAILED") << std::endl;
        }

        // Clean up FFT object with proper cast
        delete static_cast<FFT3dKokkos<DeviceType>*>(fft);
        fft = nullptr;

        EXPECT_TRUE(valid) << "Delta function test failed";
    }

    // Member variables
    void *fft;  // Type-erased pointer (actual type is FFT3dKokkos<DeviceType>*)
    int nfast, nmid, nslow;
    bool kokkos_initialized_here;
};

// =============================================================================
// Test 1: Backend Detection
// =============================================================================

TEST_F(FFT3DKokkosTest, BackendDetection)
{
    // Verify KOKKOS configuration is accessible via has_accelerator_feature
    EXPECT_TRUE(Info::has_package("KOKKOS")) << "KOKKOS package should be available";

    // Verify at least one backend is enabled
    bool has_backend = Info::has_accelerator_feature("KOKKOS", "api", "cuda") ||
                       Info::has_accelerator_feature("KOKKOS", "api", "hip") ||
                       Info::has_accelerator_feature("KOKKOS", "api", "sycl") ||
                       Info::has_accelerator_feature("KOKKOS", "api", "openmp") ||
                       Info::has_accelerator_feature("KOKKOS", "api", "serial") ||
                       Info::has_accelerator_feature("KOKKOS", "api", "pthreads");

    EXPECT_TRUE(has_backend) << "At least one KOKKOS backend should be enabled";

    // Verify FFT info is accessible
    std::string fft_info = Info::get_fft_info();
    EXPECT_FALSE(fft_info.empty()) << "FFT info should not be empty";
}

// =============================================================================
// cuFFT Tests (NVIDIA CUDA Backend)
// =============================================================================

#if defined(KOKKOS_ENABLE_CUDA) && defined(FFT_KOKKOS_CUFFT)

TEST_F(FFT3DKokkosTest, RoundTrip_cuFFT_32x32x32)
{


    // Use LMPDeviceType which is Kokkos::Cuda when CUDA is enabled
    typedef Kokkos::Cuda DeviceType;
    run_roundtrip_test<DeviceType>(32, 32, 32);
}

TEST_F(FFT3DKokkosTest, RoundTrip_cuFFT_64x64x64)
{


    typedef Kokkos::Cuda DeviceType;
    run_roundtrip_test<DeviceType>(64, 64, 64);
}

TEST_F(FFT3DKokkosTest, KnownAnswer_cuFFT_DeltaFunction)
{


    typedef Kokkos::Cuda DeviceType;
    run_delta_test<DeviceType>(32, 32, 32);
}

#endif  // KOKKOS_ENABLE_CUDA && FFT_KOKKOS_CUFFT

// =============================================================================
// hipFFT Tests (AMD HIP Backend)
// =============================================================================

#if defined(KOKKOS_ENABLE_HIP) && defined(FFT_KOKKOS_HIPFFT)

TEST_F(FFT3DKokkosTest, RoundTrip_hipFFT_32x32x32)
{


    typedef Kokkos::HIP DeviceType;
    run_roundtrip_test<DeviceType>(32, 32, 32);
}

TEST_F(FFT3DKokkosTest, RoundTrip_hipFFT_64x64x64)
{


    typedef Kokkos::HIP DeviceType;
    run_roundtrip_test<DeviceType>(64, 64, 64);
}

TEST_F(FFT3DKokkosTest, KnownAnswer_hipFFT_DeltaFunction)
{


    typedef Kokkos::HIP DeviceType;
    run_delta_test<DeviceType>(32, 32, 32);
}

#endif  // KOKKOS_ENABLE_HIP && FFT_KOKKOS_HIPFFT

// =============================================================================
// MKL_GPU Tests (Intel SYCL Backend)
// =============================================================================

#if defined(KOKKOS_ENABLE_SYCL) && defined(FFT_KOKKOS_MKL_GPU)

TEST_F(FFT3DKokkosTest, RoundTrip_MKL_GPU_32x32x32)
{


    typedef Kokkos::Experimental::SYCL DeviceType;
    run_roundtrip_test<DeviceType>(32, 32, 32);
}

TEST_F(FFT3DKokkosTest, RoundTrip_MKL_GPU_64x64x64)
{


    typedef Kokkos::Experimental::SYCL DeviceType;
    run_roundtrip_test<DeviceType>(64, 64, 64);
}

TEST_F(FFT3DKokkosTest, KnownAnswer_MKL_GPU_DeltaFunction)
{


    typedef Kokkos::Experimental::SYCL DeviceType;
    run_delta_test<DeviceType>(32, 32, 32);
}

#endif  // KOKKOS_ENABLE_SYCL && FFT_KOKKOS_MKL_GPU

// =============================================================================
// CPU Backend Tests (Task 4.3)
// =============================================================================
// These tests use LMPHostType which maps to:
//   - Kokkos::Serial (default, always available)
//   - Kokkos::OpenMP (if KOKKOS_ENABLE_OPENMP)
//   - Kokkos::Threads (if KOKKOS_ENABLE_THREADS)

TEST_F(FFT3DKokkosTest, RoundTrip_Kokkos_Serial_32x32x32)
{

    // LMPHostType is the CPU execution space (Serial, OpenMP, or Threads)
    run_roundtrip_test<LMPHostType>(32, 32, 32);
}

TEST_F(FFT3DKokkosTest, RoundTrip_Kokkos_Serial_64x64x64)
{

    run_roundtrip_test<LMPHostType>(64, 64, 64);
}

#if defined(KOKKOS_ENABLE_OPENMP)
TEST_F(FFT3DKokkosTest, RoundTrip_Kokkos_OpenMP_32x32x32)
{
    // When OpenMP is enabled, LMPHostType == Kokkos::OpenMP
    run_roundtrip_test<LMPHostType>(32, 32, 32);
}
#endif  // KOKKOS_ENABLE_OPENMP

#if defined(KOKKOS_ENABLE_THREADS)
TEST_F(FFT3DKokkosTest, RoundTrip_Kokkos_Threads_32x32x32)
{
    // When Threads is enabled, LMPHostType == Kokkos::Threads
    run_roundtrip_test<LMPHostType>(32, 32, 32);
}
#endif  // KOKKOS_ENABLE_THREADS

TEST_F(FFT3DKokkosTest, KnownAnswer_Kokkos_DeltaFunction)
{

    run_delta_test<LMPHostType>(32, 32, 32);
}

TEST_F(FFT3DKokkosTest, KnownAnswer_Kokkos_Sine)
{

    // Create FFT
    typedef LMPHostType DeviceType;
    create_serial_fft<DeviceType>(32, 32, 32);

    int nsize = nfast * nmid * nslow;

    // Create Kokkos views
    typedef FFTArrayTypes<DeviceType> FFT_AT;
    typename FFT_AT::t_FFT_SCALAR_1d d_input("fft_input", 2 * nsize);
    typename FFT_AT::t_FFT_SCALAR_1d d_output("fft_output", 2 * nsize);

    // Create host mirror and generate sine wave
    auto h_input = Kokkos::create_mirror_view(d_input);

    SineWaveGenerator gen(2, 0, 0, 1.0);
    std::vector<FFT_SCALAR> temp_data(2 * nsize);
    gen.generate(temp_data.data(), nfast, nmid, nslow);

    for (int i = 0; i < 2 * nsize; ++i) {
        h_input(i) = temp_data[i];
    }

    // Copy to device
    Kokkos::deep_copy(d_input, h_input);

    // Perform forward FFT
    auto fft_ptr = static_cast<FFT3dKokkos<DeviceType>*>(fft);
    fft_ptr->compute(d_input, d_output, FFT3dKokkos<DeviceType>::FORWARD);

    // Copy result back to host
    auto h_output = Kokkos::create_mirror_view(d_output);
    Kokkos::deep_copy(h_output, d_output);

    // Convert output to vector
    std::vector<FFT_SCALAR> output_vec(2 * nsize);
    for (int i = 0; i < 2 * nsize; ++i) {
        output_vec[i] = h_output(i);
    }

    // Validate: FFT(sin(k·x)) should have spikes at ±k
    FFT_SCALAR spike_amplitude = static_cast<FFT_SCALAR>(nsize) / 2.0;

    // Create expected output
    std::vector<FFT_SCALAR> expected_data(2 * nsize, 0.0);
    set_complex(expected_data.data(), 2, 0, 0, nfast, nmid,
                std::complex<FFT_SCALAR>(0.0, -spike_amplitude));
    set_complex(expected_data.data(), nfast - 2, 0, 0, nfast, nmid,
                std::complex<FFT_SCALAR>(0.0, spike_amplitude));

    KnownAnswerValidator validator(output_vec.data(), expected_data.data(), nfast, nmid, nslow,
                                    1e-10, verbose);
    bool valid = validator.validate();

    if (verbose || !valid) {
        std::cout << "Sine wave test:" << std::endl;
        std::cout << "  Expected spikes at k=(2,0,0) and k=(" << (nfast - 2) << ",0,0)"
                  << std::endl;
        std::cout << "  Spike amplitude: " << spike_amplitude << std::endl;
        std::cout << "  Max error: " << validator.get_error_stats().max() << std::endl;
        std::cout << "  Status: " << (valid ? "PASSED" : "FAILED") << std::endl;
    }

    EXPECT_TRUE(valid) << "Sine wave test failed";
}

// =============================================================================
// TASK 4.5: Threading Tests (OpenMP, Threads, Thread Safety)
// =============================================================================

#if defined(KOKKOS_ENABLE_OPENMP) || defined(KOKKOS_ENABLE_THREADS)

// ----------------------------------------------------------------------------
// Test: Threading_OpenMP_Concurrent
// Purpose: Test multiple concurrent FFTs with OpenMP backend
// ----------------------------------------------------------------------------

TEST_F(FFT3DKokkosTest, Threading_OpenMP_Concurrent)
{
#if !defined(KOKKOS_ENABLE_OPENMP)
    GTEST_SKIP() << "Test requires KOKKOS OpenMP backend";
#else
    // Grid dimensions
    const int grid_size = 32;
    const int nsize = grid_size * grid_size * grid_size;
    const int num_ffts = 4;


    // Create multiple FFT instances
    std::vector<FFT3dKokkos<LMPDeviceType> *> ffts;
    std::vector<std::vector<FFT_SCALAR>> input_buffers(num_ffts);
    std::vector<std::vector<FFT_SCALAR>> original_buffers(num_ffts);

    // Initialize FFT instances and data
    for (int n = 0; n < num_ffts; n++) {
        input_buffers[n].resize(2 * nsize);
        original_buffers[n].resize(2 * nsize);

        // Generate random data (unique seed for each FFT)
        RandomComplexGenerator gen(12345 + n);
        gen.generate(input_buffers[n].data(), grid_size, grid_size, grid_size);
        gen.generate(original_buffers[n].data(), grid_size, grid_size, grid_size);

        // Create FFT instance
        int in_ilo = 0, in_ihi = grid_size - 1;
        int in_jlo = 0, in_jhi = grid_size - 1;
        int in_klo = 0, in_khi = grid_size - 1;
        int out_ilo = in_ilo, out_ihi = in_ihi;
        int out_jlo = in_jlo, out_jhi = in_jhi;
        int out_klo = in_klo, out_khi = in_khi;
        int scaled = 0, permute = 0, nbuf = 0;
        int usecollective = 0, usegpu_aware = 0;

        BEGIN_HIDE_OUTPUT();
        auto fft = new FFT3dKokkos<LMPDeviceType>(
            lmp, MPI_COMM_WORLD, grid_size, grid_size, grid_size, in_ilo, in_ihi, in_jlo, in_jhi,
            in_klo, in_khi, out_ilo, out_ihi, out_jlo, out_jhi, out_klo, out_khi, scaled, permute,
            &nbuf, usecollective, usegpu_aware);
        END_HIDE_OUTPUT();

        ffts.push_back(fft);
    }

    // Run FFTs and validate
    bool all_passed = true;

    for (int n = 0; n < num_ffts; n++) {
        typedef FFTArrayTypes<LMPDeviceType> FFT_AT;
        typename FFT_AT::t_FFT_SCALAR_1d d_in("fft_input", 2 * nsize);
        typename FFT_AT::t_FFT_SCALAR_1d d_out("fft_output", 2 * nsize);

        // Copy input data to device
        auto h_in = Kokkos::create_mirror_view(d_in);
        for (int i = 0; i < 2 * nsize; i++) {
            h_in(i) = input_buffers[n][i];
        }
        Kokkos::deep_copy(d_in, h_in);

        // Forward FFT
        BEGIN_HIDE_OUTPUT();
        ffts[n]->compute(d_in, d_out, FFT3dKokkos<LMPDeviceType>::FORWARD);
        END_HIDE_OUTPUT();

        // Backward FFT
        BEGIN_HIDE_OUTPUT();
        ffts[n]->compute(d_out, d_in, FFT3dKokkos<LMPDeviceType>::BACKWARD);
        END_HIDE_OUTPUT();

        // Copy result back to host
        Kokkos::deep_copy(h_in, d_in);

        // Apply normalization
        FFT_SCALAR norm = 1.0 / static_cast<FFT_SCALAR>(nsize);
        for (int i = 0; i < 2 * nsize; i++) {
            input_buffers[n][i] = h_in(i) * norm;
        }

        // Validate round-trip
        RoundTripValidator validator(original_buffers[n].data(), input_buffers[n].data(),
                                      grid_size, grid_size, grid_size, ROUNDTRIP_TOLERANCE, verbose);
        if (!validator.validate()) {
            all_passed = false;
            std::cout << "  FFT instance " << n << " FAILED" << std::endl;
        }
    }

    // Clean up
    for (auto fft : ffts) {
        delete fft;
    }

    EXPECT_TRUE(all_passed) << "One or more concurrent FFTs failed validation";
#endif
}

// ----------------------------------------------------------------------------
// Test: Threading_Threads_Concurrent
// Purpose: Test multiple concurrent FFTs with Threads backend
// ----------------------------------------------------------------------------

TEST_F(FFT3DKokkosTest, Threading_Threads_Concurrent)
{
#if !defined(KOKKOS_ENABLE_THREADS)
    GTEST_SKIP() << "Test requires KOKKOS Threads backend";
#else
    // Grid dimensions
    const int grid_size = 32;
    const int nsize = grid_size * grid_size * grid_size;
    const int num_ffts = 4;


    // Create multiple FFT instances
    std::vector<FFT3dKokkos<LMPDeviceType> *> ffts;
    std::vector<std::vector<FFT_SCALAR>> input_buffers(num_ffts);
    std::vector<std::vector<FFT_SCALAR>> original_buffers(num_ffts);

    // Initialize FFT instances and data
    for (int n = 0; n < num_ffts; n++) {
        input_buffers[n].resize(2 * nsize);
        original_buffers[n].resize(2 * nsize);

        // Generate random data
        RandomComplexGenerator gen(54321 + n);
        gen.generate(input_buffers[n].data(), grid_size, grid_size, grid_size);
        gen.generate(original_buffers[n].data(), grid_size, grid_size, grid_size);

        // Create FFT instance
        int in_ilo = 0, in_ihi = grid_size - 1;
        int in_jlo = 0, in_jhi = grid_size - 1;
        int in_klo = 0, in_khi = grid_size - 1;
        int out_ilo = in_ilo, out_ihi = in_ihi;
        int out_jlo = in_jlo, out_jhi = in_jhi;
        int out_klo = in_klo, out_khi = in_khi;
        int scaled = 0, permute = 0, nbuf = 0;
        int usecollective = 0, usegpu_aware = 0;

        BEGIN_HIDE_OUTPUT();
        auto fft = new FFT3dKokkos<LMPDeviceType>(
            lmp, MPI_COMM_WORLD, grid_size, grid_size, grid_size, in_ilo, in_ihi, in_jlo, in_jhi,
            in_klo, in_khi, out_ilo, out_ihi, out_jlo, out_jhi, out_klo, out_khi, scaled, permute,
            &nbuf, usecollective, usegpu_aware);
        END_HIDE_OUTPUT();

        ffts.push_back(fft);
    }

    // Run FFTs and validate
    bool all_passed = true;

    for (int n = 0; n < num_ffts; n++) {
        typedef FFTArrayTypes<LMPDeviceType> FFT_AT;
        typename FFT_AT::t_FFT_SCALAR_1d d_in("fft_input", 2 * nsize);
        typename FFT_AT::t_FFT_SCALAR_1d d_out("fft_output", 2 * nsize);

        auto h_in = Kokkos::create_mirror_view(d_in);
        for (int i = 0; i < 2 * nsize; i++) {
            h_in(i) = input_buffers[n][i];
        }
        Kokkos::deep_copy(d_in, h_in);

        BEGIN_HIDE_OUTPUT();
        ffts[n]->compute(d_in, d_out, FFT3dKokkos<LMPDeviceType>::FORWARD);
        ffts[n]->compute(d_out, d_in, FFT3dKokkos<LMPDeviceType>::BACKWARD);
        END_HIDE_OUTPUT();

        Kokkos::deep_copy(h_in, d_in);

        FFT_SCALAR norm = 1.0 / static_cast<FFT_SCALAR>(nsize);
        for (int i = 0; i < 2 * nsize; i++) {
            input_buffers[n][i] = h_in(i) * norm;
        }

        RoundTripValidator validator(original_buffers[n].data(), input_buffers[n].data(),
                                      grid_size, grid_size, grid_size, ROUNDTRIP_TOLERANCE, verbose);
        if (!validator.validate()) {
            all_passed = false;
        }
    }

    for (auto fft : ffts) {
        delete fft;
    }

    EXPECT_TRUE(all_passed) << "One or more concurrent FFTs failed validation";
#endif
}

// ----------------------------------------------------------------------------
// Test: Threading_Safety
// Purpose: Validate thread-safety of FFT operations
// ----------------------------------------------------------------------------

TEST_F(FFT3DKokkosTest, Threading_Safety)
{
    // Test runs with any CPU backend (OpenMP, Threads, Serial)

    const int grid_size = 32;
    const int nsize = grid_size * grid_size * grid_size;

    // Create FFT instance
    int in_ilo = 0, in_ihi = grid_size - 1;
    int in_jlo = 0, in_jhi = grid_size - 1;
    int in_klo = 0, in_khi = grid_size - 1;
    int out_ilo = in_ilo, out_ihi = in_ihi;
    int out_jlo = in_jlo, out_jhi = in_jhi;
    int out_klo = in_klo, out_khi = in_khi;
    int scaled = 0, permute = 0, nbuf = 0;
    int usecollective = 0, usegpu_aware = 0;

    BEGIN_HIDE_OUTPUT();
    auto fft_device = new FFT3dKokkos<LMPDeviceType>(
        lmp, MPI_COMM_WORLD, grid_size, grid_size, grid_size, in_ilo, in_ihi, in_jlo, in_jhi,
        in_klo, in_khi, out_ilo, out_ihi, out_jlo, out_jhi, out_klo, out_khi, scaled, permute,
        &nbuf, usecollective, usegpu_aware);
    END_HIDE_OUTPUT();

    ASSERT_NE(fft_device, nullptr);

    // Generate test data
    std::vector<FFT_SCALAR> input_data(2 * nsize);
    std::vector<FFT_SCALAR> original_data(2 * nsize);

    RandomComplexGenerator gen(98765);
    gen.generate(input_data.data(), grid_size, grid_size, grid_size);
    gen.generate(original_data.data(), grid_size, grid_size, grid_size);

    // Create Kokkos views
    typedef FFTArrayTypes<LMPDeviceType> FFT_AT;
    typename FFT_AT::t_FFT_SCALAR_1d d_in("fft_input", 2 * nsize);
    typename FFT_AT::t_FFT_SCALAR_1d d_out("fft_output", 2 * nsize);

    // Copy input data
    auto h_in = Kokkos::create_mirror_view(d_in);
    for (int i = 0; i < 2 * nsize; i++) {
        h_in(i) = input_data[i];
    }
    Kokkos::deep_copy(d_in, h_in);

    // Run multiple FFT operations to check for data corruption
    const int num_iterations = 10;
    bool all_passed = true;

    for (int iter = 0; iter < num_iterations; iter++) {
        BEGIN_HIDE_OUTPUT();
        fft_device->compute(d_in, d_out, FFT3dKokkos<LMPDeviceType>::FORWARD);
        fft_device->compute(d_out, d_in, FFT3dKokkos<LMPDeviceType>::BACKWARD);
        END_HIDE_OUTPUT();

        Kokkos::deep_copy(h_in, d_in);

        FFT_SCALAR norm = 1.0 / static_cast<FFT_SCALAR>(nsize);
        for (int i = 0; i < 2 * nsize; i++) {
            input_data[i] = h_in(i) * norm;
        }

        RoundTripValidator validator(original_data.data(), input_data.data(),
                                      grid_size, grid_size, grid_size, ROUNDTRIP_TOLERANCE, verbose);
        if (!validator.validate()) {
            all_passed = false;
            std::cout << "  Iteration " << iter << " FAILED" << std::endl;
            break;
        }

        // Reset input for next iteration
        for (int i = 0; i < 2 * nsize; i++) {
            h_in(i) = original_data[i];
        }
        Kokkos::deep_copy(d_in, h_in);
    }

    delete fft_device;

    EXPECT_TRUE(all_passed) << "Thread safety validation failed";
}

#endif  // KOKKOS_ENABLE_OPENMP || KOKKOS_ENABLE_THREADS

// =============================================================================
// TASK 4.6: MPI Tests (2 procs, 4 procs, GPU+MPI)
// =============================================================================

// ----------------------------------------------------------------------------
// Test: RoundTrip_Kokkos_MPI_2proc_32x32x32
// Purpose: Test MPI parallel FFT with 2 processes
// Pattern: Based on test_fft3d.cpp lines 615-757
// ----------------------------------------------------------------------------

TEST_F(FFT3DKokkosTest, RoundTrip_Kokkos_MPI_2proc_32x32x32)
{
    // Check MPI environment
    int nprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (nprocs != 2) {
        GTEST_SKIP() << "Test requires exactly 2 MPI processes, got " << nprocs;
    }

    if (rank == 0) {
    }

    // Grid dimensions
    const int grid_size = 32;
    const int nsize_global = grid_size * grid_size * grid_size;

    // Domain decomposition: split along slow (z) dimension
    // Process 0: z = 0..15
    // Process 1: z = 16..31
    int in_ilo = 0, in_ihi = grid_size - 1;
    int in_jlo = 0, in_jhi = grid_size - 1;
    int in_klo, in_khi;

    if (rank == 0) {
        in_klo = 0;
        in_khi = grid_size / 2 - 1;  // 0..15
    } else {
        in_klo = grid_size / 2;
        in_khi = grid_size - 1;  // 16..31
    }

    int out_ilo = in_ilo, out_ihi = in_ihi;
    int out_jlo = in_jlo, out_jhi = in_jhi;
    int out_klo = in_klo, out_khi = in_khi;

    // Calculate local size
    int local_nslow = in_khi - in_klo + 1;
    int local_size = grid_size * grid_size * local_nslow;

    // Allocate local data buffers
    std::vector<FFT_SCALAR> input_data(2 * local_size);
    std::vector<FFT_SCALAR> original_data(2 * local_size);

    // Generate random complex data using global coordinates for consistency
    RandomComplexGenerator base_gen(12345);
    for (int k = 0; k < local_nslow; k++) {
        for (int j = 0; j < grid_size; j++) {
            for (int i = 0; i < grid_size; i++) {
                int global_k = in_klo + k;
                int global_idx = global_k * grid_size * grid_size + j * grid_size + i;

                // Generate deterministic random data for this grid point
                std::mt19937 rng(12345 + global_idx);
                std::uniform_real_distribution<FFT_SCALAR> dist(-1.0, 1.0);
                FFT_SCALAR re = dist(rng);
                FFT_SCALAR im = dist(rng);

                int local_idx = k * grid_size * grid_size + j * grid_size + i;
                input_data[2 * local_idx] = re;
                input_data[2 * local_idx + 1] = im;
                original_data[2 * local_idx] = re;
                original_data[2 * local_idx + 1] = im;
            }
        }
    }

    // FFT parameters
    int scaled = 0, permute = 0, nbuf = 0;
    int usecollective = 0, usegpu_aware = 0;

    // Create MPI-aware FFT3dKokkos object
    BEGIN_HIDE_OUTPUT();
    auto fft_mpi = new FFT3dKokkos<LMPDeviceType>(
        lmp, MPI_COMM_WORLD, grid_size, grid_size, grid_size, in_ilo, in_ihi, in_jlo, in_jhi,
        in_klo, in_khi, out_ilo, out_ihi, out_jlo, out_jhi, out_klo, out_khi, scaled, permute,
        &nbuf, usecollective, usegpu_aware);
    END_HIDE_OUTPUT();

    ASSERT_NE(fft_mpi, nullptr);

    // Create Kokkos views
    typedef FFTArrayTypes<LMPDeviceType> FFT_AT;
    typename FFT_AT::t_FFT_SCALAR_1d d_in("fft_input", 2 * local_size);
    typename FFT_AT::t_FFT_SCALAR_1d d_out("fft_output", 2 * local_size);

    // Copy input data to device
    auto h_in = Kokkos::create_mirror_view(d_in);
    for (int i = 0; i < 2 * local_size; i++) {
        h_in(i) = input_data[i];
    }
    Kokkos::deep_copy(d_in, h_in);

    // Forward FFT
    BEGIN_HIDE_OUTPUT();
    fft_mpi->compute(d_in, d_out, FFT3dKokkos<LMPDeviceType>::FORWARD);
    END_HIDE_OUTPUT();

    // Backward FFT
    BEGIN_HIDE_OUTPUT();
    fft_mpi->compute(d_out, d_in, FFT3dKokkos<LMPDeviceType>::BACKWARD);
    END_HIDE_OUTPUT();

    // Copy result back
    Kokkos::deep_copy(h_in, d_in);

    // Apply normalization
    FFT_SCALAR norm = 1.0 / static_cast<FFT_SCALAR>(nsize_global);
    for (int i = 0; i < 2 * local_size; i++) {
        input_data[i] = h_in(i) * norm;
    }

    // Validate round-trip on local data
    RoundTripValidator validator(original_data.data(), input_data.data(),
                                  grid_size, grid_size, local_nslow, ROUNDTRIP_TOLERANCE, verbose);
    bool passed = validator.validate();

    // Gather validation results from all ranks
    int passed_int = passed ? 1 : 0;
    int all_passed;
    MPI_Allreduce(&passed_int, &all_passed, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    // Gather max error from all ranks
    double local_max_error = validator.get_error_stats().max();
    double global_max_error;
    MPI_Allreduce(&local_max_error, &global_max_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    delete fft_mpi;

    EXPECT_TRUE(all_passed) << "Round-trip validation failed on rank " << rank;
}

// ----------------------------------------------------------------------------
// Test: RoundTrip_Kokkos_MPI_4proc_64x64x64
// Purpose: Test MPI parallel FFT with 4 processes
// Pattern: Based on test_fft3d.cpp lines 763-900
// ----------------------------------------------------------------------------

TEST_F(FFT3DKokkosTest, RoundTrip_Kokkos_MPI_4proc_64x64x64)
{
    // Check MPI environment
    int nprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (nprocs != 4) {
        GTEST_SKIP() << "Test requires exactly 4 MPI processes, got " << nprocs;
    }

    if (rank == 0) {
    }

    // Grid dimensions (larger for better load balancing)
    const int grid_size = 64;
    const int nsize_global = grid_size * grid_size * grid_size;

    // Domain decomposition: split along slow (z) dimension
    int in_ilo = 0, in_ihi = grid_size - 1;
    int in_jlo = 0, in_jhi = grid_size - 1;
    int in_klo, in_khi;

    int slices_per_proc = grid_size / nprocs;
    in_klo = rank * slices_per_proc;
    in_khi = (rank + 1) * slices_per_proc - 1;

    int out_ilo = in_ilo, out_ihi = in_ihi;
    int out_jlo = in_jlo, out_jhi = in_jhi;
    int out_klo = in_klo, out_khi = in_khi;

    // Calculate local size
    int local_nslow = in_khi - in_klo + 1;
    int local_size = grid_size * grid_size * local_nslow;

    // Allocate local data buffers
    std::vector<FFT_SCALAR> input_data(2 * local_size);
    std::vector<FFT_SCALAR> original_data(2 * local_size);

    // Generate random complex data
    for (int k = 0; k < local_nslow; k++) {
        for (int j = 0; j < grid_size; j++) {
            for (int i = 0; i < grid_size; i++) {
                int global_k = in_klo + k;
                int global_idx = global_k * grid_size * grid_size + j * grid_size + i;

                // Generate deterministic random data for this grid point
                std::mt19937 rng(54321 + global_idx);
                std::uniform_real_distribution<FFT_SCALAR> dist(-1.0, 1.0);
                FFT_SCALAR re = dist(rng);
                FFT_SCALAR im = dist(rng);

                int local_idx = k * grid_size * grid_size + j * grid_size + i;
                input_data[2 * local_idx] = re;
                input_data[2 * local_idx + 1] = im;
                original_data[2 * local_idx] = re;
                original_data[2 * local_idx + 1] = im;
            }
        }
    }

    // FFT parameters
    int scaled = 0, permute = 0, nbuf = 0;
    int usecollective = 0, usegpu_aware = 0;

    // Create MPI-aware FFT3dKokkos object
    BEGIN_HIDE_OUTPUT();
    auto fft_mpi = new FFT3dKokkos<LMPDeviceType>(
        lmp, MPI_COMM_WORLD, grid_size, grid_size, grid_size, in_ilo, in_ihi, in_jlo, in_jhi,
        in_klo, in_khi, out_ilo, out_ihi, out_jlo, out_jhi, out_klo, out_khi, scaled, permute,
        &nbuf, usecollective, usegpu_aware);
    END_HIDE_OUTPUT();

    ASSERT_NE(fft_mpi, nullptr);

    // Create Kokkos views
    typedef FFTArrayTypes<LMPDeviceType> FFT_AT;
    typename FFT_AT::t_FFT_SCALAR_1d d_in("fft_input", 2 * local_size);
    typename FFT_AT::t_FFT_SCALAR_1d d_out("fft_output", 2 * local_size);

    // Copy input data to device
    auto h_in = Kokkos::create_mirror_view(d_in);
    for (int i = 0; i < 2 * local_size; i++) {
        h_in(i) = input_data[i];
    }
    Kokkos::deep_copy(d_in, h_in);

    // Forward FFT
    BEGIN_HIDE_OUTPUT();
    fft_mpi->compute(d_in, d_out, FFT3dKokkos<LMPDeviceType>::FORWARD);
    END_HIDE_OUTPUT();

    // Backward FFT
    BEGIN_HIDE_OUTPUT();
    fft_mpi->compute(d_out, d_in, FFT3dKokkos<LMPDeviceType>::BACKWARD);
    END_HIDE_OUTPUT();

    // Copy result back
    Kokkos::deep_copy(h_in, d_in);

    // Apply normalization
    FFT_SCALAR norm = 1.0 / static_cast<FFT_SCALAR>(nsize_global);
    for (int i = 0; i < 2 * local_size; i++) {
        input_data[i] = h_in(i) * norm;
    }

    // Validate round-trip on local data
    RoundTripValidator validator(original_data.data(), input_data.data(),
                                  grid_size, grid_size, local_nslow, ROUNDTRIP_TOLERANCE, verbose);
    bool passed = validator.validate();

    // Gather validation results from all ranks
    int passed_int = passed ? 1 : 0;
    int all_passed;
    MPI_Allreduce(&passed_int, &all_passed, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    // Gather max error from all ranks
    double local_max_error = validator.get_error_stats().max();
    double global_max_error;
    MPI_Allreduce(&local_max_error, &global_max_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    delete fft_mpi;

    EXPECT_TRUE(all_passed) << "Round-trip validation failed on rank " << rank;
}

// ----------------------------------------------------------------------------
// Test: RoundTrip_Kokkos_MPI_GPU_2proc
// Purpose: Test MPI + GPU combination (if available)
// ----------------------------------------------------------------------------

TEST_F(FFT3DKokkosTest, RoundTrip_Kokkos_MPI_GPU_2proc)
{
    // Check MPI environment
    int nprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (nprocs != 2) {
        GTEST_SKIP() << "Test requires exactly 2 MPI processes, got " << nprocs;
    }

    // Note: GPU tests are skipped in SetUp() - no safe GPU detection available

    // Grid dimensions
    const int grid_size = 32;
    const int nsize_global = grid_size * grid_size * grid_size;

    // Domain decomposition
    int in_ilo = 0, in_ihi = grid_size - 1;
    int in_jlo = 0, in_jhi = grid_size - 1;
    int in_klo, in_khi;

    if (rank == 0) {
        in_klo = 0;
        in_khi = grid_size / 2 - 1;
    } else {
        in_klo = grid_size / 2;
        in_khi = grid_size - 1;
    }

    int out_ilo = in_ilo, out_ihi = in_ihi;
    int out_jlo = in_jlo, out_jhi = in_jhi;
    int out_klo = in_klo, out_khi = in_khi;

    int local_nslow = in_khi - in_klo + 1;
    int local_size = grid_size * grid_size * local_nslow;

    // Allocate local data buffers
    std::vector<FFT_SCALAR> input_data(2 * local_size);
    std::vector<FFT_SCALAR> original_data(2 * local_size);

    // Generate random complex data
    for (int k = 0; k < local_nslow; k++) {
        for (int j = 0; j < grid_size; j++) {
            for (int i = 0; i < grid_size; i++) {
                int global_k = in_klo + k;
                int global_idx = global_k * grid_size * grid_size + j * grid_size + i;

                // Generate deterministic random data for this grid point
                std::mt19937 rng(99999 + global_idx);
                std::uniform_real_distribution<FFT_SCALAR> dist(-1.0, 1.0);
                FFT_SCALAR re = dist(rng);
                FFT_SCALAR im = dist(rng);

                int local_idx = k * grid_size * grid_size + j * grid_size + i;
                input_data[2 * local_idx] = re;
                input_data[2 * local_idx + 1] = im;
                original_data[2 * local_idx] = re;
                original_data[2 * local_idx + 1] = im;
            }
        }
    }

    // FFT parameters (disable GPU-aware MPI for now)
    int scaled = 0, permute = 0, nbuf = 0;
    int usecollective = 0;
    int usegpu_aware = 0;  // Would check lmp->kokkos->gpu_aware_flag if KokkosLMP was complete

    // Create MPI+GPU FFT3dKokkos object
    BEGIN_HIDE_OUTPUT();
    auto fft_mpi = new FFT3dKokkos<LMPDeviceType>(
        lmp, MPI_COMM_WORLD, grid_size, grid_size, grid_size, in_ilo, in_ihi, in_jlo, in_jhi,
        in_klo, in_khi, out_ilo, out_ihi, out_jlo, out_jhi, out_klo, out_khi, scaled, permute,
        &nbuf, usecollective, usegpu_aware);
    END_HIDE_OUTPUT();

    ASSERT_NE(fft_mpi, nullptr);

    // Create Kokkos views (will be on GPU)
    typedef FFTArrayTypes<LMPDeviceType> FFT_AT;
    typename FFT_AT::t_FFT_SCALAR_1d d_in("fft_input", 2 * local_size);
    typename FFT_AT::t_FFT_SCALAR_1d d_out("fft_output", 2 * local_size);

    // Copy input data to device
    auto h_in = Kokkos::create_mirror_view(d_in);
    for (int i = 0; i < 2 * local_size; i++) {
        h_in(i) = input_data[i];
    }
    Kokkos::deep_copy(d_in, h_in);

    // Forward FFT (on GPU)
    BEGIN_HIDE_OUTPUT();
    fft_mpi->compute(d_in, d_out, FFT3dKokkos<LMPDeviceType>::FORWARD);
    END_HIDE_OUTPUT();

    // Backward FFT (on GPU)
    BEGIN_HIDE_OUTPUT();
    fft_mpi->compute(d_out, d_in, FFT3dKokkos<LMPDeviceType>::BACKWARD);
    END_HIDE_OUTPUT();

    // Copy result back to host
    Kokkos::deep_copy(h_in, d_in);

    // Apply normalization
    FFT_SCALAR norm = 1.0 / static_cast<FFT_SCALAR>(nsize_global);
    for (int i = 0; i < 2 * local_size; i++) {
        input_data[i] = h_in(i) * norm;
    }

    // Validate round-trip
    RoundTripValidator validator(original_data.data(), input_data.data(),
                                  grid_size, grid_size, local_nslow, ROUNDTRIP_TOLERANCE, verbose);
    bool passed = validator.validate();

    // Gather results
    int passed_int = passed ? 1 : 0;
    int all_passed;
    MPI_Allreduce(&passed_int, &all_passed, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    double local_max_error = validator.get_error_stats().max();
    double global_max_error;
    MPI_Allreduce(&local_max_error, &global_max_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    delete fft_mpi;

    EXPECT_TRUE(all_passed) << "Round-trip validation failed on rank " << rank;
}

#endif  // LMP_KOKKOS

// =============================================================================
// Main
// =============================================================================

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    // Check for verbose flag
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose = true;
        }
    }

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
