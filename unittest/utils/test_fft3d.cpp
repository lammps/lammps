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

// Test suite for standard FFT3d wrapper (KISS, FFTW3, MKL, NVPL)
//
// Tests:
//   1. BackendDetection - Verify FFT library and precision detection
//   2. RoundTrip_Serial_32x32x32 - Forward + backward FFT on 32³ grid
//   3. RoundTrip_Serial_64x64x64 - Forward + backward FFT on 64³ grid
//   4. RoundTrip_Serial_48x48x48 - Non-power-of-2 grid (48³)
//   5. KnownAnswer_DeltaFunction - Delta function: FFT(δ) = constant
//   6. KnownAnswer_Constant - Constant field: FFT(1) = N³·δ(k=0)
//   7. KnownAnswer_SineWave - Sine wave: FFT(sin) = spikes at ±k
//   8. ParsevalsTheorem_EnergyConservation - Energy conservation (Parseval's theorem)
//   9. RoundTrip_MPI_2proc_32x32x32 - MPI parallel FFT with 2 processes
//  10. RoundTrip_MPI_4proc_64x64x64 - MPI parallel FFT with 4 processes
//  11. FFTW3_Threading - FFTW3 with threading support (conditional)
//  12. MKL_Optimized - MKL library optimizations (conditional)
//  13. KISS_NonPowerOf2 - KISS FFT with various non-power-of-2 sizes (conditional)
//  14. HeFFTe_Distributed - HeFFTe distributed FFT (conditional)

#include "lmpfftsettings.h"
#include "KSPACE/fft3d_wrap.h"
#include "info.h"
#include "lammps.h"

#include "../testing/core.h"
#include "fft_test_helpers.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <cstring>
#include <mpi.h>
#include <random>
#include <vector>

using namespace LAMMPS_NS;
using namespace FFTTestHelpers;
using namespace FFTTestData;
using namespace FFTValidation;

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

class FFT3DTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        // Check if KSPACE package is available
        if (!Info::has_package("KSPACE")) {
            GTEST_SKIP() << "Test requires KSPACE package";
        }

        // Initialize test parameters
        LAMMPSTest::SetUp();

        // Initialize FFT-related members
        fft = nullptr;
        input_data = nullptr;
        output_data = nullptr;

        // Default grid size (will be set by individual tests)
        nfast = 0;
        nmid = 0;
        nslow = 0;
    }

    void TearDown() override
    {
        // Clean up FFT data
        if (input_data) delete[] input_data;
        if (output_data) delete[] output_data;
        if (fft) delete fft;

        // Clean up LAMMPS instance
        LAMMPSTest::TearDown();
    }

    // Helper: Create serial FFT3d object (no MPI decomposition)
    void create_serial_fft(int nfast_in, int nmid_in, int nslow_in)
    {
        nfast = nfast_in;
        nmid = nmid_in;
        nslow = nslow_in;

        // Total grid size
        int nsize = nfast * nmid * nslow;

        // Allocate data buffers (complex data: 2 * nsize)
        input_data = new FFT_SCALAR[2 * nsize];
        output_data = new FFT_SCALAR[2 * nsize];

        // Zero buffers
        std::memset(input_data, 0, 2 * nsize * sizeof(FFT_SCALAR));
        std::memset(output_data, 0, 2 * nsize * sizeof(FFT_SCALAR));

        // Serial FFT: entire grid on one processor
        int in_ilo = 0, in_ihi = nfast - 1;
        int in_jlo = 0, in_jhi = nmid - 1;
        int in_klo = 0, in_khi = nslow - 1;

        int out_ilo = 0, out_ihi = nfast - 1;
        int out_jlo = 0, out_jhi = nmid - 1;
        int out_klo = 0, out_khi = nslow - 1;

        // FFT parameters
        int scaled = 0;      // No scaling
        int permute = 0;     // No permutation
        int nbuf = 0;        // Buffer size (output)
        int usecollective = 0;  // Use point-to-point communication

        // Create FFT3d object
        BEGIN_HIDE_OUTPUT();
        fft = new FFT3d(lmp, MPI_COMM_WORLD, nfast, nmid, nslow, in_ilo, in_ihi, in_jlo, in_jhi,
                        in_klo, in_khi, out_ilo, out_ihi, out_jlo, out_jhi, out_klo, out_khi,
                        scaled, permute, &nbuf, usecollective);
        END_HIDE_OUTPUT();

        ASSERT_NE(fft, nullptr);
    }

    // Member variables
    FFT3d *fft;
    FFT_SCALAR *input_data;
    FFT_SCALAR *output_data;
    int nfast, nmid, nslow;
};

// ============================================================================
// Test 1: Backend Detection
// ============================================================================

TEST_F(FFT3DTest, BackendDetection)
{
    // Verify FFT configuration is accessible
    std::string fft_info = Info::get_fft_info();
    EXPECT_FALSE(fft_info.empty()) << "FFT info should not be empty";

    // Check FFT library macro is defined
#if defined(FFT_KISS) || defined(FFT_FFTW3) || defined(FFT_MKL) || defined(FFT_NVPL) || defined(FFT_HEFFTE)
    SUCCEED() << "FFT library: " << LMP_FFT_LIB;
#else
    FAIL() << "No FFT library defined";
#endif
}

// ============================================================================
// Test 2: Round-Trip Serial KISS FFT (32x32x32)
// ============================================================================

TEST_F(FFT3DTest, RoundTrip_Serial_32x32x32)
{
    // Create 32x32x32 grid
    create_serial_fft(32, 32, 32);

    int nsize = nfast * nmid * nslow;

    // Generate random complex data using RandomComplexGenerator
    FFTTestData::RandomComplexGenerator generator(12345);
    generator.generate(input_data, nfast, nmid, nslow);

    // Copy input data for later comparison
    std::vector<FFT_SCALAR> original_data(2 * nsize);
    std::copy(input_data, input_data + 2 * nsize, original_data.data());

    // Forward FFT
    BEGIN_HIDE_OUTPUT();
    fft->compute(input_data, output_data, FFT3d::FORWARD);
    END_HIDE_OUTPUT();

    // Backward FFT (into input_data to complete round-trip)
    BEGIN_HIDE_OUTPUT();
    fft->compute(output_data, input_data, FFT3d::BACKWARD);
    END_HIDE_OUTPUT();

    // Apply normalization: LAMMPS backward FFT does not include 1/N³ scaling
    // For round-trip correctness: IFFT(FFT(x)) = N³ × x, so we divide by N³
    FFT_SCALAR norm = 1.0 / static_cast<FFT_SCALAR>(nsize);
    for (int i = 0; i < 2 * nsize; i++) {
        input_data[i] *= norm;
    }

    // Validate round-trip: input_data should match original_data
    FFTValidation::RoundTripValidator validator(original_data.data(), input_data, nfast, nmid, nslow,
                                                 ROUNDTRIP_TOLERANCE, verbose);
    bool passed = validator.validate();

    // Report results
    if (verbose || !passed) {
        std::cout << "Round-trip test (32x32x32):" << std::endl;
        std::cout << "  Max error: " << validator.get_error_stats().max() << " (at index "
                  << validator.get_error_stats().idx() << ")" << std::endl;
        std::cout << "  Avg error: " << validator.get_error_stats().avg() << std::endl;
        std::cout << "  Tolerance: " << ROUNDTRIP_TOLERANCE << std::endl;
        std::cout << "  Status: " << (passed ? "PASSED" : "FAILED") << std::endl;
    }

    EXPECT_TRUE(passed) << "Round-trip validation failed";
    EXPECT_LT(validator.get_error_stats().max(), ROUNDTRIP_TOLERANCE);
}

// ============================================================================
// Test 3: Round-Trip Serial (64x64x64) - Larger Grid
// ============================================================================

TEST_F(FFT3DTest, RoundTrip_Serial_64x64x64)
{
    // Create 64x64x64 grid
    create_serial_fft(64, 64, 64);

    int nsize = nfast * nmid * nslow;

    // Generate random complex data using RandomComplexGenerator
    FFTTestData::RandomComplexGenerator generator(54321);
    generator.generate(input_data, nfast, nmid, nslow);

    // Copy input data for later comparison
    std::vector<FFT_SCALAR> original_data(2 * nsize);
    std::copy(input_data, input_data + 2 * nsize, original_data.data());

    // Forward FFT
    BEGIN_HIDE_OUTPUT();
    fft->compute(input_data, output_data, FFT3d::FORWARD);
    END_HIDE_OUTPUT();

    // Backward FFT (into input_data to complete round-trip)
    BEGIN_HIDE_OUTPUT();
    fft->compute(output_data, input_data, FFT3d::BACKWARD);
    END_HIDE_OUTPUT();

    // Apply normalization: LAMMPS backward FFT does not include 1/N³ scaling
    // For round-trip correctness: IFFT(FFT(x)) = N³ × x, so we divide by N³
    FFT_SCALAR norm = 1.0 / static_cast<FFT_SCALAR>(nsize);
    for (int i = 0; i < 2 * nsize; i++) {
        input_data[i] *= norm;
    }

    // Validate round-trip
    FFTValidation::RoundTripValidator validator(original_data.data(), input_data, nfast, nmid, nslow,
                                                 ROUNDTRIP_TOLERANCE, verbose);
    bool passed = validator.validate();

    // Report results
    if (verbose || !passed) {
        std::cout << "Round-trip test (64x64x64):" << std::endl;
        std::cout << "  Max error: " << validator.get_error_stats().max() << " (at index "
                  << validator.get_error_stats().idx() << ")" << std::endl;
        std::cout << "  Avg error: " << validator.get_error_stats().avg() << std::endl;
        std::cout << "  Tolerance: " << ROUNDTRIP_TOLERANCE << std::endl;
        std::cout << "  Status: " << (passed ? "PASSED" : "FAILED") << std::endl;
    }

    EXPECT_TRUE(passed) << "Round-trip validation failed";
    EXPECT_LT(validator.get_error_stats().max(), ROUNDTRIP_TOLERANCE);
}

// ============================================================================
// Test 4: Non-power-of-2 Grid (48x48x48)
// ============================================================================

TEST_F(FFT3DTest, RoundTrip_Serial_48x48x48)
{
    // Create 48x48x48 grid (non-power-of-2)
    create_serial_fft(48, 48, 48);

    int nsize = nfast * nmid * nslow;

    // Generate random complex data using RandomComplexGenerator
    FFTTestData::RandomComplexGenerator generator(98765);
    generator.generate(input_data, nfast, nmid, nslow);

    // Copy input data
    std::vector<FFT_SCALAR> original_data(2 * nsize);
    std::copy(input_data, input_data + 2 * nsize, original_data.data());

    // Forward FFT
    BEGIN_HIDE_OUTPUT();
    fft->compute(input_data, output_data, FFT3d::FORWARD);
    END_HIDE_OUTPUT();

    // Backward FFT
    BEGIN_HIDE_OUTPUT();
    fft->compute(output_data, input_data, FFT3d::BACKWARD);
    END_HIDE_OUTPUT();

    // Apply normalization: LAMMPS backward FFT does not include 1/N³ scaling
    // For round-trip correctness: IFFT(FFT(x)) = N³ × x, so we divide by N³
    FFT_SCALAR norm = 1.0 / static_cast<FFT_SCALAR>(nsize);
    for (int i = 0; i < 2 * nsize; i++) {
        input_data[i] *= norm;
    }

    // Validate round-trip
    FFTValidation::RoundTripValidator validator(original_data.data(), input_data, nfast, nmid, nslow,
                                                 ROUNDTRIP_TOLERANCE, verbose);
    bool passed = validator.validate();

    if (verbose || !passed) {
        std::cout << "Round-trip test (48x48x48, non-power-of-2):" << std::endl;
        std::cout << "  Max error: " << validator.get_error_stats().max() << std::endl;
        std::cout << "  Avg error: " << validator.get_error_stats().avg() << std::endl;
        std::cout << "  Status: " << (passed ? "PASSED" : "FAILED") << std::endl;
    }

    EXPECT_TRUE(passed) << "Round-trip validation failed for non-power-of-2 grid";
    EXPECT_LT(validator.get_error_stats().max(), ROUNDTRIP_TOLERANCE);
}
// ============================================================================
// Test 5: Known Answer - Delta Function
// ============================================================================
// FFT property: FFT(δ(x=0)) = constant in all frequency bins
// A spike at the origin transforms to a constant spectrum

TEST_F(FFT3DTest, KnownAnswer_DeltaFunction)
{
    // Create 32x32x32 grid
    create_serial_fft(32, 32, 32);

    int nsize = nfast * nmid * nslow;

    // Generate delta function (spike at origin)
    FFTTestData::DeltaFunctionGenerator generator(1.0);
    generator.generate(input_data, nfast, nmid, nslow);

    // Forward FFT
    BEGIN_HIDE_OUTPUT();
    fft->compute(input_data, output_data, FFT3d::FORWARD);
    END_HIDE_OUTPUT();

    // FFT property: FFT(δ(x=0)) should produce constant magnitude in all frequency bins
    // Expected value: all bins should have the same magnitude (equal to 1.0 for unnormalized FFT)
    // Note: LAMMPS FFT is unnormalized, so FFT(δ) = 1.0 everywhere

    // Create expected output: all bins = (1.0, 0.0)
    std::vector<FFT_SCALAR> expected_data(2 * nsize);
    for (int i = 0; i < nsize; i++) {
        set_complex_linear(expected_data.data(), i, std::complex<FFT_SCALAR>(1.0, 0.0));
    }

    // Validate against expected result
    FFTValidation::KnownAnswerValidator validator(output_data, expected_data.data(), nfast, nmid,
                                                   nslow, 1e-10, verbose);
    bool passed = validator.validate();

    if (verbose || !passed) {
        std::cout << "Known answer test (Delta function):" << std::endl;
        std::cout << "  Expected: FFT(δ(x=0)) = constant = 1.0 everywhere" << std::endl;
        std::cout << "  Max error: " << validator.get_error_stats().max() << " (at index "
                  << validator.get_error_stats().idx() << ")" << std::endl;
        std::cout << "  Avg error: " << validator.get_error_stats().avg() << std::endl;
        std::cout << "  Tolerance: 1e-10" << std::endl;
        std::cout << "  Status: " << (passed ? "PASSED" : "FAILED") << std::endl;

        // Show sample values
        std::cout << "  Sample FFT values:" << std::endl;
        for (int i = 0; i < std::min(5, nsize); i++) {
            auto val = get_complex_linear(output_data, i);
            std::cout << "    Bin " << i << ": " << val << std::endl;
        }
    }

    EXPECT_TRUE(passed) << "Delta function known answer validation failed";
    EXPECT_LT(validator.get_error_stats().max(), 1e-10);
}

// ============================================================================
// Test 6: Known Answer - Constant Field
// ============================================================================
// FFT property: FFT(constant) = N³ at k=0, zero elsewhere
// A constant field has all energy in the DC component

TEST_F(FFT3DTest, KnownAnswer_Constant)
{
    // Create 32x32x32 grid
    create_serial_fft(32, 32, 32);

    int nsize = nfast * nmid * nslow;

    // Generate constant field (all values = 1.0)
    FFTTestData::ConstantGenerator generator(1.0);
    generator.generate(input_data, nfast, nmid, nslow);

    // Forward FFT
    BEGIN_HIDE_OUTPUT();
    fft->compute(input_data, output_data, FFT3d::FORWARD);
    END_HIDE_OUTPUT();

    // FFT property: FFT(constant=1.0) should produce:
    // - Bin (0,0,0): N³ (sum of all input values)
    // - All other bins: 0
    // Note: LAMMPS FFT is unnormalized

    // Create expected output
    std::vector<FFT_SCALAR> expected_data(2 * nsize, 0.0);
    FFT_SCALAR dc_value = static_cast<FFT_SCALAR>(nsize);  // N³
    set_complex(expected_data.data(), 0, 0, 0, nfast, nmid,
                std::complex<FFT_SCALAR>(dc_value, 0.0));

    // Validate against expected result
    FFTValidation::KnownAnswerValidator validator(output_data, expected_data.data(), nfast, nmid,
                                                   nslow, 1e-10, verbose);
    bool passed = validator.validate();

    if (verbose || !passed) {
        std::cout << "Known answer test (Constant field):" << std::endl;
        std::cout << "  Expected: FFT(1.0) = " << dc_value << " at k=0, zero elsewhere"
                  << std::endl;
        std::cout << "  Max error: " << validator.get_error_stats().max() << " (at index "
                  << validator.get_error_stats().idx() << ")" << std::endl;
        std::cout << "  Avg error: " << validator.get_error_stats().avg() << std::endl;
        std::cout << "  Tolerance: 1e-10" << std::endl;
        std::cout << "  Status: " << (passed ? "PASSED" : "FAILED") << std::endl;

        // Show DC component
        auto dc_component = get_complex(output_data, 0, 0, 0, nfast, nmid);
        std::cout << "  DC component: " << dc_component << " (expected: " << dc_value << ")"
                  << std::endl;

        // Show sample of other bins (should be ~0)
        std::cout << "  Sample non-DC values:" << std::endl;
        auto val1 = get_complex(output_data, 1, 0, 0, nfast, nmid);
        auto val2 = get_complex(output_data, 0, 1, 0, nfast, nmid);
        auto val3 = get_complex(output_data, 0, 0, 1, nfast, nmid);
        std::cout << "    Bin (1,0,0): " << val1 << std::endl;
        std::cout << "    Bin (0,1,0): " << val2 << std::endl;
        std::cout << "    Bin (0,0,1): " << val3 << std::endl;
    }

    EXPECT_TRUE(passed) << "Constant field known answer validation failed";
    EXPECT_LT(validator.get_error_stats().max(), 1e-10);
}

// ============================================================================
// Test 7: Known Answer - Sine Wave
// ============================================================================
// FFT property: FFT(sin(k₀·x)) = spikes at ±k₀
// A single-frequency sine wave produces two symmetric peaks in frequency space

TEST_F(FFT3DTest, KnownAnswer_SineWave)
{
    // Create 32x32x32 grid
    create_serial_fft(32, 32, 32);

    int nsize = nfast * nmid * nslow;

    // Generate sine wave: sin(2πkx·x) with kx=2
    // This corresponds to mode (2,0,0) in frequency space
    FFTTestData::SineWaveGenerator generator(2, 0, 0, 1.0);
    generator.generate(input_data, nfast, nmid, nslow);

    // Forward FFT
    BEGIN_HIDE_OUTPUT();
    fft->compute(input_data, output_data, FFT3d::FORWARD);
    END_HIDE_OUTPUT();

    // FFT property: sin(k·x) = [exp(ik·x) - exp(-ik·x)] / (2i)
    // FFT(sin(k₀·x)) produces spikes at k=+k₀ and k=-k₀
    // For sin, the spikes are purely imaginary with opposite signs
    // Spike at k=(2,0,0): -i·N³/2
    // Spike at k=(-2,0,0) = (nfast-2,0,0): +i·N³/2
    // All other bins: 0

    // Expected amplitude for unnormalized FFT of sin wave
    FFT_SCALAR spike_amplitude = static_cast<FFT_SCALAR>(nsize) / 2.0;

    // Create expected output: zeros everywhere except at ±k₀
    std::vector<FFT_SCALAR> expected_data(2 * nsize, 0.0);

    // Positive frequency component (2,0,0): -i·N³/2 = (0, -N³/2)
    set_complex(expected_data.data(), 2, 0, 0, nfast, nmid,
                std::complex<FFT_SCALAR>(0.0, -spike_amplitude));

    // Negative frequency component (-2,0,0) = (nfast-2,0,0): +i·N³/2 = (0, +N³/2)
    set_complex(expected_data.data(), nfast - 2, 0, 0, nfast, nmid,
                std::complex<FFT_SCALAR>(0.0, spike_amplitude));

    // Validate against expected result
    FFTValidation::KnownAnswerValidator validator(output_data, expected_data.data(), nfast, nmid,
                                                   nslow, 1e-10, verbose);
    bool passed = validator.validate();

    if (verbose || !passed) {
        std::cout << "Known answer test (Sine wave):" << std::endl;
        std::cout << "  Input: sin(2π·2·x/N) - mode (2,0,0)" << std::endl;
        std::cout << "  Expected: Spikes at k=(2,0,0) and k=(-2,0,0)=(30,0,0)" << std::endl;
        std::cout << "  Spike amplitude: " << spike_amplitude << std::endl;
        std::cout << "  Max error: " << validator.get_error_stats().max() << " (at index "
                  << validator.get_error_stats().idx() << ")" << std::endl;
        std::cout << "  Avg error: " << validator.get_error_stats().avg() << std::endl;
        std::cout << "  Tolerance: 1e-10" << std::endl;
        std::cout << "  Status: " << (passed ? "PASSED" : "FAILED") << std::endl;

        // Show the spike values
        auto spike_pos = get_complex(output_data, 2, 0, 0, nfast, nmid);
        auto spike_neg = get_complex(output_data, nfast - 2, 0, 0, nfast, nmid);
        std::cout << "  Spike at (2,0,0): " << spike_pos << " (expected: (0, "
                  << -spike_amplitude << "))" << std::endl;
        std::cout << "  Spike at (" << (nfast - 2) << ",0,0): " << spike_neg
                  << " (expected: (0, " << spike_amplitude << "))" << std::endl;

        // Show sample of bins that should be zero
        std::cout << "  Sample zero bins:" << std::endl;
        auto val1 = get_complex(output_data, 0, 0, 0, nfast, nmid);
        auto val2 = get_complex(output_data, 1, 0, 0, nfast, nmid);
        auto val3 = get_complex(output_data, 3, 0, 0, nfast, nmid);
        std::cout << "    Bin (0,0,0): " << val1 << std::endl;
        std::cout << "    Bin (1,0,0): " << val2 << std::endl;
        std::cout << "    Bin (3,0,0): " << val3 << std::endl;
    }

    EXPECT_TRUE(passed) << "Sine wave known answer validation failed";
    EXPECT_LT(validator.get_error_stats().max(), 1e-10);
}


// ============================================================================
// Test 8: Parseval's Theorem - Energy Conservation (32x32x32)
// ============================================================================

TEST_F(FFT3DTest, ParsevalsTheorem_EnergyConservation)
{
    // Create 32x32x32 grid
    create_serial_fft(32, 32, 32);

    int nsize = nfast * nmid * nslow;

    // Generate random complex data using RandomComplexGenerator
    // This ensures we have non-trivial data with both real and imaginary components
    FFTTestData::RandomComplexGenerator generator(42424);
    generator.generate(input_data, nfast, nmid, nslow);

    // Copy input data (spatial domain)
    std::vector<FFT_SCALAR> spatial_data(2 * nsize);
    std::copy(input_data, input_data + 2 * nsize, spatial_data.data());

    // Compute energy in spatial domain: E_spatial = Σ|x[i]|²
    double spatial_energy = 0.0;
    for (int i = 0; i < nsize; i++) {
        std::complex<FFT_SCALAR> val = FFTTestHelpers::get_complex_linear(spatial_data.data(), i);
        spatial_energy += std::norm(val);  // |z|² = re² + im²
    }

    // Forward FFT
    BEGIN_HIDE_OUTPUT();
    fft->compute(input_data, output_data, FFT3d::FORWARD);
    END_HIDE_OUTPUT();

    // Compute energy in frequency domain: E_freq = Σ|X[i]|²
    double frequency_energy = 0.0;
    for (int i = 0; i < nsize; i++) {
        std::complex<FFT_SCALAR> val = FFTTestHelpers::get_complex_linear(output_data, i);
        frequency_energy += std::norm(val);  // |Z|² = Re² + Im²
    }

    // Apply Parseval's theorem normalization: E_spatial = (1/N³) × E_freq
    double n_cubed = static_cast<double>(nsize);
    double frequency_energy_normalized = frequency_energy / n_cubed;

    // Calculate relative error
    double abs_error = std::abs(spatial_energy - frequency_energy_normalized);
    double relative_error = (spatial_energy > 1e-14) ? abs_error / spatial_energy : abs_error;

    // Validate using ParsevalValidator
    FFTValidation::ParsevalValidator validator(spatial_data.data(), output_data, nfast, nmid, nslow,
                                                PARSEVAL_TOLERANCE, verbose);
    bool passed = validator.validate();

    // Report results
    if (verbose || !passed) {
        std::cout << "Parseval's theorem test (32x32x32):" << std::endl;
        std::cout << "  Grid size (N³):                 " << nsize << std::endl;
        std::cout << "  Spatial energy (Σ|x|²):         " << spatial_energy << std::endl;
        std::cout << "  Frequency energy (Σ|X|²):       " << frequency_energy << std::endl;
        std::cout << "  Normalized frequency (Σ|X|²/N³): " << frequency_energy_normalized << std::endl;
        std::cout << "  Relative error:                 " << relative_error << std::endl;
        std::cout << "  Tolerance:                      " << PARSEVAL_TOLERANCE << std::endl;
        std::cout << "  Status: " << (passed ? "PASSED" : "FAILED") << std::endl;
    }

    EXPECT_TRUE(passed) << "Parseval's theorem validation failed";
    EXPECT_LT(relative_error, PARSEVAL_TOLERANCE);
}

// ============================================================================
// Test 9: Round-Trip MPI (2 processes, 32x32x32)
// ============================================================================

TEST_F(FFT3DTest, RoundTrip_MPI_2proc_32x32x32)
{
    // Check MPI environment
    int nprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (nprocs != 2) {
        GTEST_SKIP() << "Test requires exactly 2 MPI processes, got " << nprocs;
    }

    // Grid dimensions
    nfast = 32;
    nmid = 32;
    nslow = 32;

    // Domain decomposition: split along slow (z) dimension
    // Each process owns half of the z-slices
    // Process 0: z = 0..15
    // Process 1: z = 16..31
    int in_ilo = 0, in_ihi = nfast - 1;
    int in_jlo = 0, in_jhi = nmid - 1;
    int in_klo, in_khi;

    if (rank == 0) {
        in_klo = 0;
        in_khi = nslow / 2 - 1;  // 0..15
    } else {
        in_klo = nslow / 2;
        in_khi = nslow - 1;      // 16..31
    }

    // Output decomposition: same as input for simplicity
    int out_ilo = in_ilo, out_ihi = in_ihi;
    int out_jlo = in_jlo, out_jhi = in_jhi;
    int out_klo = in_klo, out_khi = in_khi;

    // Calculate local size: each rank owns a slab
    int local_nfast = nfast;
    int local_nmid = nmid;
    int local_nslow = in_khi - in_klo + 1;
    int local_size = local_nfast * local_nmid * local_nslow;

    if (verbose) {
        std::cout << "Rank " << rank << ": local grid = " << local_nfast << "x" << local_nmid
                  << "x" << local_nslow << " (z-range: " << in_klo << ".." << in_khi << ")"
                  << std::endl;
    }

    // Allocate local data buffers
    input_data = new FFT_SCALAR[2 * local_size];
    output_data = new FFT_SCALAR[2 * local_size];
    std::memset(input_data, 0, 2 * local_size * sizeof(FFT_SCALAR));
    std::memset(output_data, 0, 2 * local_size * sizeof(FFT_SCALAR));

    // FFT parameters
    int scaled = 0;      // No scaling
    int permute = 0;     // No permutation
    int nbuf = 0;        // Buffer size (output)
    int usecollective = 0;  // Use point-to-point communication

    // Create MPI-aware FFT3d object
    BEGIN_HIDE_OUTPUT();
    fft = new FFT3d(lmp, MPI_COMM_WORLD, nfast, nmid, nslow, in_ilo, in_ihi, in_jlo, in_jhi,
                    in_klo, in_khi, out_ilo, out_ihi, out_jlo, out_jhi, out_klo, out_khi,
                    scaled, permute, &nbuf, usecollective);
    END_HIDE_OUTPUT();

    ASSERT_NE(fft, nullptr);

    // Generate random complex data in local portion
    // Use deterministic function based on global coordinates for consistency across ranks
    for (int k = 0; k < local_nslow; k++) {
        for (int j = 0; j < local_nmid; j++) {
            for (int i = 0; i < local_nfast; i++) {
                // Calculate global indices for this point
                int global_k = in_klo + k;

                // Generate value based on global position (for consistency)
                int global_idx = global_k * nmid * nfast + j * nfast + i;
                std::mt19937 rng(12345 + global_idx);
                std::uniform_real_distribution<FFT_SCALAR> dist(-1.0, 1.0);

                FFT_SCALAR re = dist(rng);
                FFT_SCALAR im = dist(rng);

                int local_idx = k * nmid * nfast + j * nfast + i;
                input_data[2 * local_idx] = re;
                input_data[2 * local_idx + 1] = im;
            }
        }
    }

    // Copy input data for later comparison
    std::vector<FFT_SCALAR> original_data(2 * local_size);
    std::copy(input_data, input_data + 2 * local_size, original_data.data());

    // Forward FFT
    BEGIN_HIDE_OUTPUT();
    fft->compute(input_data, output_data, FFT3d::FORWARD);
    END_HIDE_OUTPUT();

    // Backward FFT (into input_data to complete round-trip)
    BEGIN_HIDE_OUTPUT();
    fft->compute(output_data, input_data, FFT3d::BACKWARD);
    END_HIDE_OUTPUT();

    // Apply normalization: LAMMPS backward FFT does not include 1/N³ scaling
    int total_size = nfast * nmid * nslow;
    FFT_SCALAR norm = 1.0 / static_cast<FFT_SCALAR>(total_size);
    for (int i = 0; i < 2 * local_size; i++) {
        input_data[i] *= norm;
    }

    // Validate round-trip on local data
    FFTValidation::RoundTripValidator validator(original_data.data(), input_data, local_nfast,
                                                 local_nmid, local_nslow, ROUNDTRIP_TOLERANCE,
                                                 verbose);
    bool passed = validator.validate();

    // Gather validation results from all ranks
    int passed_int = passed ? 1 : 0;
    int all_passed;
    MPI_Allreduce(&passed_int, &all_passed, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    // Gather max error from all ranks
    double local_max_error = validator.get_error_stats().max();
    double global_max_error;
    MPI_Allreduce(&local_max_error, &global_max_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // Report results (rank 0 only)
    if (rank == 0 && (verbose || !all_passed)) {
        std::cout << "Round-trip test (MPI 2 procs, 32x32x32):" << std::endl;
        std::cout << "  Global grid: " << nfast << "x" << nmid << "x" << nslow << std::endl;
        std::cout << "  MPI processes: " << nprocs << std::endl;
        std::cout << "  Global max error: " << global_max_error << std::endl;
        std::cout << "  Tolerance: " << ROUNDTRIP_TOLERANCE << std::endl;
        std::cout << "  Status: " << (all_passed ? "PASSED" : "FAILED") << std::endl;
    }

    EXPECT_TRUE(all_passed) << "Round-trip validation failed on rank " << rank;
    EXPECT_LT(global_max_error, ROUNDTRIP_TOLERANCE);
}

// ============================================================================
// Test 10: Round-Trip MPI (4 processes, 64x64x64)
// ============================================================================

TEST_F(FFT3DTest, RoundTrip_MPI_4proc_64x64x64)
{
    // Check MPI environment
    int nprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (nprocs != 4) {
        GTEST_SKIP() << "Test requires exactly 4 MPI processes, got " << nprocs;
    }

    // Grid dimensions (larger for better load balancing with 4 procs)
    nfast = 64;
    nmid = 64;
    nslow = 64;

    // Domain decomposition: split along slow (z) dimension
    // Each process owns 1/4 of the z-slices
    // Process 0: z = 0..15
    // Process 1: z = 16..31
    // Process 2: z = 32..47
    // Process 3: z = 48..63
    int in_ilo = 0, in_ihi = nfast - 1;
    int in_jlo = 0, in_jhi = nmid - 1;
    int in_klo, in_khi;

    int slices_per_proc = nslow / nprocs;
    in_klo = rank * slices_per_proc;
    in_khi = (rank + 1) * slices_per_proc - 1;

    // Output decomposition: same as input
    int out_ilo = in_ilo, out_ihi = in_ihi;
    int out_jlo = in_jlo, out_jhi = in_jhi;
    int out_klo = in_klo, out_khi = in_khi;

    // Calculate local size
    int local_nfast = nfast;
    int local_nmid = nmid;
    int local_nslow = in_khi - in_klo + 1;
    int local_size = local_nfast * local_nmid * local_nslow;

    if (verbose) {
        std::cout << "Rank " << rank << ": local grid = " << local_nfast << "x" << local_nmid
                  << "x" << local_nslow << " (z-range: " << in_klo << ".." << in_khi << ")"
                  << std::endl;
    }

    // Allocate local data buffers
    input_data = new FFT_SCALAR[2 * local_size];
    output_data = new FFT_SCALAR[2 * local_size];
    std::memset(input_data, 0, 2 * local_size * sizeof(FFT_SCALAR));
    std::memset(output_data, 0, 2 * local_size * sizeof(FFT_SCALAR));

    // FFT parameters
    int scaled = 0;      // No scaling
    int permute = 0;     // No permutation
    int nbuf = 0;        // Buffer size (output)
    int usecollective = 0;  // Use point-to-point communication

    // Create MPI-aware FFT3d object
    BEGIN_HIDE_OUTPUT();
    fft = new FFT3d(lmp, MPI_COMM_WORLD, nfast, nmid, nslow, in_ilo, in_ihi, in_jlo, in_jhi,
                    in_klo, in_khi, out_ilo, out_ihi, out_jlo, out_jhi, out_klo, out_khi,
                    scaled, permute, &nbuf, usecollective);
    END_HIDE_OUTPUT();

    ASSERT_NE(fft, nullptr);

    // Generate random complex data in local portion
    // Use deterministic function based on global coordinates for consistency across ranks
    for (int k = 0; k < local_nslow; k++) {
        for (int j = 0; j < local_nmid; j++) {
            for (int i = 0; i < local_nfast; i++) {
                // Calculate global indices
                int global_k = in_klo + k;

                // Generate value based on global position
                int global_idx = global_k * nmid * nfast + j * nfast + i;
                std::mt19937 rng(54321 + global_idx);
                std::uniform_real_distribution<FFT_SCALAR> dist(-1.0, 1.0);

                FFT_SCALAR re = dist(rng);
                FFT_SCALAR im = dist(rng);

                int local_idx = k * nmid * nfast + j * nfast + i;
                input_data[2 * local_idx] = re;
                input_data[2 * local_idx + 1] = im;
            }
        }
    }

    // Copy input data for later comparison
    std::vector<FFT_SCALAR> original_data(2 * local_size);
    std::copy(input_data, input_data + 2 * local_size, original_data.data());

    // Forward FFT
    BEGIN_HIDE_OUTPUT();
    fft->compute(input_data, output_data, FFT3d::FORWARD);
    END_HIDE_OUTPUT();

    // Backward FFT (into input_data to complete round-trip)
    BEGIN_HIDE_OUTPUT();
    fft->compute(output_data, input_data, FFT3d::BACKWARD);
    END_HIDE_OUTPUT();

    // Apply normalization
    int total_size = nfast * nmid * nslow;
    FFT_SCALAR norm = 1.0 / static_cast<FFT_SCALAR>(total_size);
    for (int i = 0; i < 2 * local_size; i++) {
        input_data[i] *= norm;
    }

    // Validate round-trip on local data
    FFTValidation::RoundTripValidator validator(original_data.data(), input_data, local_nfast,
                                                 local_nmid, local_nslow, ROUNDTRIP_TOLERANCE,
                                                 verbose);
    bool passed = validator.validate();

    // Gather validation results from all ranks
    int passed_int = passed ? 1 : 0;
    int all_passed;
    MPI_Allreduce(&passed_int, &all_passed, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    // Gather max error from all ranks
    double local_max_error = validator.get_error_stats().max();
    double global_max_error;
    MPI_Allreduce(&local_max_error, &global_max_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // Report results (rank 0 only)
    if (rank == 0 && (verbose || !all_passed)) {
        std::cout << "Round-trip test (MPI 4 procs, 64x64x64):" << std::endl;
        std::cout << "  Global grid: " << nfast << "x" << nmid << "x" << nslow << std::endl;
        std::cout << "  MPI processes: " << nprocs << std::endl;
        std::cout << "  Slices per process: " << slices_per_proc << std::endl;
        std::cout << "  Global max error: " << global_max_error << std::endl;
        std::cout << "  Tolerance: " << ROUNDTRIP_TOLERANCE << std::endl;
        std::cout << "  Status: " << (all_passed ? "PASSED" : "FAILED") << std::endl;
    }

    EXPECT_TRUE(all_passed) << "Round-trip validation failed on rank " << rank;
    EXPECT_LT(global_max_error, ROUNDTRIP_TOLERANCE);
}

// ============================================================================
// Library-Specific Tests
// ============================================================================
// These tests exercise features specific to particular FFT libraries
// (FFTW3, MKL, KISS, etc.). Tests are conditionally compiled and skip
// if the required library is not available.

// ============================================================================
// Test 11: FFTW3 Threading Support
// ============================================================================
// Test FFTW3 with threading enabled (if configured)
// Validates that threaded FFTW3 produces correct results

TEST_F(FFT3DTest, FFTW3_Threading)
{
#ifndef FFT_FFTW3
    GTEST_SKIP() << "Test requires FFTW3, built with: " << LMP_FFT_LIB;
#endif

#ifndef FFT_FFTW_THREADS
    GTEST_SKIP() << "Test requires FFTW3 with threading enabled (FFT_FFTW_THREADS)";
#endif

    // Create 64x64x64 grid (large enough to benefit from threading)
    create_serial_fft(64, 64, 64);

    int nsize = nfast * nmid * nslow;

    // Generate random complex data
    FFTTestData::RandomComplexGenerator generator(77777);
    generator.generate(input_data, nfast, nmid, nslow);

    // Copy input data for later comparison
    std::vector<FFT_SCALAR> original_data(2 * nsize);
    std::copy(input_data, input_data + 2 * nsize, original_data.data());

    // Forward FFT (threaded FFTW3)
    BEGIN_HIDE_OUTPUT();
    fft->compute(input_data, output_data, FFT3d::FORWARD);
    END_HIDE_OUTPUT();

    // Backward FFT
    BEGIN_HIDE_OUTPUT();
    fft->compute(output_data, input_data, FFT3d::BACKWARD);
    END_HIDE_OUTPUT();

    // Apply normalization
    FFT_SCALAR norm = 1.0 / static_cast<FFT_SCALAR>(nsize);
    for (int i = 0; i < 2 * nsize; i++) {
        input_data[i] *= norm;
    }

    // Validate round-trip (threaded FFTW3 should give identical results)
    FFTValidation::RoundTripValidator validator(original_data.data(), input_data, nfast, nmid,
                                                 nslow, ROUNDTRIP_TOLERANCE, verbose);
    bool passed = validator.validate();

    EXPECT_TRUE(passed) << "FFTW3 threaded round-trip validation failed"
                        << "\n  Max error: " << validator.get_error_stats().max()
                        << "\n  Tolerance: " << ROUNDTRIP_TOLERANCE;
    EXPECT_LT(validator.get_error_stats().max(), ROUNDTRIP_TOLERANCE);
}

// ============================================================================
// Test 12: MKL Optimized FFT
// ============================================================================
// Test MKL library-specific optimizations
// Validates that MKL FFT produces correct results

TEST_F(FFT3DTest, MKL_Optimized)
{
#ifndef FFT_MKL
    GTEST_SKIP() << "Test requires MKL FFT, built with: " << LMP_FFT_LIB;
#endif

    // Create 64x64x64 grid (tests MKL optimizations for larger grids)
    create_serial_fft(64, 64, 64);

    int nsize = nfast * nmid * nslow;

    // Generate random complex data
    FFTTestData::RandomComplexGenerator generator(88888);
    generator.generate(input_data, nfast, nmid, nslow);

    // Copy input data
    std::vector<FFT_SCALAR> original_data(2 * nsize);
    std::copy(input_data, input_data + 2 * nsize, original_data.data());

    // Forward FFT (MKL optimized)
    BEGIN_HIDE_OUTPUT();
    fft->compute(input_data, output_data, FFT3d::FORWARD);
    END_HIDE_OUTPUT();

    // Test known answer: delta function
    // This validates MKL's correctness, not just round-trip
    std::vector<FFT_SCALAR> delta_input(2 * nsize, 0.0);
    std::vector<FFT_SCALAR> delta_output(2 * nsize, 0.0);

    // Create delta function
    FFTTestData::DeltaFunctionGenerator delta_gen(1.0);
    delta_gen.generate(delta_input.data(), nfast, nmid, nslow);

    // Compute FFT of delta (should be constant)
    BEGIN_HIDE_OUTPUT();
    fft->compute(delta_input.data(), delta_output.data(), FFT3d::FORWARD);
    END_HIDE_OUTPUT();

    // Validate: FFT(δ) = constant
    std::vector<FFT_SCALAR> expected_data(2 * nsize);
    for (int i = 0; i < nsize; i++) {
        set_complex_linear(expected_data.data(), i, std::complex<FFT_SCALAR>(1.0, 0.0));
    }

    FFTValidation::KnownAnswerValidator validator(delta_output.data(), expected_data.data(),
                                                   nfast, nmid, nslow, 1e-10, verbose);
    bool passed = validator.validate();

    EXPECT_TRUE(passed) << "MKL delta function validation failed"
                        << "\n  Max error: " << validator.get_error_stats().max()
                        << "\n  Tolerance: " << 1e-10;
    EXPECT_LT(validator.get_error_stats().max(), 1e-10);
}

// ============================================================================
// Test 13: KISS FFT Non-Power-of-2 Sizes
// ============================================================================
// Test KISS FFT with various non-power-of-2 grid sizes
// KISS FFT handles arbitrary sizes, unlike some optimized libraries

TEST_F(FFT3DTest, KISS_NonPowerOf2)
{
#ifndef FFT_KISS
    GTEST_SKIP() << "Test requires KISS FFT, built with: " << LMP_FFT_LIB;
#endif

    // Test multiple non-power-of-2 sizes
    std::vector<int> test_sizes = {30, 50, 100, 128};

    for (int size : test_sizes) {
        // Clean up previous FFT
        if (fft) {
            delete fft;
            fft = nullptr;
        }
        if (input_data) {
            delete[] input_data;
            input_data = nullptr;
        }
        if (output_data) {
            delete[] output_data;
            output_data = nullptr;
        }

        // Create grid with current size
        create_serial_fft(size, size, size);

        int nsize = nfast * nmid * nslow;

        if (verbose) {
            std::cout << "  Testing size: " << size << "x" << size << "x" << size << " (N³="
                      << nsize << ")" << std::endl;
        }

        // Generate random complex data
        FFTTestData::RandomComplexGenerator generator(99999 + size);
        generator.generate(input_data, nfast, nmid, nslow);

        // Copy input data
        std::vector<FFT_SCALAR> original_data(2 * nsize);
        std::copy(input_data, input_data + 2 * nsize, original_data.data());

        // Forward FFT
        BEGIN_HIDE_OUTPUT();
        fft->compute(input_data, output_data, FFT3d::FORWARD);
        END_HIDE_OUTPUT();

        // Backward FFT
        BEGIN_HIDE_OUTPUT();
        fft->compute(output_data, input_data, FFT3d::BACKWARD);
        END_HIDE_OUTPUT();

        // Apply normalization
        FFT_SCALAR norm = 1.0 / static_cast<FFT_SCALAR>(nsize);
        for (int i = 0; i < 2 * nsize; i++) {
            input_data[i] *= norm;
        }

        // Validate round-trip
        FFTValidation::RoundTripValidator validator(original_data.data(), input_data, nfast, nmid,
                                                     nslow, ROUNDTRIP_TOLERANCE, verbose);
        bool passed = validator.validate();

        if (verbose || !passed) {
            std::cout << "    Size " << size << ": Max error = " << validator.get_error_stats().max()
                      << " (tolerance = " << ROUNDTRIP_TOLERANCE << ")" << std::endl;
            std::cout << "    Status: " << (passed ? "PASSED" : "FAILED") << std::endl;
        }

        EXPECT_TRUE(passed) << "KISS FFT round-trip failed for size " << size;
        EXPECT_LT(validator.get_error_stats().max(), ROUNDTRIP_TOLERANCE);
    }
}

// ============================================================================
// Test 14: HeFFTe Distributed FFT (Optional)
// ============================================================================
// Test HeFFTe distributed FFT capabilities
// HeFFTe is designed for distributed-memory parallel FFT

TEST_F(FFT3DTest, HeFFTe_Distributed)
{
#ifndef FFT_HEFFTE
    GTEST_SKIP() << "Test requires HeFFTe, built with: " << LMP_FFT_LIB;
#endif

    // Note: HeFFTe is designed for distributed-memory parallelism
    // In serial mode, we just validate basic correctness
    // Full distributed testing would require MPI with multiple processes

    // Create 32x32x32 grid
    create_serial_fft(32, 32, 32);

    int nsize = nfast * nmid * nslow;

    // Generate random complex data
    FFTTestData::RandomComplexGenerator generator(55555);
    generator.generate(input_data, nfast, nmid, nslow);

    // Copy input data
    std::vector<FFT_SCALAR> original_data(2 * nsize);
    std::copy(input_data, input_data + 2 * nsize, original_data.data());

    // Forward FFT (HeFFTe)
    BEGIN_HIDE_OUTPUT();
    fft->compute(input_data, output_data, FFT3d::FORWARD);
    END_HIDE_OUTPUT();

    // Backward FFT
    BEGIN_HIDE_OUTPUT();
    fft->compute(output_data, input_data, FFT3d::BACKWARD);
    END_HIDE_OUTPUT();

    // Apply normalization
    FFT_SCALAR norm = 1.0 / static_cast<FFT_SCALAR>(nsize);
    for (int i = 0; i < 2 * nsize; i++) {
        input_data[i] *= norm;
    }

    // Validate round-trip
    FFTValidation::RoundTripValidator validator(original_data.data(), input_data, nfast, nmid,
                                                 nslow, ROUNDTRIP_TOLERANCE, verbose);
    bool passed = validator.validate();

    EXPECT_TRUE(passed) << "HeFFTe round-trip validation failed"
                        << "\n  Max error: " << validator.get_error_stats().max()
                        << "\n  Tolerance: " << ROUNDTRIP_TOLERANCE
                        << "\n  Note: Running in serial mode";
    EXPECT_LT(validator.get_error_stats().max(), ROUNDTRIP_TOLERANCE);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);

    // Handle command line options
    if (argc > 1 && strcmp(argv[1], "-v") == 0) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
