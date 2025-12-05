/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// FFT Test Helpers - Test data generators, validators, and utilities for FFT testing
// Implements comprehensive validation strategies for LAMMPS FFT wrappers
//
// Documentation: See doc/src/Developer_unittest.rst
//                (Section: FFT Testing Infrastructure)

#ifndef LMP_FFT_TEST_HELPERS_H
#define LMP_FFT_TEST_HELPERS_H

#include "../force-styles/error_stats.h"
#include "lmpfftsettings.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <random>
#include <string>
#include <vector>

// =============================================================================
// NAMESPACE: FFTTestHelpers
// =============================================================================

namespace FFTTestHelpers {

// =============================================================================
// SECTION 1: Utility Functions and Helpers
// =============================================================================

// Constants
constexpr double PI     = 3.14159265358979323846;
constexpr double TWO_PI = 2.0 * PI;

// Precision-aware tolerances (based on FFT_SCALAR type)
#ifdef FFT_SINGLE
constexpr FFT_SCALAR ROUNDTRIP_TOLERANCE    = 1e-5;
constexpr FFT_SCALAR KNOWN_ANSWER_TOLERANCE = 1e-6;
constexpr FFT_SCALAR PARSEVAL_TOLERANCE     = 1e-5;
constexpr FFT_SCALAR HERMITIAN_TOLERANCE    = 1e-6;
constexpr FFT_SCALAR LINEARITY_TOLERANCE    = 1e-5;
#else
constexpr FFT_SCALAR ROUNDTRIP_TOLERANCE    = 1e-11;
constexpr FFT_SCALAR KNOWN_ANSWER_TOLERANCE = 1e-12;
constexpr FFT_SCALAR PARSEVAL_TOLERANCE     = 1e-10;
constexpr FFT_SCALAR HERMITIAN_TOLERANCE    = 1e-12;
constexpr FFT_SCALAR LINEARITY_TOLERANCE    = 1e-11;
#endif

// Grid size scaling factor for tolerances (larger grids accumulate more error)
inline FFT_SCALAR grid_scale_factor(int nfast, int nmid, int nslow)
{
    int n_total = nfast * nmid * nslow;
    // Scale tolerance by sqrt(N) for accumulated floating-point errors
    return std::sqrt(static_cast<FFT_SCALAR>(n_total) / static_cast<FFT_SCALAR>(1000.0));
}

// Calculate scaled tolerance based on grid size
inline FFT_SCALAR scaled_tolerance(FFT_SCALAR base_tol, int nfast, int nmid, int nslow)
{
    return base_tol * grid_scale_factor(nfast, nmid, nslow);
}

// =============================================================================
// SUBSECTION 1.1: Index and Complex Number Utilities
// =============================================================================

// Convert 3D grid indices to 1D array index
// FFT data is stored as: [re0, im0, re1, im1, ..., re(N-1), im(N-1)]
// Layout: i varies fastest, then j, then k (C-style row-major)
inline int idx3d(int i, int j, int k, int nfast, int nmid)
{
    return 2 * (k * nmid * nfast + j * nfast + i);
}

// Get complex value at grid point (i, j, k)
inline std::complex<FFT_SCALAR> get_complex(const FFT_SCALAR *data, int i, int j, int k, int nfast,
                                            int nmid)
{
    int idx = idx3d(i, j, k, nfast, nmid);
    return std::complex<FFT_SCALAR>(data[idx], data[idx + 1]);
}

// Set complex value at grid point (i, j, k)
inline void set_complex(FFT_SCALAR *data, int i, int j, int k, int nfast, int nmid,
                        const std::complex<FFT_SCALAR> &value)
{
    int idx       = idx3d(i, j, k, nfast, nmid);
    data[idx]     = value.real();
    data[idx + 1] = value.imag();
}

// Get complex value at linear index (for direct array access)
inline std::complex<FFT_SCALAR> get_complex_linear(const FFT_SCALAR *data, int linear_idx)
{
    return std::complex<FFT_SCALAR>(data[2 * linear_idx], data[2 * linear_idx + 1]);
}

// Set complex value at linear index
inline void set_complex_linear(FFT_SCALAR *data, int linear_idx,
                               const std::complex<FFT_SCALAR> &value)
{
    data[2 * linear_idx]     = value.real();
    data[2 * linear_idx + 1] = value.imag();
}

// =============================================================================
// SUBSECTION 1.2: FFTBuffer - RAII Memory Management
// =============================================================================

// RAII wrapper for FFT data buffers
class FFTBuffer {
public:
    FFTBuffer(int nfast, int nmid, int nslow) :
        nfast_(nfast), nmid_(nmid), nslow_(nslow), size_(2 * nfast * nmid * nslow)
    {
        data_ = new FFT_SCALAR[size_];
        std::fill(data_, data_ + size_, 0.0);
    }

    ~FFTBuffer() { delete[] data_; }

    // Non-copyable (avoid accidental copies of large arrays)
    FFTBuffer(const FFTBuffer &)            = delete;
    FFTBuffer &operator=(const FFTBuffer &) = delete;

    // Movable
    FFTBuffer(FFTBuffer &&other) noexcept :
        data_(other.data_), nfast_(other.nfast_), nmid_(other.nmid_), nslow_(other.nslow_),
        size_(other.size_)
    {
        other.data_ = nullptr;
        other.size_ = 0;
    }

    FFTBuffer &operator=(FFTBuffer &&other) noexcept
    {
        if (this != &other) {
            delete[] data_;
            data_       = other.data_;
            nfast_      = other.nfast_;
            nmid_       = other.nmid_;
            nslow_      = other.nslow_;
            size_       = other.size_;
            other.data_ = nullptr;
            other.size_ = 0;
        }
        return *this;
    }

    // Accessors
    FFT_SCALAR *data() { return data_; }
    const FFT_SCALAR *data() const { return data_; }
    int size() const { return size_; }
    int nfast() const { return nfast_; }
    int nmid() const { return nmid_; }
    int nslow() const { return nslow_; }

    // Zero all data
    void zero() { std::fill(data_, data_ + size_, 0.0); }

    // Copy from another buffer
    void copy_from(const FFT_SCALAR *src) { std::copy(src, src + size_, data_); }

    // Copy to another buffer
    void copy_to(FFT_SCALAR *dest) const { std::copy(data_, data_ + size_, dest); }

private:
    FFT_SCALAR *data_;
    int nfast_, nmid_, nslow_;
    int size_;
};

// Report FFT configuration to output stream
inline void report_fft_config(std::ostream &out)
{
    out << "=== FFT Configuration ===" << std::endl;
    out << "Library:   " << LMP_FFT_LIB << std::endl;
    out << "Precision: " << LMP_FFT_PREC << std::endl;

#if defined(FFT_FFTW_THREADS)
    out << "Threading: FFTW threads enabled" << std::endl;
#elif defined(FFT_MKL_THREADS)
    out << "Threading: MKL threads enabled" << std::endl;
#else
    out << "Threading: disabled" << std::endl;
#endif

#if defined(FFT_KOKKOS_CUFFT)
    out << "Backend:   CUDA (cuFFT)" << std::endl;
#elif defined(FFT_KOKKOS_HIPFFT)
    out << "Backend:   HIP (hipFFT)" << std::endl;
#elif defined(FFT_KOKKOS_MKL_GPU)
    out << "Backend:   SYCL (MKL GPU)" << std::endl;
#elif defined(FFT_HEFFTE)
    out << "Backend:   heFFTe (distributed FFT)" << std::endl;
#else
    out << "Backend:   CPU" << std::endl;
#endif

    out << "=========================" << std::endl;
}

} // namespace FFTTestHelpers

// =============================================================================
// NAMESPACE: FFTTestData
// =============================================================================

namespace FFTTestData {

using namespace FFTTestHelpers;

// =============================================================================
// SECTION 2: Test Data Generators
// =============================================================================

// Base class for all test data generators
class DataGenerator {
public:
    virtual ~DataGenerator() {}

    // Generate test data in FFT_SCALAR array (complex interleaved format)
    virtual void generate(FFT_SCALAR *data, int nfast, int nmid, int nslow) = 0;

    // Get total energy (for Parseval tests): sum |data[i]|^2
    virtual double get_energy(const FFT_SCALAR *data, int nfast, int nmid, int nslow) const
    {
        double energy = 0.0;
        int n_points  = nfast * nmid * nslow;
        for (int i = 0; i < n_points; i++) {
            std::complex<FFT_SCALAR> val = get_complex_linear(data, i);
            energy += std::norm(val); // |z|^2 = re^2 + im^2
        }
        return energy;
    }

    // Check if input is purely real (for Hermitian symmetry tests)
    virtual bool is_real_input() const = 0;
};

// =============================================================================
// SUBSECTION 2.1: DeltaFunctionGenerator
// =============================================================================

// Delta function: spike at origin (0,0,0)
// FFT properties: δ(x) → constant spectrum (all frequencies equal amplitude)
class DeltaFunctionGenerator : public DataGenerator {
public:
    DeltaFunctionGenerator(FFT_SCALAR amplitude = 1.0) : amplitude_(amplitude) {}

    void generate(FFT_SCALAR *data, int nfast, int nmid, int nslow) override
    {
        // Zero all data
        int size = 2 * nfast * nmid * nslow;
        std::fill(data, data + size, 0.0);

        // Set spike at origin (0,0,0)
        set_complex(data, 0, 0, 0, nfast, nmid, std::complex<FFT_SCALAR>(amplitude_, 0.0));
    }

    bool is_real_input() const override { return true; }

private:
    FFT_SCALAR amplitude_;
};

// =============================================================================
// SUBSECTION 2.2: ConstantGenerator
// =============================================================================

// Constant: all grid points = same value
// FFT properties: constant → spike at k=0 (DC component only)
class ConstantGenerator : public DataGenerator {
public:
    ConstantGenerator(FFT_SCALAR value = 1.0) : value_(value) {}

    void generate(FFT_SCALAR *data, int nfast, int nmid, int nslow) override
    {
        int n_points = nfast * nmid * nslow;
        for (int i = 0; i < n_points; i++) {
            set_complex_linear(data, i, std::complex<FFT_SCALAR>(value_, 0.0));
        }
    }

    bool is_real_input() const override { return true; }

private:
    FFT_SCALAR value_;
};

// =============================================================================
// SUBSECTION 2.3: SineWaveGenerator
// =============================================================================

// Single-frequency sine wave: sin(k·x)
// FFT properties: sin(k₀·x) → spikes at ±k₀ in frequency space
class SineWaveGenerator : public DataGenerator {
public:
    // Constructor: specify wave vector (kx, ky, kz) as integer multiples of fundamental modes
    SineWaveGenerator(int kx = 1, int ky = 0, int kz = 0, FFT_SCALAR amplitude = 1.0) :
        kx_(kx), ky_(ky), kz_(kz), amplitude_(amplitude)
    {
    }

    void generate(FFT_SCALAR *data, int nfast, int nmid, int nslow) override
    {
        for (int k = 0; k < nslow; k++) {
            for (int j = 0; j < nmid; j++) {
                for (int i = 0; i < nfast; i++) {
                    // Calculate phase: k·x = (kx*i/nfast + ky*j/nmid + kz*k/nslow) * 2π
                    double phase = TWO_PI * (kx_ * static_cast<double>(i) / nfast +
                                             ky_ * static_cast<double>(j) / nmid +
                                             kz_ * static_cast<double>(k) / nslow);

                    // sin(phase) in real part, 0 in imaginary part
                    FFT_SCALAR value = amplitude_ * std::sin(phase);
                    set_complex(data, i, j, k, nfast, nmid, std::complex<FFT_SCALAR>(value, 0.0));
                }
            }
        }
    }

    bool is_real_input() const override { return true; }

private:
    int kx_, ky_, kz_;
    FFT_SCALAR amplitude_;
};

// =============================================================================
// SUBSECTION 2.4: GaussianGenerator
// =============================================================================

// 3D Gaussian: exp(-r²/(2σ²))
// FFT properties: Gaussian → Gaussian in frequency space
class GaussianGenerator : public DataGenerator {
public:
    // Constructor: sigma in grid units
    GaussianGenerator(FFT_SCALAR sigma = 1.0, FFT_SCALAR amplitude = 1.0) :
        sigma_(sigma), amplitude_(amplitude)
    {
    }

    void generate(FFT_SCALAR *data, int nfast, int nmid, int nslow) override
    {
        FFT_SCALAR sigma_sq = sigma_ * sigma_;

        for (int k = 0; k < nslow; k++) {
            for (int j = 0; j < nmid; j++) {
                for (int i = 0; i < nfast; i++) {
                    // Calculate distance from center (periodic boundary wrapping)
                    double dx   = std::min(i, nfast - i);
                    double dy   = std::min(j, nmid - j);
                    double dz   = std::min(k, nslow - k);
                    double r_sq = dx * dx + dy * dy + dz * dz;

                    // Gaussian: exp(-r²/(2σ²))
                    FFT_SCALAR value = amplitude_ * std::exp(-r_sq / (2.0 * sigma_sq));
                    set_complex(data, i, j, k, nfast, nmid, std::complex<FFT_SCALAR>(value, 0.0));
                }
            }
        }
    }

    bool is_real_input() const override { return true; }

private:
    FFT_SCALAR sigma_;
    FFT_SCALAR amplitude_;
};

// =============================================================================
// SUBSECTION 2.5: RandomComplexGenerator
// =============================================================================

// Random complex data (both real and imaginary parts)
// Used for round-trip testing without specific analytical properties
class RandomComplexGenerator : public DataGenerator {
public:
    RandomComplexGenerator(unsigned int seed = 12345, FFT_SCALAR max_val = 1.0) :
        seed_(seed), max_val_(max_val)
    {
    }

    void generate(FFT_SCALAR *data, int nfast, int nmid, int nslow) override
    {
        std::mt19937 rng(seed_);
        std::uniform_real_distribution<FFT_SCALAR> dist(-max_val_, max_val_);

        int n_points = nfast * nmid * nslow;
        for (int i = 0; i < n_points; i++) {
            FFT_SCALAR re = dist(rng);
            FFT_SCALAR im = dist(rng);
            set_complex_linear(data, i, std::complex<FFT_SCALAR>(re, im));
        }
    }

    bool is_real_input() const override { return false; }

private:
    unsigned int seed_;
    FFT_SCALAR max_val_;
};

// =============================================================================
// SUBSECTION 2.6: MixedModesGenerator
// =============================================================================

// Superposition of multiple sine waves (tests linearity)
// FFT(a·x + b·y) = a·FFT(x) + b·FFT(y)
class MixedModesGenerator : public DataGenerator {
public:
    struct Mode {
        int kx, ky, kz;
        FFT_SCALAR amplitude;
        FFT_SCALAR phase;

        Mode(int kx_, int ky_, int kz_, FFT_SCALAR amp = 1.0, FFT_SCALAR ph = 0.0) :
            kx(kx_), ky(ky_), kz(kz_), amplitude(amp), phase(ph)
        {
        }
    };

    MixedModesGenerator() {}

    // Add a mode to the superposition
    void add_mode(int kx, int ky, int kz, FFT_SCALAR amplitude = 1.0, FFT_SCALAR phase = 0.0)
    {
        modes_.emplace_back(kx, ky, kz, amplitude, phase);
    }

    void generate(FFT_SCALAR *data, int nfast, int nmid, int nslow) override
    {
        // Zero initialization
        int size = 2 * nfast * nmid * nslow;
        std::fill(data, data + size, 0.0);

        // Add each mode
        for (const auto &mode : modes_) {
            for (int k = 0; k < nslow; k++) {
                for (int j = 0; j < nmid; j++) {
                    for (int i = 0; i < nfast; i++) {
                        double wave_phase = TWO_PI * (mode.kx * static_cast<double>(i) / nfast +
                                                      mode.ky * static_cast<double>(j) / nmid +
                                                      mode.kz * static_cast<double>(k) / nslow);

                        // Add cos(phase) for real part, sin(phase) for imaginary part
                        // This creates cos(k·x + φ) = cos(k·x)cos(φ) - sin(k·x)sin(φ)
                        FFT_SCALAR cos_val = mode.amplitude * std::cos(wave_phase + mode.phase);

                        auto current = get_complex(data, i, j, k, nfast, nmid);
                        current += std::complex<FFT_SCALAR>(cos_val, 0.0);
                        set_complex(data, i, j, k, nfast, nmid, current);
                    }
                }
            }
        }
    }

    bool is_real_input() const override { return true; }

    const std::vector<Mode> &get_modes() const { return modes_; }

private:
    std::vector<Mode> modes_;
};

} // namespace FFTTestData

// =============================================================================
// NAMESPACE: FFTValidation
// =============================================================================

namespace FFTValidation {

using namespace FFTTestHelpers;
using namespace FFTTestData;

// =============================================================================
// SECTION 3: Validators
// =============================================================================

// Base class for all validators
class Validator {
public:
    Validator(bool verbose = false) : verbose_(verbose) {}
    virtual ~Validator() {}

    // Main validation method - returns true if validation passes
    virtual bool validate() = 0;

    // Get error statistics from last validation
    virtual const ErrorStats &get_error_stats() const { return error_stats_; }

    // Enable/disable verbose output
    void set_verbose(bool verbose) { verbose_ = verbose; }

protected:
    bool verbose_;
    ErrorStats error_stats_;
};

// =============================================================================
// SUBSECTION 3.1: RoundTripValidator
// =============================================================================

// Validates IFFT(FFT(x)) ≈ x (with proper normalization)
// Tests: Forward + backward FFT should recover original data
class RoundTripValidator : public Validator {
public:
    RoundTripValidator(const FFT_SCALAR *original_data, const FFT_SCALAR *recovered_data, int nfast,
                       int nmid, int nslow, FFT_SCALAR tolerance = ROUNDTRIP_TOLERANCE,
                       bool verbose = false) :
        Validator(verbose), original_data_(original_data), recovered_data_(recovered_data),
        nfast_(nfast), nmid_(nmid), nslow_(nslow), tolerance_(tolerance)
    {
    }

    bool validate() override
    {
        error_stats_.reset();
        int n_points = nfast_ * nmid_ * nslow_;

        FFT_SCALAR max_original = 0.0;
        for (int i = 0; i < n_points; i++) {
            auto orig_val = get_complex_linear(original_data_, i);
            max_original  = std::max(max_original, std::abs(orig_val));
        }

        // Check if data is essentially zero (avoid division by zero)
        bool is_zero_data = (max_original < 1e-12);

        for (int i = 0; i < n_points; i++) {
            auto orig = get_complex_linear(original_data_, i);
            auto recv = get_complex_linear(recovered_data_, i);

            // Calculate absolute error
            double abs_error = std::abs(orig - recv);

            // Calculate relative error (if data is not zero)
            double error = is_zero_data ? abs_error : abs_error / max_original;
            error_stats_.add(error);

            if (verbose_ && error > tolerance_) {
                std::cout << "Round-trip error at point " << i << ": " << error
                          << " (original=" << orig << ", recovered=" << recv << ")" << std::endl;
            }
        }

        bool passed = (error_stats_.max() < tolerance_);

        if (verbose_) {
            std::cout << "Round-trip validation: " << (passed ? "PASSED" : "FAILED") << std::endl;
            std::cout << "Error stats: " << error_stats_ << std::endl;
        }

        return passed;
    }

private:
    const FFT_SCALAR *original_data_;
    const FFT_SCALAR *recovered_data_;
    int nfast_, nmid_, nslow_;
    FFT_SCALAR tolerance_;
};

// =============================================================================
// SUBSECTION 3.2: KnownAnswerValidator
// =============================================================================

// Validates FFT output against known analytical result
// Example: FFT(constant) = δ(k=0), FFT(δ(x=0)) = constant
class KnownAnswerValidator : public Validator {
public:
    KnownAnswerValidator(const FFT_SCALAR *computed_fft, const FFT_SCALAR *expected_fft, int nfast,
                         int nmid, int nslow, FFT_SCALAR tolerance = KNOWN_ANSWER_TOLERANCE,
                         bool verbose = false) :
        Validator(verbose), computed_fft_(computed_fft), expected_fft_(expected_fft), nfast_(nfast),
        nmid_(nmid), nslow_(nslow), tolerance_(tolerance)
    {
    }

    bool validate() override
    {
        error_stats_.reset();
        int n_points = nfast_ * nmid_ * nslow_;

        // Find maximum expected value for relative error calculation
        FFT_SCALAR max_expected = 0.0;
        for (int i = 0; i < n_points; i++) {
            auto exp_val = get_complex_linear(expected_fft_, i);
            max_expected = std::max(max_expected, std::abs(exp_val));
        }

        bool is_zero_expected = (max_expected < 1e-12);

        for (int i = 0; i < n_points; i++) {
            auto computed = get_complex_linear(computed_fft_, i);
            auto expected = get_complex_linear(expected_fft_, i);

            double abs_error = std::abs(computed - expected);
            double error     = is_zero_expected ? abs_error : abs_error / max_expected;
            error_stats_.add(error);

            if (verbose_ && error > tolerance_) {
                std::cout << "Known answer error at point " << i << ": " << error
                          << " (computed=" << computed << ", expected=" << expected << ")"
                          << std::endl;
            }
        }

        bool passed = (error_stats_.max() < tolerance_);

        if (verbose_) {
            std::cout << "Known answer validation: " << (passed ? "PASSED" : "FAILED") << std::endl;
            std::cout << "Error stats: " << error_stats_ << std::endl;
        }

        return passed;
    }

private:
    const FFT_SCALAR *computed_fft_;
    const FFT_SCALAR *expected_fft_;
    int nfast_, nmid_, nslow_;
    FFT_SCALAR tolerance_;
};

// =============================================================================
// SUBSECTION 3.3: ParsevalValidator
// =============================================================================

// Validates Parseval's theorem: sum|x|² = (1/N³) sum|X|²
// Tests energy conservation in FFT
class ParsevalValidator : public Validator {
public:
    ParsevalValidator(const FFT_SCALAR *spatial_data, const FFT_SCALAR *frequency_data, int nfast,
                      int nmid, int nslow, FFT_SCALAR tolerance = PARSEVAL_TOLERANCE,
                      bool verbose = false) :
        Validator(verbose), spatial_data_(spatial_data), frequency_data_(frequency_data),
        nfast_(nfast), nmid_(nmid), nslow_(nslow), tolerance_(tolerance)
    {
    }

    bool validate() override
    {
        error_stats_.reset();
        int n_points = nfast_ * nmid_ * nslow_;

        // Calculate spatial energy: sum|x|²
        double spatial_energy = 0.0;
        for (int i = 0; i < n_points; i++) {
            auto val = get_complex_linear(spatial_data_, i);
            spatial_energy += std::norm(val);
        }

        // Calculate frequency energy: sum|X|²
        double frequency_energy = 0.0;
        for (int i = 0; i < n_points; i++) {
            auto val = get_complex_linear(frequency_data_, i);
            frequency_energy += std::norm(val);
        }

        // Parseval's theorem: E_spatial = (1/N³) * E_frequency
        double expected_frequency_energy = spatial_energy * n_points;
        double abs_error                 = std::abs(frequency_energy - expected_frequency_energy);

        // Relative error
        double error = (expected_frequency_energy > 1e-14) ? abs_error / expected_frequency_energy
                                                           : abs_error;

        error_stats_.add(error);

        bool passed = (error < tolerance_);

        if (verbose_) {
            std::cout << "Parseval validation: " << (passed ? "PASSED" : "FAILED") << std::endl;
            std::cout << "Spatial energy:   " << spatial_energy << std::endl;
            std::cout << "Frequency energy: " << frequency_energy << std::endl;
            std::cout << "Expected (N³×E):  " << expected_frequency_energy << std::endl;
            std::cout << "Relative error:   " << error << std::endl;
        }

        return passed;
    }

private:
    const FFT_SCALAR *spatial_data_;
    const FFT_SCALAR *frequency_data_;
    int nfast_, nmid_, nslow_;
    FFT_SCALAR tolerance_;
};

// =============================================================================
// SUBSECTION 3.4: HermitianSymmetryValidator
// =============================================================================

// Validates Hermitian symmetry: X(k) = X*(-k) for real input
// Tests: Real input → conjugate symmetry in frequency domain
class HermitianSymmetryValidator : public Validator {
public:
    HermitianSymmetryValidator(const FFT_SCALAR *frequency_data, int nfast, int nmid, int nslow,
                               FFT_SCALAR tolerance = HERMITIAN_TOLERANCE, bool verbose = false) :
        Validator(verbose), frequency_data_(frequency_data), nfast_(nfast), nmid_(nmid),
        nslow_(nslow), tolerance_(tolerance)
    {
    }

    bool validate() override
    {
        error_stats_.reset();

        // Check symmetry: X(kx,ky,kz) = conj(X(-kx,-ky,-kz))
        // Note: negative indices wrap around in FFT (periodic boundary)
        for (int k = 0; k < nslow_; k++) {
            for (int j = 0; j < nmid_; j++) {
                for (int i = 0; i < nfast_; i++) {
                    // Get X(kx, ky, kz)
                    auto val_pos = get_complex(frequency_data_, i, j, k, nfast_, nmid_);

                    // Get X(-kx, -ky, -kz) with periodic wrapping
                    int i_neg    = (i == 0) ? 0 : nfast_ - i;
                    int j_neg    = (j == 0) ? 0 : nmid_ - j;
                    int k_neg    = (k == 0) ? 0 : nslow_ - k;
                    auto val_neg = get_complex(frequency_data_, i_neg, j_neg, k_neg, nfast_, nmid_);

                    // Check if val_pos ≈ conj(val_neg)
                    auto expected    = std::conj(val_neg);
                    double abs_error = std::abs(val_pos - expected);

                    // Relative error (normalize by larger magnitude)
                    double magnitude = std::max(std::abs(val_pos), std::abs(expected));
                    double error     = (magnitude > 1e-14) ? abs_error / magnitude : abs_error;

                    error_stats_.add(error);

                    if (verbose_ && error > tolerance_) {
                        std::cout << "Hermitian symmetry error at (" << i << "," << j << "," << k
                                  << "): " << error << std::endl;
                        std::cout << "  X(k) = " << val_pos << ", X*(-k) = " << expected
                                  << std::endl;
                    }
                }
            }
        }

        bool passed = (error_stats_.max() < tolerance_);

        if (verbose_) {
            std::cout << "Hermitian symmetry validation: " << (passed ? "PASSED" : "FAILED")
                      << std::endl;
            std::cout << "Error stats: " << error_stats_ << std::endl;
        }

        return passed;
    }

private:
    const FFT_SCALAR *frequency_data_;
    int nfast_, nmid_, nslow_;
    FFT_SCALAR tolerance_;
};

// =============================================================================
// SUBSECTION 3.5: LinearityValidator
// =============================================================================

// Validates FFT linearity: FFT(a·x + b·y) = a·FFT(x) + b·FFT(y)
// Tests superposition principle
class LinearityValidator : public Validator {
public:
    LinearityValidator(const FFT_SCALAR *fft_sum,  // FFT(a·x + b·y)
                       const FFT_SCALAR *fft_x,    // FFT(x)
                       const FFT_SCALAR *fft_y,    // FFT(y)
                       FFT_SCALAR a, FFT_SCALAR b, // scaling coefficients
                       int nfast, int nmid, int nslow, FFT_SCALAR tolerance = LINEARITY_TOLERANCE,
                       bool verbose = false) :
        Validator(verbose), fft_sum_(fft_sum), fft_x_(fft_x), fft_y_(fft_y), a_(a), b_(b),
        nfast_(nfast), nmid_(nmid), nslow_(nslow), tolerance_(tolerance)
    {
    }

    bool validate() override
    {
        error_stats_.reset();
        int n_points = nfast_ * nmid_ * nslow_;

        for (int i = 0; i < n_points; i++) {
            // FFT(a·x + b·y)
            auto fft_sum = get_complex_linear(fft_sum_, i);

            // a·FFT(x) + b·FFT(y)
            auto fft_x    = get_complex_linear(fft_x_, i);
            auto fft_y    = get_complex_linear(fft_y_, i);
            auto expected = a_ * fft_x + b_ * fft_y;

            double abs_error = std::abs(fft_sum - expected);
            double magnitude = std::max(std::abs(fft_sum), std::abs(expected));
            double error     = (magnitude > 1e-14) ? abs_error / magnitude : abs_error;

            error_stats_.add(error);

            if (verbose_ && error > tolerance_) {
                std::cout << "Linearity error at point " << i << ": " << error << std::endl;
                std::cout << "  FFT(a·x+b·y) = " << fft_sum << std::endl;
                std::cout << "  a·FFT(x)+b·FFT(y) = " << expected << std::endl;
            }
        }

        bool passed = (error_stats_.max() < tolerance_);

        if (verbose_) {
            std::cout << "Linearity validation: " << (passed ? "PASSED" : "FAILED") << std::endl;
            std::cout << "Coefficients: a=" << a_ << ", b=" << b_ << std::endl;
            std::cout << "Error stats: " << error_stats_ << std::endl;
        }

        return passed;
    }

private:
    const FFT_SCALAR *fft_sum_;
    const FFT_SCALAR *fft_x_;
    const FFT_SCALAR *fft_y_;
    FFT_SCALAR a_, b_;
    int nfast_, nmid_, nslow_;
    FFT_SCALAR tolerance_;
};

} // namespace FFTValidation

#endif // LMP_FFT_TEST_HELPERS_H
