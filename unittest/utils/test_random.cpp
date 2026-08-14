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

#include "lammps.h"
#include "math_const.h"
#include "random_mars.h"
#include "random_park.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <vector>

using namespace LAMMPS_NS;
using MathConst::MY_PI2;

// =========================================================================
// Test fixture that creates a minimal LAMMPS instance for RNG tests
// =========================================================================

class RNGTest : public ::testing::Test {
protected:
    LAMMPS *lmp;

    RNGTest() : lmp(nullptr) {}

    void SetUp() override
    {
        int flag;
        MPI_Initialized(&flag);
        if (!flag) {
            int argc     = 1;
            char *args[] = {(char *)"RNGTest", nullptr};
            MPI_Init(&argc, (char ***)&args);
        }

        LAMMPS::argv args = {"RNGTest", "-log",    "none", "-echo",
                             "none",    "-screen", "none", "-nocite"};

        lmp = new LAMMPS(args, MPI_COMM_WORLD);
    }

    void TearDown() override
    {
        delete lmp;
        lmp = nullptr;
    }
};

// =========================================================================
// RanMars tests
// =========================================================================

TEST_F(RNGTest, RanMars_uniform_range)
{
    RanMars rng(lmp, 12345);
    for (int i = 0; i < 10000; i++) {
        double val = rng.uniform();
        EXPECT_GE(val, 0.0) << "uniform() produced value < 0 at iteration " << i;
        EXPECT_LE(val, 1.0) << "uniform() produced value > 1 at iteration " << i;
    }
}

TEST_F(RNGTest, RanMars_uniform_not_constant)
{
    RanMars rng(lmp, 54321);
    double first   = rng.uniform();
    bool different = false;
    for (int i = 0; i < 100; i++) {
        if (rng.uniform() != first) {
            different = true;
            break;
        }
    }
    EXPECT_TRUE(different) << "uniform() returned the same value 100 times";
}

TEST_F(RNGTest, RanMars_uniform_reproducible)
{
    std::vector<double> seq1(100), seq2(100);

    {
        RanMars rng(lmp, 42);
        for (int i = 0; i < 100; i++)
            seq1[i] = rng.uniform();
    }
    {
        RanMars rng(lmp, 42);
        for (int i = 0; i < 100; i++)
            seq2[i] = rng.uniform();
    }

    for (int i = 0; i < 100; i++) {
        EXPECT_DOUBLE_EQ(seq1[i], seq2[i]) << "Mismatch at index " << i;
    }
}

TEST_F(RNGTest, RanMars_different_seeds)
{
    RanMars rng1(lmp, 111);
    RanMars rng2(lmp, 222);

    bool different = false;
    for (int i = 0; i < 100; i++) {
        if (rng1.uniform() != rng2.uniform()) {
            different = true;
            break;
        }
    }
    EXPECT_TRUE(different) << "Different seeds produced identical sequences";
}

TEST_F(RNGTest, RanMars_gaussian_statistics)
{
    RanMars rng(lmp, 98765);
    const int N = 50000;
    double sum = 0.0, sumsq = 0.0;

    for (int i = 0; i < N; i++) {
        double val = rng.gaussian();
        sum += val;
        sumsq += val * val;
    }

    double mean     = sum / N;
    double variance = sumsq / N - mean * mean;

    // Mean should be close to 0, variance close to 1
    EXPECT_NEAR(mean, 0.0, 0.05) << "Gaussian mean should be near 0";
    EXPECT_NEAR(variance, 1.0, 0.05) << "Gaussian variance should be near 1";
}

TEST_F(RNGTest, RanMars_gaussian_with_params)
{
    RanMars rng(lmp, 11111);
    const int N = 50000;
    double mu = 5.0, sigma = 2.0;
    double sum = 0.0, sumsq = 0.0;

    for (int i = 0; i < N; i++) {
        double val = rng.gaussian(mu, sigma);
        sum += val;
        sumsq += val * val;
    }

    double mean     = sum / N;
    double variance = sumsq / N - mean * mean;

    EXPECT_NEAR(mean, mu, 0.1) << "Gaussian(mu,sigma) mean should be near mu";
    EXPECT_NEAR(variance, sigma * sigma, 0.2)
        << "Gaussian(mu,sigma) variance should be near sigma^2";
}

TEST_F(RNGTest, RanMars_rayleigh_positive)
{
    RanMars rng(lmp, 33333);
    for (int i = 0; i < 1000; i++) {
        double val = rng.rayleigh(1.0);
        EXPECT_GT(val, 0.0) << "Rayleigh should produce positive values";
    }
}

TEST_F(RNGTest, RanMars_rayleigh_statistics)
{
    RanMars rng(lmp, 44444);
    const int N  = 50000;
    double sigma = 2.0;
    double sum   = 0.0;

    for (int i = 0; i < N; i++) {
        sum += rng.rayleigh(sigma);
    }

    double mean = sum / N;
    // Rayleigh mean = sigma * sqrt(pi/2)
    double expected_mean = sigma * std::sqrt(MY_PI2);
    EXPECT_NEAR(mean, expected_mean, 0.1) << "Rayleigh mean should be sigma*sqrt(pi/2)";
}

TEST_F(RNGTest, RanMars_state_save_restore)
{
    RanMars rng(lmp, 55555);

    // Advance the RNG some steps
    for (int i = 0; i < 100; i++)
        rng.uniform();

    // Save state
    double state[103];
    rng.get_state(state);

    // Generate a sequence
    std::vector<double> seq1(50);
    for (int i = 0; i < 50; i++)
        seq1[i] = rng.uniform();

    // Restore state
    rng.set_state(state);

    // Should reproduce the same sequence
    for (int i = 0; i < 50; i++) {
        EXPECT_DOUBLE_EQ(rng.uniform(), seq1[i]) << "State restore mismatch at index " << i;
    }
}

// =========================================================================
// RanPark tests
// =========================================================================

TEST_F(RNGTest, RanPark_uniform_range)
{
    RanPark rng(lmp, 12345);
    for (int i = 0; i < 10000; i++) {
        double val = rng.uniform();
        EXPECT_GE(val, 0.0) << "uniform() produced value < 0 at iteration " << i;
        EXPECT_LE(val, 1.0) << "uniform() produced value > 1 at iteration " << i;
    }
}

TEST_F(RNGTest, RanPark_uniform_reproducible)
{
    std::vector<double> seq1(100), seq2(100);

    {
        RanPark rng(lmp, 42);
        for (int i = 0; i < 100; i++)
            seq1[i] = rng.uniform();
    }
    {
        RanPark rng(lmp, 42);
        for (int i = 0; i < 100; i++)
            seq2[i] = rng.uniform();
    }

    for (int i = 0; i < 100; i++) {
        EXPECT_DOUBLE_EQ(seq1[i], seq2[i]) << "Mismatch at index " << i;
    }
}

TEST_F(RNGTest, RanPark_different_seeds)
{
    RanPark rng1(lmp, 111);
    RanPark rng2(lmp, 222);

    bool different = false;
    for (int i = 0; i < 100; i++) {
        if (rng1.uniform() != rng2.uniform()) {
            different = true;
            break;
        }
    }
    EXPECT_TRUE(different) << "Different seeds produced identical sequences";
}

TEST_F(RNGTest, RanPark_gaussian_statistics)
{
    RanPark rng(lmp, 98765);
    const int N = 50000;
    double sum = 0.0, sumsq = 0.0;

    for (int i = 0; i < N; i++) {
        double val = rng.gaussian();
        sum += val;
        sumsq += val * val;
    }

    double mean     = sum / N;
    double variance = sumsq / N - mean * mean;

    EXPECT_NEAR(mean, 0.0, 0.05) << "Gaussian mean should be near 0";
    EXPECT_NEAR(variance, 1.0, 0.05) << "Gaussian variance should be near 1";
}

TEST_F(RNGTest, RanPark_reset_int)
{
    RanPark rng(lmp, 12345);
    // Generate some values
    for (int i = 0; i < 50; i++)
        rng.uniform();

    // Reset with a new seed
    rng.reset(67890);

    // Create a fresh RNG with same seed
    RanPark rng2(lmp, 67890);

    // They should produce the same sequence
    for (int i = 0; i < 100; i++) {
        EXPECT_DOUBLE_EQ(rng.uniform(), rng2.uniform()) << "Reset sequence mismatch at index " << i;
    }
}

TEST_F(RNGTest, RanPark_reset_with_coords)
{
    RanPark rng1(lmp, 12345);
    RanPark rng2(lmp, 99999);

    double coord[3] = {1.5, 2.5, 3.5};

    // Reset both with same seed and coords
    rng1.reset(42, coord);
    rng2.reset(42, coord);

    for (int i = 0; i < 50; i++) {
        EXPECT_DOUBLE_EQ(rng1.uniform(), rng2.uniform())
            << "Coord-reset sequence mismatch at index " << i;
    }
}

TEST_F(RNGTest, RanPark_reset_with_different_coords)
{
    RanPark rng1(lmp, 12345);
    RanPark rng2(lmp, 12345);

    double coord1[3] = {1.0, 2.0, 3.0};
    double coord2[3] = {4.0, 5.0, 6.0};

    rng1.reset(42, coord1);
    rng2.reset(42, coord2);

    bool different = false;
    for (int i = 0; i < 20; i++) {
        if (rng1.uniform() != rng2.uniform()) {
            different = true;
            break;
        }
    }
    EXPECT_TRUE(different) << "Different coords should produce different sequences";
}

TEST_F(RNGTest, RanPark_state)
{
    RanPark rng(lmp, 12345);
    int initial_state = rng.state();
    EXPECT_EQ(initial_state, 12345);

    // After generating some numbers, state should change
    rng.uniform();
    int new_state = rng.state();
    EXPECT_NE(new_state, initial_state);
    EXPECT_GT(new_state, 0);
}
