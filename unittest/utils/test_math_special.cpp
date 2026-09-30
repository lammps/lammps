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

#include "math_const.h"
#include "math_special.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <limits>

using namespace LAMMPS_NS;
using MathConst::MY_PI;
using MathConst::MY_PI2;

// -------------------------------------------------------------------------
// factorial()
// -------------------------------------------------------------------------

TEST(MathSpecial, factorial_zero)
{
    EXPECT_DOUBLE_EQ(MathSpecial::factorial(0), 1.0);
}

TEST(MathSpecial, factorial_small_values)
{
    EXPECT_DOUBLE_EQ(MathSpecial::factorial(1), 1.0);
    EXPECT_DOUBLE_EQ(MathSpecial::factorial(2), 2.0);
    EXPECT_DOUBLE_EQ(MathSpecial::factorial(3), 6.0);
    EXPECT_DOUBLE_EQ(MathSpecial::factorial(4), 24.0);
    EXPECT_DOUBLE_EQ(MathSpecial::factorial(5), 120.0);
    EXPECT_DOUBLE_EQ(MathSpecial::factorial(6), 720.0);
    EXPECT_DOUBLE_EQ(MathSpecial::factorial(10), 3628800.0);
}

TEST(MathSpecial, factorial_medium_values)
{
    EXPECT_DOUBLE_EQ(MathSpecial::factorial(20), 2.43290200817664e+18);
    EXPECT_DOUBLE_EQ(MathSpecial::factorial(50), 3.04140932017134e+64);
}

TEST(MathSpecial, factorial_max_valid)
{
    // 167! is the largest factorial representable as double
    double val = MathSpecial::factorial(167);
    EXPECT_FALSE(std::isnan(val));
    EXPECT_GT(val, 0.0);
    EXPECT_TRUE(std::isfinite(val));
}

TEST(MathSpecial, factorial_out_of_range)
{
    // n < 0 should return NaN
    EXPECT_TRUE(std::isnan(MathSpecial::factorial(-1)));
    EXPECT_TRUE(std::isnan(MathSpecial::factorial(-100)));
    // n > 167 should return NaN
    EXPECT_TRUE(std::isnan(MathSpecial::factorial(168)));
    EXPECT_TRUE(std::isnan(MathSpecial::factorial(1000)));
}

// -------------------------------------------------------------------------
// square(), cube()
// -------------------------------------------------------------------------

TEST(MathSpecial, square)
{
    EXPECT_DOUBLE_EQ(MathSpecial::square(0.0), 0.0);
    EXPECT_DOUBLE_EQ(MathSpecial::square(1.0), 1.0);
    EXPECT_DOUBLE_EQ(MathSpecial::square(3.0), 9.0);
    EXPECT_DOUBLE_EQ(MathSpecial::square(-4.0), 16.0);
    EXPECT_DOUBLE_EQ(MathSpecial::square(0.5), 0.25);
    EXPECT_DOUBLE_EQ(MathSpecial::square(-2.5), 6.25);
}

TEST(MathSpecial, cube)
{
    EXPECT_DOUBLE_EQ(MathSpecial::cube(0.0), 0.0);
    EXPECT_DOUBLE_EQ(MathSpecial::cube(1.0), 1.0);
    EXPECT_DOUBLE_EQ(MathSpecial::cube(2.0), 8.0);
    EXPECT_DOUBLE_EQ(MathSpecial::cube(-2.0), -8.0);
    EXPECT_DOUBLE_EQ(MathSpecial::cube(0.5), 0.125);
    EXPECT_DOUBLE_EQ(MathSpecial::cube(-3.0), -27.0);
}

// -------------------------------------------------------------------------
// powsign()
// -------------------------------------------------------------------------

TEST(MathSpecial, powsign)
{
    EXPECT_DOUBLE_EQ(MathSpecial::powsign(0), 1.0);
    EXPECT_DOUBLE_EQ(MathSpecial::powsign(1), -1.0);
    EXPECT_DOUBLE_EQ(MathSpecial::powsign(2), 1.0);
    EXPECT_DOUBLE_EQ(MathSpecial::powsign(3), -1.0);
    EXPECT_DOUBLE_EQ(MathSpecial::powsign(100), 1.0);
    EXPECT_DOUBLE_EQ(MathSpecial::powsign(101), -1.0);
    EXPECT_DOUBLE_EQ(MathSpecial::powsign(-1), -1.0);
    EXPECT_DOUBLE_EQ(MathSpecial::powsign(-2), 1.0);
}

// -------------------------------------------------------------------------
// powint()
// -------------------------------------------------------------------------

TEST(MathSpecial, powint_zero_exponent)
{
    EXPECT_DOUBLE_EQ(MathSpecial::powint(5.0, 0), 1.0);
    EXPECT_DOUBLE_EQ(MathSpecial::powint(0.0, 0), 1.0);
    EXPECT_DOUBLE_EQ(MathSpecial::powint(-3.0, 0), 1.0);
}

TEST(MathSpecial, powint_zero_base)
{
    EXPECT_DOUBLE_EQ(MathSpecial::powint(0.0, 1), 0.0);
    EXPECT_DOUBLE_EQ(MathSpecial::powint(0.0, 5), 0.0);
}

TEST(MathSpecial, powint_positive_exponents)
{
    EXPECT_DOUBLE_EQ(MathSpecial::powint(2.0, 1), 2.0);
    EXPECT_DOUBLE_EQ(MathSpecial::powint(2.0, 2), 4.0);
    EXPECT_DOUBLE_EQ(MathSpecial::powint(2.0, 3), 8.0);
    EXPECT_DOUBLE_EQ(MathSpecial::powint(2.0, 10), 1024.0);
    EXPECT_DOUBLE_EQ(MathSpecial::powint(3.0, 4), 81.0);
    EXPECT_DOUBLE_EQ(MathSpecial::powint(-2.0, 3), -8.0);
    EXPECT_DOUBLE_EQ(MathSpecial::powint(-2.0, 4), 16.0);
}

TEST(MathSpecial, powint_negative_exponents)
{
    EXPECT_DOUBLE_EQ(MathSpecial::powint(2.0, -1), 0.5);
    EXPECT_DOUBLE_EQ(MathSpecial::powint(2.0, -2), 0.25);
    EXPECT_DOUBLE_EQ(MathSpecial::powint(4.0, -1), 0.25);
    EXPECT_NEAR(MathSpecial::powint(3.0, -3), 1.0 / 27.0, 1e-15);
}

TEST(MathSpecial, powint_consistency_with_pow)
{
    // Check that powint matches std::pow for various cases
    for (int n = -5; n <= 10; n++) {
        double base = 1.7;
        EXPECT_NEAR(MathSpecial::powint(base, n), std::pow(base, n),
                    std::abs(std::pow(base, n)) * 1e-14)
            << "powint(" << base << ", " << n << ")";
    }
}

// -------------------------------------------------------------------------
// powsinxx()
// -------------------------------------------------------------------------

TEST(MathSpecial, powsinxx_zero_argument)
{
    // sin(0)/0 = 1, so (sin(0)/0)^n = 1
    EXPECT_DOUBLE_EQ(MathSpecial::powsinxx(0.0, 1), 1.0);
    EXPECT_DOUBLE_EQ(MathSpecial::powsinxx(0.0, 5), 1.0);
}

TEST(MathSpecial, powsinxx_known_values)
{
    // (sin(pi/2)/(pi/2))^1 = 2/pi
    double x        = MY_PI2;
    double expected = sin(x) / x;
    EXPECT_NEAR(MathSpecial::powsinxx(x, 1), expected, 1e-15);

    // (sin(x)/x)^2
    EXPECT_NEAR(MathSpecial::powsinxx(x, 2), expected * expected, 1e-15);
}

// -------------------------------------------------------------------------
// fm_exp() and exp2_x86()
// -------------------------------------------------------------------------

TEST(MathSpecial, fm_exp_known_values)
{
    // Compare against standard exp() within reasonable tolerance
    EXPECT_NEAR(MathSpecial::fm_exp(0.0), 1.0, 1e-14);
    EXPECT_NEAR(MathSpecial::fm_exp(1.0), std::exp(1.0), std::exp(1.0) * 1e-12);
    EXPECT_NEAR(MathSpecial::fm_exp(-1.0), std::exp(-1.0), 1e-12);
    EXPECT_NEAR(MathSpecial::fm_exp(2.0), std::exp(2.0), std::exp(2.0) * 1e-12);
}

TEST(MathSpecial, fm_exp_range_of_values)
{
    // Test a range of values match std::exp closely
    for (double x = -10.0; x <= 10.0; x += 0.5) {
        double expected = std::exp(x);
        EXPECT_NEAR(MathSpecial::fm_exp(x), expected, expected * 1e-10) << "fm_exp(" << x << ")";
    }
}

TEST(MathSpecial, fm_exp_large_negative)
{
    // Very large negative => ~0
    double val = MathSpecial::fm_exp(-800.0);
    EXPECT_GE(val, 0.0);
    EXPECT_LT(val, 1e-300);
}

TEST(MathSpecial, exp2_x86_known_values)
{
    EXPECT_NEAR(MathSpecial::exp2_x86(0.0), 1.0, 1e-14);
    EXPECT_NEAR(MathSpecial::exp2_x86(1.0), 2.0, 2.0 * 1e-12);
    EXPECT_NEAR(MathSpecial::exp2_x86(10.0), 1024.0, 1024.0 * 1e-12);
    EXPECT_NEAR(MathSpecial::exp2_x86(-1.0), 0.5, 0.5 * 1e-12);
}

// -------------------------------------------------------------------------
// expmsq() - exp(-x*x)
// -------------------------------------------------------------------------

TEST(MathSpecial, expmsq_known_values)
{
    EXPECT_NEAR(MathSpecial::expmsq(0.0), 1.0, 1e-14);
    EXPECT_NEAR(MathSpecial::expmsq(1.0), std::exp(-1.0), 1e-12);
    EXPECT_NEAR(MathSpecial::expmsq(-1.0), std::exp(-1.0), 1e-12); // symmetric
    EXPECT_NEAR(MathSpecial::expmsq(2.0), std::exp(-4.0), std::exp(-4.0) * 1e-10);
}

TEST(MathSpecial, expmsq_large_argument)
{
    // Very large x => exp(-x^2) ~ 0
    EXPECT_NEAR(MathSpecial::expmsq(40.0), 0.0, 1e-300);
    EXPECT_NEAR(MathSpecial::expmsq(-40.0), 0.0, 1e-300);
}

TEST(MathSpecial, expmsq_symmetry)
{
    // exp(-x^2) should be symmetric
    for (double x = 0.1; x < 5.0; x += 0.3) {
        EXPECT_NEAR(MathSpecial::expmsq(x), MathSpecial::expmsq(-x), 1e-15)
            << "expmsq symmetry at x=" << x;
    }
}

// -------------------------------------------------------------------------
// my_erfcx() - scaled complementary error function exp(x^2)*erfc(x)
// -------------------------------------------------------------------------

TEST(MathSpecial, my_erfcx_known_values)
{
    // erfcx(0) = erfc(0) = 1.0
    EXPECT_NEAR(MathSpecial::my_erfcx(0.0), 1.0, 1e-12);

    // For large positive x, erfcx(x) ~ 1/(x*sqrt(pi))
    double x      = 10.0;
    double approx = 1.0 / (x * std::sqrt(MY_PI));
    EXPECT_NEAR(MathSpecial::my_erfcx(x), approx, approx * 0.01);
}

TEST(MathSpecial, my_erfcx_positive_values)
{
    // erfcx should be positive for all real arguments
    for (double x = -5.0; x <= 20.0; x += 0.5) {
        EXPECT_GT(MathSpecial::my_erfcx(x), 0.0) << "erfcx(" << x << ") should be positive";
    }
}

TEST(MathSpecial, my_erfcx_monotone_positive)
{
    // For x > 0, erfcx is monotonically decreasing
    double prev = MathSpecial::my_erfcx(0.1);
    for (double x = 0.2; x <= 10.0; x += 0.1) {
        double curr = MathSpecial::my_erfcx(x);
        EXPECT_LT(curr, prev) << "erfcx should decrease for positive x at x=" << x;
        prev = curr;
    }
}

// -------------------------------------------------------------------------
// mdftaper() - MDF taper function
// -------------------------------------------------------------------------

TEST(MathSpecial, mdftaper_below_rmin)
{
    double f, df;
    MathSpecial::mdftaper(0.5, 1.0, 2.0, f, df);
    EXPECT_DOUBLE_EQ(f, 1.0);
    EXPECT_DOUBLE_EQ(df, 0.0);
}

TEST(MathSpecial, mdftaper_at_rmin)
{
    double f, df;
    MathSpecial::mdftaper(1.0, 1.0, 2.0, f, df);
    EXPECT_DOUBLE_EQ(f, 1.0);
    EXPECT_DOUBLE_EQ(df, 0.0);
}

TEST(MathSpecial, mdftaper_above_rmax)
{
    double f, df;
    MathSpecial::mdftaper(3.0, 1.0, 2.0, f, df);
    EXPECT_DOUBLE_EQ(f, 0.0);
    EXPECT_DOUBLE_EQ(df, 0.0);
}

TEST(MathSpecial, mdftaper_at_rmax)
{
    double f, df;
    MathSpecial::mdftaper(2.0, 1.0, 2.0, f, df);
    EXPECT_DOUBLE_EQ(f, 0.0);
    EXPECT_DOUBLE_EQ(df, 0.0);
}

TEST(MathSpecial, mdftaper_midpoint)
{
    double f, df;
    MathSpecial::mdftaper(1.5, 1.0, 2.0, f, df);
    // At midpoint x=0.5: f = (1-0.5)^3 * (1 + 3*0.5 + 6*0.25) = 0.125 * 4 = 0.5
    EXPECT_NEAR(f, 0.5, 1e-14);
    // f should be between 0 and 1
    EXPECT_GT(f, 0.0);
    EXPECT_LT(f, 1.0);
}

TEST(MathSpecial, mdftaper_monotone_decrease)
{
    // The taper function should monotonically decrease from rmin to rmax
    double rmin = 2.0, rmax = 5.0;
    double prev_f, prev_df;
    MathSpecial::mdftaper(rmin + 0.01, rmin, rmax, prev_f, prev_df);

    for (double r = rmin + 0.1; r < rmax; r += 0.1) {
        double f, df;
        MathSpecial::mdftaper(r, rmin, rmax, f, df);
        EXPECT_LE(f, prev_f) << "taper should decrease at r=" << r;
        EXPECT_LE(df, 0.0) << "derivative should be non-positive at r=" << r;
        prev_f = f;
    }
}
