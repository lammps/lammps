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

#include "tabular_function.h"

#include "math_const.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <vector>

using namespace LAMMPS_NS;
using MathConst::MY_PI;

// =========================================================================
// TabularFunction tests
// =========================================================================

class TabularFunctionTest : public ::testing::Test {
protected:
    TabularFunction tab;
};

TEST_F(TabularFunctionTest, default_construction)
{
    EXPECT_DOUBLE_EQ(tab.get_xmax(), 0.0);
    EXPECT_DOUBLE_EQ(tab.get_xmaxsq(), 0.0);
    EXPECT_DOUBLE_EQ(tab.get_vmax(), 0.0);
}

TEST_F(TabularFunctionTest, linear_function)
{
    // Tabulate f(x) = x on [0, 10]
    const int n = 101;
    std::vector<double> values(n);
    for (int i = 0; i < n; i++) values[i] = 10.0 * i / (n - 1);

    tab.set_values(n, 0.0, 10.0, values.data());

    EXPECT_DOUBLE_EQ(tab.get_xmax(), 10.0);
    EXPECT_DOUBLE_EQ(tab.get_xmaxsq(), 100.0);
    EXPECT_DOUBLE_EQ(tab.get_vmax(), 10.0);

    // Check interpolated values
    double y, y1;
    tab.value(0.0, y, 1, y1, 1);
    EXPECT_NEAR(y, 0.0, 1e-10);
    EXPECT_NEAR(y1, 1.0, 0.05);    // derivative of x is 1

    tab.value(5.0, y, 1, y1, 1);
    EXPECT_NEAR(y, 5.0, 0.05);
    EXPECT_NEAR(y1, 1.0, 0.05);

    tab.value(10.0, y, 1, y1, 1);
    EXPECT_NEAR(y, 10.0, 0.05);
}

TEST_F(TabularFunctionTest, quadratic_function)
{
    // Tabulate f(x) = x^2 on [0, 5]
    const int n = 201;
    std::vector<double> values(n);
    double xmin = 0.0, xmax = 5.0;
    for (int i = 0; i < n; i++) {
        double x = xmin + (xmax - xmin) * i / (n - 1);
        values[i] = x * x;
    }

    tab.set_values(n, xmin, xmax, values.data());

    EXPECT_DOUBLE_EQ(tab.get_xmax(), 5.0);
    EXPECT_DOUBLE_EQ(tab.get_xmaxsq(), 25.0);
    EXPECT_DOUBLE_EQ(tab.get_vmax(), 25.0);

    // Check value and derivative at x=2: f(2)=4, f'(2)=4
    double y, y1;
    tab.value(2.0, y, 1, y1, 1);
    EXPECT_NEAR(y, 4.0, 0.01);
    EXPECT_NEAR(y1, 4.0, 0.15);    // derivative tolerance for spline

    // Check at x=3: f(3)=9, f'(3)=6
    tab.value(3.0, y, 1, y1, 1);
    EXPECT_NEAR(y, 9.0, 0.01);
    EXPECT_NEAR(y1, 6.0, 0.15);
}

TEST_F(TabularFunctionTest, constant_function)
{
    // Tabulate f(x) = 5.0 on [0, 10]
    const int n = 51;
    std::vector<double> values(n, 5.0);

    tab.set_values(n, 0.0, 10.0, values.data());

    EXPECT_DOUBLE_EQ(tab.get_vmax(), 5.0);

    double y, y1;
    tab.value(3.0, y, 1, y1, 1);
    EXPECT_NEAR(y, 5.0, 1e-10);
    EXPECT_NEAR(y1, 0.0, 1e-8);    // derivative of constant is 0
}

TEST_F(TabularFunctionTest, sinusoidal_function)
{
    // Tabulate f(x) = sin(x) on [0, pi]
    const int n = 501;
    std::vector<double> values(n);
    double xmin = 0.0, xmax = MY_PI;
    for (int i = 0; i < n; i++) {
        double x = xmin + (xmax - xmin) * i / (n - 1);
        values[i] = std::sin(x);
    }

    tab.set_values(n, xmin, xmax, values.data());

    // Check at midpoint: sin(pi/2) = 1
    double y, y1;
    tab.value(MY_PI / 2.0, y, 1, y1, 1);
    EXPECT_NEAR(y, 1.0, 0.001);
    EXPECT_NEAR(y1, 0.0, 0.02);    // cos(pi/2) = 0
}

TEST_F(TabularFunctionTest, value_only_mode)
{
    const int n = 51;
    std::vector<double> values(n);
    for (int i = 0; i < n; i++) values[i] = 2.0 * i / (n - 1);

    tab.set_values(n, 0.0, 2.0, values.data());

    // Request only value, not derivative (ny1=0)
    double y = 0.0, y1 = -999.0;
    tab.value(1.0, y, 1, y1, 0);
    EXPECT_NEAR(y, 1.0, 0.05);
    EXPECT_DOUBLE_EQ(y1, -999.0);    // y1 should not be modified
}

TEST_F(TabularFunctionTest, derivative_only_mode)
{
    const int n = 51;
    std::vector<double> values(n);
    for (int i = 0; i < n; i++) values[i] = 2.0 * i / (n - 1);

    tab.set_values(n, 0.0, 2.0, values.data());

    // Request only derivative, not value (ny=0)
    double y = -999.0, y1 = 0.0;
    tab.value(1.0, y, 0, y1, 1);
    EXPECT_DOUBLE_EQ(y, -999.0);    // y should not be modified
    EXPECT_NEAR(y1, 1.0, 0.05);
}

TEST_F(TabularFunctionTest, boundary_indices)
{
    // Test that indices are clamped to valid range, even though
    // the spline may extrapolate beyond boundary values
    const int n = 51;
    std::vector<double> values(n);
    for (int i = 0; i < n; i++) values[i] = 1.0 * i / (n - 1);

    tab.set_values(n, 0.0, 1.0, values.data());

    // Querying at boundaries should work without crashing
    double y, y1;
    tab.value(0.0, y, 1, y1, 0);
    EXPECT_NEAR(y, 0.0, 0.01);

    tab.value(1.0, y, 1, y1, 0);
    EXPECT_NEAR(y, 1.0, 0.01);

    // Values outside [xmin, xmax] should not crash (index clamping)
    tab.value(-1.0, y, 1, y1, 0);
    EXPECT_TRUE(std::isfinite(y));

    tab.value(2.0, y, 1, y1, 0);
    EXPECT_TRUE(std::isfinite(y));
}

TEST_F(TabularFunctionTest, reset_with_new_values)
{
    // First set of values
    const int n1 = 51;
    std::vector<double> values1(n1, 1.0);
    tab.set_values(n1, 0.0, 1.0, values1.data());
    EXPECT_DOUBLE_EQ(tab.get_vmax(), 1.0);

    // Reset with different values
    const int n2 = 101;
    std::vector<double> values2(n2, 3.0);
    tab.set_values(n2, 0.0, 2.0, values2.data());
    EXPECT_DOUBLE_EQ(tab.get_vmax(), 3.0);
    EXPECT_DOUBLE_EQ(tab.get_xmax(), 2.0);
}
