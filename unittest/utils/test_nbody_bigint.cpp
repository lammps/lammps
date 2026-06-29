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

// Unit tests for nbody bigint type capacity
// Tests that bigint can handle values > INT_MAX without allocating actual rigid bodies

#include "lmptype.h"
#include "gtest/gtest.h"

#include <climits>
#include <cstdint>

using namespace LAMMPS_NS;

namespace {

// Test that bigint can store values > INT_MAX
TEST(NbodyBigint, CanStoreLargeValues)
{
    // INT_MAX = 2,147,483,647
    bigint nbody = static_cast<bigint>(INT_MAX) + 1;
    EXPECT_GT(nbody, INT_MAX);
    EXPECT_EQ(nbody, 2147483648LL);

    // Test much larger value (3 billion)
    nbody = 3000000000LL;
    EXPECT_GT(nbody, INT_MAX);
    EXPECT_EQ(nbody, 3000000000LL);
}

// Test arithmetic operations used in fix_rigid_small.cpp:3640
// double tfactor = force->mvv2e / ((6.0*nbody - nlinear) * force->boltz);
TEST(NbodyBigint, ArithmeticOperations)
{
    bigint nbody = static_cast<bigint>(INT_MAX) + 1000;
    int nlinear = 5;

    // Test the actual expression from the code
    double result = 6.0 * nbody - nlinear;

    // Expected: 6 * 2,147,484,647 - 5 = 12,884,907,877
    double expected = 6.0 * (static_cast<bigint>(INT_MAX) + 1000) - 5.0;
    EXPECT_DOUBLE_EQ(result, expected);
    EXPECT_GT(result, 6.0 * INT_MAX);  // Verify it exceeds int range
}

// Test increment operation (used in fix_rigid_small.cpp:2853)
TEST(NbodyBigint, IncrementOperation)
{
    bigint nbody = static_cast<bigint>(INT_MAX) + 100;
    bigint original = nbody;

    nbody++;
    EXPECT_EQ(nbody, original + 1);
    EXPECT_EQ(nbody, static_cast<bigint>(INT_MAX) + 101);
}

// Test comparison operations
TEST(NbodyBigint, ComparisonOperations)
{
    bigint nbody_large = static_cast<bigint>(INT_MAX) + 1;
    bigint nbody_zero = 0;
    int int_max_value = INT_MAX;

    EXPECT_GT(nbody_large, 0);
    EXPECT_GT(nbody_large, int_max_value);
    EXPECT_NE(nbody_large, int_max_value);
    EXPECT_EQ(nbody_zero, 0);
}

// Test boundary values
TEST(NbodyBigint, BoundaryValues)
{
    bigint at_limit = INT_MAX;
    bigint over_limit = static_cast<bigint>(INT_MAX) + 1;
    bigint way_over = static_cast<bigint>(INT_MAX) * 2;

    EXPECT_EQ(at_limit, 2147483647);
    EXPECT_EQ(over_limit, 2147483648LL);
    EXPECT_EQ(way_over, 4294967294LL);

    // Verify they're all different
    EXPECT_NE(at_limit, over_limit);
    EXPECT_NE(over_limit, way_over);
    EXPECT_LT(at_limit, over_limit);
    EXPECT_LT(over_limit, way_over);
}

// Test that calculations don't overflow when converting to double
TEST(NbodyBigint, DoubleConversion)
{
    bigint nbody = 3000000000LL;  // 3 billion

    // This conversion happens in compute_scalar() for temperature calculation
    double as_double = static_cast<double>(nbody);
    EXPECT_DOUBLE_EQ(as_double, 3000000000.0);

    // Test arithmetic after conversion
    double result = 6.0 * as_double;
    EXPECT_DOUBLE_EQ(result, 18000000000.0);  // 18 billion
}

// Test realistic scenario: 2.5 billion rigid bodies
TEST(NbodyBigint, RealisticLargeScale)
{
    bigint nbody = 2500000000LL;  // 2.5 billion bodies
    int nlinear = 1000;           // Some linear bodies

    // Simulate the DOF calculation from compute_scalar()
    double dof_removed = 6.0 * nbody - nlinear;

    EXPECT_GT(dof_removed, 0.0);
    EXPECT_DOUBLE_EQ(dof_removed, 15000000000.0 - 1000.0);

    // Verify we can do division (for temperature calculation)
    double temperature_factor = 1.0 / dof_removed;
    EXPECT_GT(temperature_factor, 0.0);
    EXPECT_LT(temperature_factor, 1.0);
}

}  // namespace
