/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "USER-MISC/angle_fourier_simple_approx.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cmath>

using namespace LAMMPS_NS;
using ::testing::Eq;

TEST(afs_approx, cos)
{
    float angles[] = {0.0,  -0.0,   1.0,     0.25, 0.33, 3.245, 0.5,
                      -0.5, -3.245, 3.03819, 6.1,  -10,  1000};

    for (size_t i = 0; i < 13; i++) {
        float refValue = cosf(angles[i]);
        EXPECT_NEAR(refValue, fastCos(angles[i]), 3.0e-5);
    }
}

TEST(afs_approx, sin)
{
    float angles[] = {0.0,  -0.0,   1.0,     0.25, 0.33, 3.245, 0.5,
                      -0.5, -3.245, 3.03819, 6.1,  -10,  1000};

    for (size_t i = 0; i < 13; i++) {
        double refValue = sin(angles[i]);
        EXPECT_NEAR(refValue, fastSin(angles[i]), 3.0e-5);
    }
}

TEST(afs_approx, acos)
{
    double angles[] = {0.0, 0.99999, 0.25, 0.33, 0.5, -0.5, -0.99};

    for (size_t i = 0; i < 7; i++) {
        double refValue = acos(angles[i]);
        EXPECT_NEAR(refValue, fastAcos(angles[i]), 6.76e-5);
    }
}
