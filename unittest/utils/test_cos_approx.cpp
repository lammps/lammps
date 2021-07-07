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
    double angles[] = {0.0, 1.0, 0.25, 0.33, 0.5, -0.5, -3.245, -10, 1000};

    for (size_t i = 0; i < 9; i++) {
        /* code */
        double refValue = cos(angles[i]);
        EXPECT_NEAR(refValue, fastCos(angles[i]), 0.0011);
    }
}

TEST(afs_approx, acos)
{
    double angles[] = {0.0, 0.99999, 0.25, 0.33, 0.5, -0.5, -0.99};

    for (size_t i = 0; i < 7; i++) {
        /* code */
        double refValue = acos(angles[i]);
        EXPECT_NEAR(refValue, fastAcos(angles[i]), 6.76e-5);
    }
}
