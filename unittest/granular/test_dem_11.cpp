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

// DEM verification test 11: elastic Hertzian normal impact at the particle
// impact level -- two identical spheres (and a sphere on a rigid wall) colliding
// head-on with restitution 1, verifying the Hertzian peak contact mechanics
// against the closed-form solution (Chung & Ooi, Granular Matter 13:643-656,
// 2011, Tests 1 and 2).  All test logic is shared via test_dem_common; the
// reference systems live in tests/dem11-*.yaml.

#include "test_dem_common.h"
#include "test_main.h"

#include "info.h"

#include "gtest/gtest.h"

using namespace LAMMPS_NS;

TEST(Dem11, newton_on)
{
    if (!Info::has_package("GRANULAR")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
    run_dem_trajectory_test(true, "newton on");
}

TEST(Dem11, newton_off)
{
    if (!Info::has_package("GRANULAR")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
    run_dem_trajectory_test(false, "newton off");
}
