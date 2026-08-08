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

// Regression test for GitHub issue #1933: the pairwise lubrication force of
// pair style lubricate/poly must obey Newton's third law (be equal and
// opposite) also for *polydisperse* particles.  The near-field resistance
// functions (Jeffrey & Onishi, 1984) are expansions in the symmetric gap
// xi = 2*gap/(radi+radj); an earlier implementation used the per-particle
// h_sep = gap/radi in the log-order terms and dropped a beta0 = radj/radi
// factor in the squeeze term, so the force computed on particle i (using its
// own radius as the reference length) did not match minus the force computed
// on particle j.  These tests assert pairwise/global force balance, which the
// buggy version violated by ~50% for strongly unequal radii.

#include "lammps.h"

#include "atom.h"
#include "info.h"
#include "input.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {
class LubricatePolyTest : public LAMMPSTest {
protected:
    Atom *atom;

    void SetUp() override
    {
        testbinary = "LubricatePolyTest";
        LAMMPSTest::SetUp();
        if (!info->has_style("pair", "lubricate/poly")) GTEST_SKIP() << "COLLOID package not enabled";
        atom = lmp->atom;
    }

    void TearDown() override { LAMMPSTest::TearDown(); }

    void make_box()
    {
        BEGIN_HIDE_OUTPUT();
        command("units lj");
        command("atom_style sphere");
        command("atom_modify map array");
        command("boundary p p p");
        command("newton off");
        command("comm_modify mode single vel yes");
        command("region box block -20 20 -20 20 -20 20");
        command("create_box 1 box");
        END_HIDE_OUTPUT();
    }

    // net force summed over all (local) atoms; pure pairwise interactions
    // (flagfld = 0) must conserve linear momentum, i.e. sum to zero.
    void net_force(double *fnet)
    {
        fnet[0] = fnet[1] = fnet[2] = 0.0;
        double **f = atom->f;
        for (int i = 0; i < atom->nlocal; ++i) {
            fnet[0] += f[i][0];
            fnet[1] += f[i][1];
            fnet[2] += f[i][2];
        }
    }
};

// Two unequal spheres, general relative translation + spin, in the regular
// (unclamped) lubrication regime.  The force on each must be exactly opposite.
TEST_F(LubricatePolyTest, Newton3rdLawTwoSpheres)
{
    make_box();
    BEGIN_HIDE_OUTPUT();
    command("create_atoms 1 single 0.0 0.0 0.0");
    command("create_atoms 1 single 4.5 0.0 0.0");
    command("set atom 1 diameter 2.0 density 1.0");
    command("set atom 2 diameter 6.0 density 1.0");
    command("set atom 1 vx 0.5 vy 0.3 vz -0.1");
    command("set atom 2 vx -0.4 vy 0.2 vz 0.15");
    command("set atom 1 omega 0.0 0.0 1.0");
    command("set atom 2 omega 0.0 0.0 -0.5");
    // mu flaglog flagfld(=0 -> pure pairwise) cutinner cutoff
    command("pair_style lubricate/poly 1.0 1 0 4.1 6.0");
    command("pair_coeff * *");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    const int i1 = atom->map(1);
    const int i2 = atom->map(2);
    ASSERT_GE(i1, 0);
    ASSERT_GE(i2, 0);
    double **f = atom->f;
    // forces are O(20); the fixed code balances to ~1e-14, the bug to ~1e-1
    for (int k = 0; k < 3; ++k)
        EXPECT_NEAR(f[i1][k] + f[i2][k], 0.0, 1.0e-8) << "force component " << k;
}

// Same, but with the centers inside cutinner so the gap is clamped.  This is
// the parameter regime of the original bug report and must be symmetric too.
TEST_F(LubricatePolyTest, Newton3rdLawClampedGap)
{
    make_box();
    BEGIN_HIDE_OUTPUT();
    command("create_atoms 1 single 0.0 0.0 0.0");
    command("create_atoms 1 single 4.5 0.0 0.0");
    command("set atom 1 diameter 2.0 density 1.0");
    command("set atom 2 diameter 6.0 density 1.0");
    command("set atom 1 vx 0.5 vy 0.0 vz 0.0");
    command("set atom 2 vx -0.5 vy 0.0 vz 0.0");
    // cutinner = 5.5 > r = 4.5 -> the h_sep clamp branch is exercised
    command("pair_style lubricate/poly 1.0 1 0 5.5 6.0");
    command("pair_coeff * *");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    const int i1 = atom->map(1);
    const int i2 = atom->map(2);
    ASSERT_GE(i1, 0);
    ASSERT_GE(i2, 0);
    double **f = atom->f;
    for (int k = 0; k < 3; ++k)
        EXPECT_NEAR(f[i1][k] + f[i2][k], 0.0, 1.0e-8) << "force component " << k;
}

// A cluster of several polydisperse spheres: total linear momentum input from
// the pairwise lubrication forces must vanish.
TEST_F(LubricatePolyTest, Newton3rdLawCluster)
{
    make_box();
    BEGIN_HIDE_OUTPUT();
    command("create_atoms 1 single  0.0  0.0  0.0");
    command("create_atoms 1 single  3.6  0.4 -0.3");
    command("create_atoms 1 single -2.9  2.7  0.6");
    command("create_atoms 1 single  1.1 -3.2  2.4");
    command("create_atoms 1 single -1.5 -1.2 -3.0");
    command("set atom 1 diameter 1.0 density 1.0");
    command("set atom 2 diameter 5.0 density 1.0");
    command("set atom 3 diameter 2.4 density 1.0");
    command("set atom 4 diameter 3.7 density 1.0");
    command("set atom 5 diameter 1.8 density 1.0");
    command("velocity all create 1.5 5934 mom no rot no");
    command("pair_style lubricate/poly 1.0 1 0 0.9 6.0");
    command("pair_coeff * *");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    double fnet[3];
    net_force(fnet);
    for (int k = 0; k < 3; ++k) EXPECT_NEAR(fnet[k], 0.0, 1.0e-8) << "net force component " << k;
}
} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = LAMMPS_NS::utils::split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }
    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
