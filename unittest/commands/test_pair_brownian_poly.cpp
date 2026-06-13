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

// Regression test for GitHub issue #2933: the pairwise *random* Brownian force
// of pair style brownian/poly must obey Newton's third law (be equal and
// opposite) so that linear momentum is conserved and the system is not heated
// spuriously.  The original implementation ran with "newton off" on a full
// neighbor list and applied the force only to the local atom i, drawing an
// *independent* random number for the (i,j) and (j,i) visits.  The two random
// forces were therefore unrelated rather than equal and opposite, injecting
// net momentum (and kinetic energy) every step.  The fix draws each pair's
// random force from a deterministic, order- and MPI-rank-independent RNG keyed
// on the unordered atom-tag pair and the timestep, and applies it equal and
// opposite to the two particles.  These tests assert pairwise/global force
// balance, which the buggy version violated by O(1).

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
class BrownianPolyTest : public LAMMPSTest {
protected:
    Atom *atom;

    void SetUp() override
    {
        testbinary = "BrownianPolyTest";
        LAMMPSTest::SetUp();
        if (!info->has_style("pair", "brownian/poly")) GTEST_SKIP() << "COLLOID package not enabled";
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

    // net force summed over all (local) atoms; the pure pairwise Brownian force
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

// Two unequal spheres, log (full) resistance terms enabled, pure pairwise
// (flagfld = 0).  The random force on each must be exactly opposite, at every
// timestep (the RNG seed includes the timestep).
TEST_F(BrownianPolyTest, Newton3rdLawTwoSpheres)
{
    make_box();
    BEGIN_HIDE_OUTPUT();
    command("create_atoms 1 single 0.0 0.0 0.0");
    command("create_atoms 1 single 4.5 0.0 0.0");
    command("set atom 1 diameter 2.0 density 1.0");
    command("set atom 2 diameter 6.0 density 1.0");
    // mu flaglog flagfld(=0 -> pure pairwise) cutinner cutoff t_target seed
    command("pair_style brownian/poly 1.0 1 0 4.1 6.0 1.0 5827 1 0");
    command("pair_coeff * *");
    END_HIDE_OUTPUT();

    double **f = atom->f;
    for (int step = 0; step < 5; ++step) {
        BEGIN_HIDE_OUTPUT();
        command("run 1 post no");
        END_HIDE_OUTPUT();
        const int i1 = atom->map(1);
        const int i2 = atom->map(2);
        ASSERT_GE(i1, 0);
        ASSERT_GE(i2, 0);
        // forces are O(100); the bug leaves a net force of the same order
        for (int k = 0; k < 3; ++k)
            EXPECT_NEAR(f[i1][k] + f[i2][k], 0.0, 1.0e-8)
                << "step " << step << " force component " << k;
    }
}

// Same, but with the centers inside cutinner so the gap is clamped to the
// minimum.  This is the parameter regime of the original bug report.
TEST_F(BrownianPolyTest, Newton3rdLawClampedGap)
{
    make_box();
    BEGIN_HIDE_OUTPUT();
    command("create_atoms 1 single 0.0 0.0 0.0");
    command("create_atoms 1 single 4.5 0.0 0.0");
    command("set atom 1 diameter 2.0 density 1.0");
    command("set atom 2 diameter 6.0 density 1.0");
    // cutinner = 5.5 > r = 4.5 -> the h_sep clamp branch is exercised
    command("pair_style brownian/poly 1.0 1 0 5.5 6.0 1.0 91237 1 0");
    command("pair_coeff * *");
    command("run 1 post no");
    END_HIDE_OUTPUT();

    const int i1 = atom->map(1);
    const int i2 = atom->map(2);
    ASSERT_GE(i1, 0);
    ASSERT_GE(i2, 0);
    double **f = atom->f;
    for (int k = 0; k < 3; ++k)
        EXPECT_NEAR(f[i1][k] + f[i2][k], 0.0, 1.0e-8) << "force component " << k;
}

// A cluster of polydisperse spheres: total linear momentum input from the
// pairwise random Brownian forces must vanish.  The near-field lubrication
// resistance is only physical close to contact (a_sq < 0 for large gaps, where
// the sqrt in the Brownian amplitude is undefined), so the particles are placed
// in a compact square near contact and the cutoff (3.0) excludes the longer
// diagonal pairs.
TEST_F(BrownianPolyTest, Newton3rdLawCluster)
{
    make_box();
    BEGIN_HIDE_OUTPUT();
    command("create_atoms 1 single 0.0 0.0 0.0");
    command("create_atoms 1 single 2.4 0.0 0.0");
    command("create_atoms 1 single 0.0 2.4 0.0");
    command("create_atoms 1 single 2.4 2.4 0.0");
    command("set atom 1 diameter 2.0 density 1.0");
    command("set atom 2 diameter 2.2 density 1.0");
    command("set atom 3 diameter 2.1 density 1.0");
    command("set atom 4 diameter 2.3 density 1.0");
    command("pair_style brownian/poly 1.0 1 0 0.5 3.0 1.0 40591 1 0");
    command("pair_coeff * *");
    command("run 1 post no");
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
