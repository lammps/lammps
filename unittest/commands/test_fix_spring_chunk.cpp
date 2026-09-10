/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "../testing/core.h"

#include "atom.h"
#include "fix.h"
#include "modify.h"
#include "utils.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdlib>
#include <cstring>
#include <mpi.h>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

class FixSpringChunkTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "FixSpringChunkTest";
        LAMMPSTest::SetUp();

        BEGIN_HIDE_OUTPUT();
        command("units lj");
        command("atom_style atomic");
        command("boundary f f f");
        command("region box block -2 2 -2 2 -2 2");
        command("create_box 1 box");
        command("create_atoms 1 single 0.0 0.0 0.0 units box");
        command("mass 1 1.0");
        command("pair_style zero 1.0");
        command("pair_coeff * *");
        command("compute chunks all chunk/atom type");
        command("compute centers all com/chunk chunks");
        command("fix integrate all nve");
        command("fix spring all spring/chunk 2.0 chunks centers");
        command("run 0 post no");
        END_HIDE_OUTPUT();
    }
};

TEST_F(FixSpringChunkTest, HarmonicForceAndEnergy)
{
    ASSERT_EQ(lmp->atom->nlocal, 1);

    auto *spring = lmp->modify->get_fix_by_id("spring");
    ASSERT_NE(spring, nullptr);

    auto **force = lmp->atom->f;
    EXPECT_DOUBLE_EQ(spring->compute_scalar(), 0.0);
    EXPECT_DOUBLE_EQ(force[0][0], 0.0);
    EXPECT_DOUBLE_EQ(force[0][1], 0.0);
    EXPECT_DOUBLE_EQ(force[0][2], 0.0);

    BEGIN_HIDE_OUTPUT();
    command("displace_atoms all move 0.25 0.0 0.0 units box");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    EXPECT_DOUBLE_EQ(spring->compute_scalar(), 0.0625);
    EXPECT_DOUBLE_EQ(force[0][0], -0.5);
    EXPECT_DOUBLE_EQ(force[0][1], 0.0);
    EXPECT_DOUBLE_EQ(force[0][2], 0.0);

    // Independently verify that the applied force is the negative numerical
    // derivative of the energy returned by the fix.

    constexpr double delta = 1.0e-4;

    BEGIN_HIDE_OUTPUT();
    command("displace_atoms all move 0.0001 0.0 0.0 units box");
    command("run 0 post no");
    END_HIDE_OUTPUT();
    const double eplus = spring->compute_scalar();

    BEGIN_HIDE_OUTPUT();
    command("displace_atoms all move -0.0002 0.0 0.0 units box");
    command("run 0 post no");
    END_HIDE_OUTPUT();
    const double eminus = spring->compute_scalar();

    BEGIN_HIDE_OUTPUT();
    command("displace_atoms all move 0.0001 0.0 0.0 units box");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    const double numerical_force = -(eplus-eminus) / (2.0*delta);
    EXPECT_NEAR(force[0][0], numerical_force, 1.0e-12);

    BEGIN_HIDE_OUTPUT();
    command("displace_atoms all move 0.25 0.0 0.0 units box");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    EXPECT_DOUBLE_EQ(spring->compute_scalar(), 0.25);
    EXPECT_DOUBLE_EQ(force[0][0], -1.0);
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
            if (arg == "-v") verbose = true;
        }
    }

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
