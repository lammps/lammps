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

#include "lammps.h"

#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "math_const.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

using LAMMPS_NS::MathConst::MY_PI;
using LAMMPS_NS::utils::split_words;

namespace LAMMPS_NS {
using ::testing::ContainsRegex;
using ::testing::ExitedWithCode;
using ::testing::StrEq;

class SetTest : public LAMMPSTest {
protected:
    Atom *atom;
    Domain *domain;
    void SetUp() override
    {
        testbinary = "SetTest";
        args       = {"-log", "none", "-echo", "screen", "-nocite", "-v", "num", "1"};
        LAMMPSTest::SetUp();
        atom   = lmp->atom;
        domain = lmp->domain;
    }

    void TearDown() override { LAMMPSTest::TearDown(); }

    void atomic_system(const std::string &atom_style, const std::string units = "real")
    {
        BEGIN_HIDE_OUTPUT();
        command("atom_style " + atom_style);
        command("atom_modify map array");
        command("units " + units);
        command("lattice sc 1.0 origin 0.125 0.125 0.125");
        command("region box block 0 2 0 2 0 2");
        command("create_box 8 box");
        command("create_atoms 1 box");
        command("mass * 1.0");
        command("region left block 0.0 1.0 INF INF INF INF");
        command("region right block 1.0 2.0 INF INF INF INF");
        command("region top block INF INF 0.0 1.0 INF INF");
        command("region bottom block INF INF 1.0 2.0 INF INF");
        command("region front block INF INF INF INF 0.0 1.0");
        command("region back block INF INF INF 1.0 2.0 INF");
        command("group top region top");
        command("group bottom region bottom");
        END_HIDE_OUTPUT();
    }
};

TEST_F(SetTest, NoBoxNoAtoms)
{
    ASSERT_EQ(atom->natoms, 0);
    ASSERT_EQ(domain->box_exist, 0);
    TEST_FAILURE(".*ERROR: Set command before simulation box is.*", command("set type 1 x 0.0"););

    BEGIN_HIDE_OUTPUT();
    command("region box block 0 2 0 2 0 2");
    command("create_box 1 box");
    END_HIDE_OUTPUT();
    TEST_FAILURE(".*ERROR: Set command on system without atoms.*", command("set type 1 x 0.0"););

    BEGIN_HIDE_OUTPUT();
    command("create_atoms 1 single 0.5 0.5 0.5");
    END_HIDE_OUTPUT();
    TEST_FAILURE(".*ERROR: Illegal set command: need at least four.*", command("set type 1 x"););
    TEST_FAILURE(".*ERROR: Unknown set command style: xxx.*", command("set xxx 1 x 0.0"););
    TEST_FAILURE(".*ERROR: Set keyword or custom property yyy does not exist.*",
                 command("set type 1 yyy 0.0"););
}

TEST_F(SetTest, StylesTypes)
{
    atomic_system("molecular");
    ASSERT_EQ(atom->natoms, 8);

    BEGIN_HIDE_OUTPUT();
    command("set group all mol 1");
    command("set group top type 2");
    command("set region back type 3");
    command("set region left mol 2");
    END_HIDE_OUTPUT();
    ASSERT_EQ(atom->type[0], 2);
    ASSERT_EQ(atom->type[1], 2);
    ASSERT_EQ(atom->type[2], 1);
    ASSERT_EQ(atom->type[3], 1);
    ASSERT_EQ(atom->type[4], 2);
    ASSERT_EQ(atom->type[5], 2);
    ASSERT_EQ(atom->type[6], 1);
    ASSERT_EQ(atom->type[7], 1);

    BEGIN_HIDE_OUTPUT();
    command("set mol 1 type 4");
    command("set atom 4*7 type 5");
    END_HIDE_OUTPUT();
    ASSERT_EQ(atom->type[0], 2);
    ASSERT_EQ(atom->type[1], 4);
    ASSERT_EQ(atom->type[2], 1);
    ASSERT_EQ(atom->type[3], 5);
    ASSERT_EQ(atom->type[4], 5);
    ASSERT_EQ(atom->type[5], 5);
    ASSERT_EQ(atom->type[6], 5);
    ASSERT_EQ(atom->type[7], 4);

    BEGIN_HIDE_OUTPUT();
    command("variable rev atom 9-id");
    command("set group all type v_rev");
    END_HIDE_OUTPUT();
    ASSERT_EQ(atom->type[0], 8);
    ASSERT_EQ(atom->type[1], 7);
    ASSERT_EQ(atom->type[2], 6);
    ASSERT_EQ(atom->type[3], 5);
    ASSERT_EQ(atom->type[4], 4);
    ASSERT_EQ(atom->type[5], 3);
    ASSERT_EQ(atom->type[6], 2);
    ASSERT_EQ(atom->type[7], 1);

    BEGIN_HIDE_OUTPUT();
    command("set group all type 1");
    command("set group all type/fraction 2 0.5 453246");
    END_HIDE_OUTPUT();
    int sum = 0;
    for (int i = 0; i < 8; ++i)
        sum += (atom->type[i] == 2) ? 1 : 0;
    ASSERT_EQ(sum, 4);

    BEGIN_HIDE_OUTPUT();
    command("set group all type 1");
    command("set group all type/ratio 2 0.5 5784536");
    END_HIDE_OUTPUT();
    sum = 0;
    for (int i = 0; i < 8; ++i)
        sum += (atom->type[i] == 2) ? 1 : 0;
    ASSERT_EQ(sum, 4);

    BEGIN_HIDE_OUTPUT();
    command("set group all type 1");
    command("set group all type/subset 2 4 784536");
    END_HIDE_OUTPUT();
    sum = 0;
    for (int i = 0; i < 8; ++i)
        sum += (atom->type[i] == 2) ? 1 : 0;
    ASSERT_EQ(sum, 4);
}
} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (platform::mpi_vendor() == "Open MPI" && !LAMMPS_NS::Info::has_exceptions())
        std::cout << "Warning: using OpenMPI without exceptions. Death tests will be skipped\n";

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = split_words(var);
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
