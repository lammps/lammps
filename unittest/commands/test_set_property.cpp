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

#include "lammps.h"

#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "math_const.h"
#include "modify.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

using ::testing::ContainsRegex;
using ::testing::ExitedWithCode;
using ::testing::StrEq;

namespace LAMMPS_NS {
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

    void atomic_system(const std::string &atom_style, const std::string &units = "real")
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
    command("compute 0 all property/atom proc");
    END_HIDE_OUTPUT();
    auto compute = lmp->modify->get_compute_by_id("0");
    compute->compute_peratom();
    ASSERT_EQ(compute->vector_atom[0], 0);

    TEST_FAILURE(".*ERROR: Illegal set command: need at least four.*", command("set type 1 x"););
    TEST_FAILURE(".*ERROR: Unknown set command style: xxx.*", command("set xxx 1 x 0.0"););
    TEST_FAILURE(".*ERROR: Set keyword or custom property yyy does not exist.*",
                 command("set type 1 yyy 0.0"););

    TEST_FAILURE(".*ERROR: Cannot set attribute spin/atom for atom style atomic.*",
                 command("set atom * spin/atom 1.0 1.0 0.0 0.0"););
    TEST_FAILURE(".*ERROR: Cannot set attribute spin/atom/random for atom style atomic.*",
                 command("set atom * spin/atom/random 436273456 1.0"););

    TEST_FAILURE(".*ERROR: Illegal compute property/atom command: missing argument.*",
                 command("compute 1 all property/atom"););
    TEST_FAILURE(".*ERROR: Compute property/atom mol is not available.*",
                 command("compute 1 all property/atom mol"););
}

TEST_F(SetTest, StylesTypes)
{
    if (!Info::has_package("MOLECULE")) GTEST_SKIP();
    atomic_system("molecular");
    ASSERT_EQ(atom->natoms, 8);

    BEGIN_HIDE_OUTPUT();
    command("set group all mol 1");
    command("set group top type 2");
    command("set region back type 3");
    command("set region left mol 2");
    command("compute 1 all property/atom id type mol");
    END_HIDE_OUTPUT();

    auto compute = lmp->modify->get_compute_by_id("1");
    ASSERT_NE(compute, nullptr);
    compute->compute_peratom();

    ASSERT_EQ(atom->type[0], 2);
    ASSERT_EQ(atom->type[1], 2);
    ASSERT_EQ(atom->type[2], 1);
    ASSERT_EQ(atom->type[3], 1);
    ASSERT_EQ(atom->type[4], 2);
    ASSERT_EQ(atom->type[5], 2);
    ASSERT_EQ(atom->type[6], 1);
    ASSERT_EQ(atom->type[7], 1);

    ASSERT_EQ(atom->molecule[0], 2);
    ASSERT_EQ(atom->molecule[1], 1);
    ASSERT_EQ(atom->molecule[2], 2);
    ASSERT_EQ(atom->molecule[3], 1);
    ASSERT_EQ(atom->molecule[4], 2);
    ASSERT_EQ(atom->molecule[5], 1);
    ASSERT_EQ(atom->molecule[6], 2);
    ASSERT_EQ(atom->molecule[7], 1);

    // atom ID
    ASSERT_EQ(compute->array_atom[0][0], 1);
    ASSERT_EQ(compute->array_atom[1][0], 2);
    ASSERT_EQ(compute->array_atom[2][0], 3);
    ASSERT_EQ(compute->array_atom[3][0], 4);
    ASSERT_EQ(compute->array_atom[4][0], 5);
    ASSERT_EQ(compute->array_atom[5][0], 6);
    ASSERT_EQ(compute->array_atom[6][0], 7);
    ASSERT_EQ(compute->array_atom[7][0], 8);

    // atom type
    ASSERT_EQ(compute->array_atom[0][1], 2);
    ASSERT_EQ(compute->array_atom[1][1], 2);
    ASSERT_EQ(compute->array_atom[2][1], 1);
    ASSERT_EQ(compute->array_atom[3][1], 1);
    ASSERT_EQ(compute->array_atom[4][1], 2);
    ASSERT_EQ(compute->array_atom[5][1], 2);
    ASSERT_EQ(compute->array_atom[6][1], 1);
    ASSERT_EQ(compute->array_atom[7][1], 1);

    // mol ID
    ASSERT_EQ(compute->array_atom[0][2], 2);
    ASSERT_EQ(compute->array_atom[1][2], 1);
    ASSERT_EQ(compute->array_atom[2][2], 2);
    ASSERT_EQ(compute->array_atom[3][2], 1);
    ASSERT_EQ(compute->array_atom[4][2], 2);
    ASSERT_EQ(compute->array_atom[5][2], 1);
    ASSERT_EQ(compute->array_atom[6][2], 2);
    ASSERT_EQ(compute->array_atom[7][2], 1);

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
    command("labelmap atom 1 C 2 H");
    command("set group all type C");
    command("set group all type/fraction H 0.5 453246");
    END_HIDE_OUTPUT();
    sum = 0;
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

    BEGIN_HIDE_OUTPUT();
    command("set group all type C");
    command("set group all type/subset H 5 784536");
    END_HIDE_OUTPUT();
    sum = 0;
    for (int i = 0; i < 8; ++i)
        sum += (atom->type[i] == 2) ? 1 : 0;
    ASSERT_EQ(sum, 5);

    TEST_FAILURE(".*ERROR: Numeric index 9 is out of bounds .1-8.*", command("set type 9 x 0.0"););
    TEST_FAILURE(".*ERROR: Invalid range string: 3:10.*", command("set type 3:10 x 0.0"););
    TEST_FAILURE(".*ERROR: Could not find set group ID nope.*", command("set group nope x 0.0"););
}

TEST_F(SetTest, PosVelCharge)
{
    atomic_system("charge");
    ASSERT_EQ(atom->natoms, 8);

    BEGIN_HIDE_OUTPUT();
    command("set group top charge 1.0");
    command("set atom 5*8 charge -1.0");
    END_HIDE_OUTPUT();
    ASSERT_EQ(atom->q[0], 1);
    ASSERT_EQ(atom->q[1], 1);
    ASSERT_EQ(atom->q[2], 0);
    ASSERT_EQ(atom->q[3], 0);
    ASSERT_EQ(atom->q[4], -1);
    ASSERT_EQ(atom->q[5], -1);
    ASSERT_EQ(atom->q[6], -1);
    ASSERT_EQ(atom->q[7], -1);

    BEGIN_HIDE_OUTPUT();
    command("labelmap atom 1 C 2 H");
    command("set region right type H");
    END_HIDE_OUTPUT();

    ASSERT_EQ(atom->type[0], 1);
    ASSERT_EQ(atom->type[1], 2);
    ASSERT_EQ(atom->type[2], 1);
    ASSERT_EQ(atom->type[3], 2);
    ASSERT_EQ(atom->type[4], 1);
    ASSERT_EQ(atom->type[5], 2);
    ASSERT_EQ(atom->type[6], 1);
    ASSERT_EQ(atom->type[7], 2);

    BEGIN_HIDE_OUTPUT();
    command("set type C charge 1.25");
    command("set type H charge -1.25");
    END_HIDE_OUTPUT();

    ASSERT_EQ(atom->q[0], 1.25);
    ASSERT_EQ(atom->q[1], -1.25);
    ASSERT_EQ(atom->q[2], 1.25);
    ASSERT_EQ(atom->q[3], -1.25);
    ASSERT_EQ(atom->q[4], 1.25);
    ASSERT_EQ(atom->q[5], -1.25);
    ASSERT_EQ(atom->q[6], 1.25);
    ASSERT_EQ(atom->q[7], -1.25);

    BEGIN_HIDE_OUTPUT();
    command("variable xpos atom 0.5-x");
    command("variable ypos atom y*0.5");
    command("set atom * x v_xpos y v_ypos z 0.5");
    command("set group all vx v_xpos vy v_ypos vz 0.5");
    END_HIDE_OUTPUT();
    ASSERT_EQ(atom->x[0][0], 0.375);
    ASSERT_EQ(atom->x[0][1], 0.0625);
    ASSERT_EQ(atom->x[0][2], 0.5);
    ASSERT_EQ(atom->x[1][0], -0.625);
    ASSERT_EQ(atom->x[1][1], 0.0625);
    ASSERT_EQ(atom->x[1][2], 0.5);
    ASSERT_EQ(atom->x[2][0], 0.375);
    ASSERT_EQ(atom->x[2][1], 0.5625);
    ASSERT_EQ(atom->x[2][2], 0.5);
    ASSERT_EQ(atom->x[3][0], -0.625);
    ASSERT_EQ(atom->x[3][1], 0.5625);
    ASSERT_EQ(atom->x[3][2], 0.5);
    ASSERT_EQ(atom->x[4][0], 0.375);
    ASSERT_EQ(atom->x[4][1], 0.0625);
    ASSERT_EQ(atom->x[4][2], 0.5);
    ASSERT_EQ(atom->x[5][0], -0.625);
    ASSERT_EQ(atom->x[5][1], 0.0625);
    ASSERT_EQ(atom->x[5][2], 0.5);
    ASSERT_EQ(atom->x[6][0], 0.375);
    ASSERT_EQ(atom->x[6][1], 0.5625);
    ASSERT_EQ(atom->x[6][2], 0.5);
    ASSERT_EQ(atom->x[7][0], -0.625);
    ASSERT_EQ(atom->x[7][1], 0.5625);
    ASSERT_EQ(atom->x[7][2], 0.5);

    ASSERT_EQ(atom->v[0][0], 0.125);
    ASSERT_EQ(atom->v[0][1], 0.03125);
    ASSERT_EQ(atom->v[0][2], 0.5);
    ASSERT_EQ(atom->v[1][0], 1.125);
    ASSERT_EQ(atom->v[1][1], 0.03125);
    ASSERT_EQ(atom->v[1][2], 0.5);
    ASSERT_EQ(atom->v[2][0], 0.125);
    ASSERT_EQ(atom->v[2][1], 0.28125);
    ASSERT_EQ(atom->v[2][2], 0.5);
    ASSERT_EQ(atom->v[3][0], 1.125);
    ASSERT_EQ(atom->v[3][1], 0.28125);
    ASSERT_EQ(atom->v[3][2], 0.5);
    ASSERT_EQ(atom->v[4][0], 0.125);
    ASSERT_EQ(atom->v[4][1], 0.03125);
    ASSERT_EQ(atom->v[4][2], 0.5);
    ASSERT_EQ(atom->v[5][0], 1.125);
    ASSERT_EQ(atom->v[5][1], 0.03125);
    ASSERT_EQ(atom->v[5][2], 0.5);
    ASSERT_EQ(atom->v[6][0], 0.125);
    ASSERT_EQ(atom->v[6][1], 0.28125);
    ASSERT_EQ(atom->v[6][2], 0.5);
    ASSERT_EQ(atom->v[7][0], 1.125);
    ASSERT_EQ(atom->v[7][1], 0.28125);
    ASSERT_EQ(atom->v[7][2], 0.5);
}

TEST_F(SetTest, SpinPackage)
{
    if (!Info::has_package("SPIN")) GTEST_SKIP();
    atomic_system("spin");
    ASSERT_EQ(atom->natoms, 8);

    BEGIN_HIDE_OUTPUT();
    command("set atom 1*2 spin/atom 0.5 0.1 0.5 -0.1");
    command("set atom 8 spin/atom/random 23974 0.25");
    END_HIDE_OUTPUT();
    constexpr double vx = 0.1;
    constexpr double vy = 0.5;
    constexpr double vz = -0.1;
    const double norm   = 1.0 / sqrt(vx * vx + vy * vy + vz * vz);
    ASSERT_EQ(atom->sp[0][0], vx * norm);
    ASSERT_EQ(atom->sp[0][1], vy * norm);
    ASSERT_EQ(atom->sp[0][2], vz * norm);
    ASSERT_EQ(atom->sp[0][3], 0.5);
    ASSERT_EQ(atom->sp[1][0], vx * norm);
    ASSERT_EQ(atom->sp[1][1], vy * norm);
    ASSERT_EQ(atom->sp[1][2], vz * norm);
    ASSERT_EQ(atom->sp[1][3], 0.5);
    ASSERT_NE(atom->sp[7][0], 0.0);
    ASSERT_NE(atom->sp[7][1], 0.0);
    ASSERT_NE(atom->sp[7][2], 0.0);
    ASSERT_EQ(atom->sp[7][3], 0.25);

    for (int i = 2; i < 7; ++i)
        for (int j = 0; j < 4; ++j)
            ASSERT_EQ(atom->sp[i][j], 0);

    TEST_FAILURE(".*ERROR: Invalid spin magnitude -0.1 in set spin/atom command.*",
                 command("set atom * spin/atom -0.1 1.0 0.0 0.0"););
    TEST_FAILURE(".*ERROR: At least one spin vector component must be non-zero.*",
                 command("set atom * spin/atom 1.0 0.0 0.0 0.0"););
    TEST_FAILURE(".*ERROR: Invalid spin magnitude -0.2 in set spin/atom/random command.*",
                 command("set atom * spin/atom/random 436273456 -0.2"););
    TEST_FAILURE(".*ERROR: Invalid random number seed 0 in set spin/atom/random command.*",
                 command("set atom * spin/atom/random 0 1.0"););
}

TEST_F(SetTest, EffPackage)
{
    if (!Info::has_package("EFF")) GTEST_SKIP();
    atomic_system("electron");
    ASSERT_EQ(atom->natoms, 8);

    BEGIN_HIDE_OUTPUT();
    command("set atom 1*2 spin/electron -1");
    command("set atom 3*4 spin/electron 1");
    command("set atom 5 spin/electron 0");
    command("set atom 6 spin/electron 2");
    command("set atom 7* spin/electron 3");
    command("set region left radius/electron 0.5");
    command("set region right radius/electron 1.0");
    command("compute 2 all property/atom espin eradius");
    END_HIDE_OUTPUT();

    auto compute = lmp->modify->get_compute_by_id("2");
    ASSERT_NE(compute, nullptr);
    compute->compute_peratom();

    EXPECT_EQ(atom->spin[0], -1);
    EXPECT_EQ(atom->spin[1], -1);
    EXPECT_EQ(atom->spin[2], 1);
    EXPECT_EQ(atom->spin[3], 1);
    EXPECT_EQ(atom->spin[4], 0);
    EXPECT_EQ(atom->spin[5], 2);
    EXPECT_EQ(atom->spin[6], 3);
    EXPECT_EQ(atom->spin[7], 3);
    EXPECT_EQ(atom->eradius[0], 0.5);
    EXPECT_EQ(atom->eradius[1], 1.0);
    EXPECT_EQ(atom->eradius[2], 0.5);
    EXPECT_EQ(atom->eradius[3], 1.0);
    EXPECT_EQ(atom->eradius[4], 0.5);
    EXPECT_EQ(atom->eradius[5], 1.0);
    EXPECT_EQ(atom->eradius[6], 0.5);
    EXPECT_EQ(atom->eradius[7], 1.0);

    EXPECT_EQ(compute->array_atom[0][0], -1);
    EXPECT_EQ(compute->array_atom[1][0], -1);
    EXPECT_EQ(compute->array_atom[2][0], 1);
    EXPECT_EQ(compute->array_atom[3][0], 1);
    EXPECT_EQ(compute->array_atom[4][0], 0);
    EXPECT_EQ(compute->array_atom[5][0], 2);
    EXPECT_EQ(compute->array_atom[6][0], 3);
    EXPECT_EQ(compute->array_atom[7][0], 3);
    EXPECT_EQ(compute->array_atom[0][1], 0.5);
    EXPECT_EQ(compute->array_atom[1][1], 1.0);
    EXPECT_EQ(compute->array_atom[2][1], 0.5);
    EXPECT_EQ(compute->array_atom[3][1], 1.0);
    EXPECT_EQ(compute->array_atom[4][1], 0.5);
    EXPECT_EQ(compute->array_atom[5][1], 1.0);
    EXPECT_EQ(compute->array_atom[6][1], 0.5);
    EXPECT_EQ(compute->array_atom[7][1], 1.0);

    TEST_FAILURE(".*ERROR on proc 0: Incorrect value for electron spin: 0.5.*",
                 command("set atom * spin/electron 0.5"););
    TEST_FAILURE(".*ERROR on proc 0: Incorrect value for electron radius: -0.5.*",
                 command("set atom * radius/electron -0.5"););
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
