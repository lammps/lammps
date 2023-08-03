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
#include "domain.h"
#include "group.h"
#include "info.h"
#include "input.h"
#include "region.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstring>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {
using ::testing::ExitedWithCode;
using ::testing::StrEq;

class RegionTest : public LAMMPSTest {
protected:
    Atom *atom;
    Domain *domain;
    Group *group;

    void SetUp() override
    {
        testbinary = "RegionTest";
        LAMMPSTest::SetUp();
        atom   = lmp->atom;
        domain = lmp->domain;
        group  = lmp->group;
    }

    void atomic_system()
    {
        BEGIN_HIDE_OUTPUT();
        command("units real");
        command("lattice sc 1.0 origin 0.125 0.125 0.125");
        command("region box block -10 10 -10 10 -10 10");
        command("create_box 8 box");
        command("create_atoms 1 box");
        command("mass * 1.0");
        END_HIDE_OUTPUT();
    }
};

TEST_F(RegionTest, NoBox)
{
    auto list = domain->get_region_list();
    ASSERT_EQ(list.size(), 0);
    BEGIN_HIDE_OUTPUT();
    command("region reg1 block 0 1 0 1 0 1");
    command("region reg2 cone x 0 0 2 1 0 1 units box");
    command("region reg3 plane 0 0 0 0 0 1 side out");
    command("region reg4 prism 0 1 0 1 0 1 0.1 0.2 0.3");
    command("region reg5 sphere 0 0 0 1");
    command("region reg6 union 3 reg1 reg2 reg3");
    command("region reg7 intersect 3 reg1 reg2 reg4");
    command("region reg8 ellipsoid 0 0 0 2 1 2");
    command("region reg9 cylinder  y 0 0 1 0 1 open 1 units box");
    END_HIDE_OUTPUT();
    list = domain->get_region_list();
    EXPECT_EQ(list.size(), 9);

    auto reg = domain->get_region_by_id("reg1");
    EXPECT_EQ(reg->interior, 1);
    EXPECT_EQ(reg->scaleflag, 1);
    EXPECT_EQ(reg->bboxflag, 1);
    EXPECT_EQ(reg->varshape, 0);
    EXPECT_EQ(reg->dynamic, 0);
    EXPECT_EQ(reg->moveflag, 0);
    EXPECT_EQ(reg->rotateflag, 0);
    EXPECT_EQ(reg->openflag, 0);

    reg = domain->get_region_by_id("reg2");
    EXPECT_EQ(reg->interior, 1);
    EXPECT_EQ(reg->scaleflag, 0);
    EXPECT_EQ(reg->bboxflag, 1);
    EXPECT_EQ(reg->varshape, 0);
    EXPECT_EQ(reg->dynamic, 0);
    EXPECT_EQ(reg->moveflag, 0);
    EXPECT_EQ(reg->rotateflag, 0);
    EXPECT_EQ(reg->openflag, 0);

    reg = domain->get_region_by_id("reg3");
    EXPECT_EQ(reg->interior, 0);
    EXPECT_EQ(reg->scaleflag, 1);
    EXPECT_EQ(reg->bboxflag, 0);
    EXPECT_EQ(reg->varshape, 0);
    EXPECT_EQ(reg->dynamic, 0);
    EXPECT_EQ(reg->moveflag, 0);
    EXPECT_EQ(reg->rotateflag, 0);
    EXPECT_EQ(reg->openflag, 0);

    reg = domain->get_region_by_id("reg4");
    EXPECT_EQ(reg->interior, 1);
    EXPECT_EQ(reg->scaleflag, 1);
    EXPECT_EQ(reg->bboxflag, 1);
    EXPECT_EQ(reg->varshape, 0);
    EXPECT_EQ(reg->dynamic, 0);
    EXPECT_EQ(reg->moveflag, 0);
    EXPECT_EQ(reg->rotateflag, 0);
    EXPECT_EQ(reg->openflag, 0);

    reg = domain->get_region_by_id("reg5");
    EXPECT_EQ(reg->interior, 1);
    EXPECT_EQ(reg->scaleflag, 1);
    EXPECT_EQ(reg->bboxflag, 1);
    EXPECT_EQ(reg->varshape, 0);
    EXPECT_EQ(reg->dynamic, 0);
    EXPECT_EQ(reg->moveflag, 0);
    EXPECT_EQ(reg->rotateflag, 0);
    EXPECT_EQ(reg->openflag, 0);

    reg = domain->get_region_by_id("reg6");
    EXPECT_EQ(reg->interior, 1);
    EXPECT_EQ(reg->scaleflag, 1);
    EXPECT_EQ(reg->bboxflag, 0);
    EXPECT_EQ(reg->varshape, 0);
    EXPECT_EQ(reg->dynamic, 0);
    EXPECT_EQ(reg->moveflag, 0);
    EXPECT_EQ(reg->rotateflag, 0);
    EXPECT_EQ(reg->openflag, 0);

    reg = domain->get_region_by_id("reg7");
    EXPECT_EQ(reg->interior, 1);
    EXPECT_EQ(reg->scaleflag, 1);
    EXPECT_EQ(reg->bboxflag, 1);
    EXPECT_EQ(reg->varshape, 0);
    EXPECT_EQ(reg->dynamic, 0);
    EXPECT_EQ(reg->moveflag, 0);
    EXPECT_EQ(reg->rotateflag, 0);
    EXPECT_EQ(reg->openflag, 0);

    reg = domain->get_region_by_id("reg8");
    EXPECT_EQ(reg->interior, 1);
    EXPECT_EQ(reg->scaleflag, 1);
    EXPECT_EQ(reg->bboxflag, 1);
    EXPECT_EQ(reg->varshape, 0);
    EXPECT_EQ(reg->dynamic, 0);
    EXPECT_EQ(reg->moveflag, 0);
    EXPECT_EQ(reg->rotateflag, 0);
    EXPECT_EQ(reg->openflag, 0);

    reg = domain->get_region_by_id("reg9");
    EXPECT_EQ(reg->interior, 1);
    EXPECT_EQ(reg->scaleflag, 0);
    EXPECT_EQ(reg->bboxflag, 1);
    EXPECT_EQ(reg->varshape, 0);
    EXPECT_EQ(reg->dynamic, 0);
    EXPECT_EQ(reg->moveflag, 0);
    EXPECT_EQ(reg->rotateflag, 0);
    EXPECT_EQ(reg->openflag, 1);

    BEGIN_HIDE_OUTPUT();
    command("region reg3 delete");
    command("region reg5 delete");
    command("region reg6 delete");
    command("region reg1 delete");
    command("region reg9 delete");
    END_HIDE_OUTPUT();
    list = domain->get_region_list();
    EXPECT_EQ(list.size(), 4);

    reg = domain->get_region_by_id("reg7");
    TEST_FAILURE(".*ERROR: Region intersect region reg1 does not exist.*", reg->init(););
    TEST_FAILURE(".*ERROR: Delete region reg3 does not exist.*", command("region reg3 delete"););
}

TEST_F(RegionTest, DeathTests)
{
    atomic_system();

    auto list = domain->get_region_list();
    ASSERT_EQ(list.size(), 1);

    TEST_FAILURE(".*ERROR: Illegal region block xlo: 1 >= xhi: 0.*",
                 command("region reg1 block 1 0 0 1 0 1"););
    TEST_FAILURE(".*ERROR: Illegal region cone open face: 4.*",
                 command("region reg2 cone x 0 0 2 1 0 1 open 4"););
    TEST_FAILURE(".*ERROR: Illegal region plane normal vector: 0 0 0.*",
                 command("region reg3 plane 0 0 0 0 0 0 side out"););
    TEST_FAILURE(".*ERROR: Illegal region prism non-zero xy tilt with infinite x size.*",
                 command("region reg4 prism INF INF 0 1 0 1 0.1 0.2 0.3"););
    TEST_FAILURE(".*ERROR: Illegal region sphere radius: -1.*",
                 command("region reg5 sphere 0 0 0 -1"););
    TEST_FAILURE(".*ERROR: Illegal region ellipsoid c: -2.*",
                 command("region reg8 ellipsoid 0 0 0 2 1 -2"););
    TEST_FAILURE(".*ERROR: Illegal region cylinder axis: xx.*",
                 command("region reg9 cylinder  xx 0 0 1 0 1 open 1 units box"););

    TEST_FAILURE(".*ERROR: Unrecognized region style 'xxx'.*", command("region new1 xxx"););
    // TEST_FAILURE(".*ERROR: Illegal region command.*", command("region new1 block 0 1"););
    TEST_FAILURE(".*ERROR: Illegal region command: missing argument\\(s\\).*",
                 command("region new1 block 0 1"););
    TEST_FAILURE(".*ERROR: Delete region new3 does not exist.*", command("region new3 delete"););
}

TEST_F(RegionTest, Counts)
{
    atomic_system();

    BEGIN_HIDE_OUTPUT();
    command("region reg1 block 0 5 0 5 -5 5");
    command("region reg2 cone y 0 0 10 0 -5 5 units box");
    command("region reg3 plane 0 0 0 0 1 0 side out");
    command("region reg4 prism -5 5 -5 5 -5 5 0.1 0.2 0.3 units box");
    command("region reg5 sphere 1 0 0 6");
    command("region reg6 union 3 reg1 reg2 reg3");
    command("region reg7 intersect 3 reg1 reg2 reg4");
    command("region reg8 ellipsoid 0 0 0 5 2 3");
    command("region reg9 ellipsoid 1 0 0 6 6 6");           // same as sphere
    command("region reg10 prism 0 5 0 5 -5 5 0.0 0.0 0.0"); // same as block
    END_HIDE_OUTPUT();

    auto x     = atom->x;
    auto reg1  = domain->get_region_by_id("reg1");
    auto reg2  = domain->get_region_by_id("reg2");
    auto reg3  = domain->get_region_by_id("reg3");
    auto reg4  = domain->get_region_by_id("reg4");
    auto reg5  = domain->get_region_by_id("reg5");
    auto reg6  = domain->get_region_by_id("reg6");
    auto reg7  = domain->get_region_by_id("reg7");
    auto reg8  = domain->get_region_by_id("reg8");
    auto reg9  = domain->get_region_by_id("reg9");
    auto reg10 = domain->get_region_by_id("reg10");
    int count1, count2, count3, count4, count5, count6, count7, count8, count9, count10;
    count1 = count2 = count3 = count4 = count5 = count6 = count7 = count8 = count9 = count10 = 0;
    reg1->prematch();
    reg2->prematch();
    reg3->prematch();
    reg4->prematch();
    reg5->prematch();
    reg6->prematch();
    reg7->prematch();
    reg8->prematch();
    reg9->prematch();
    reg10->prematch();

    for (int i = 0; i < atom->nlocal; ++i) {
        if (reg1->match(x[i][0], x[i][1], x[i][2])) ++count1;
        if (reg2->match(x[i][0], x[i][1], x[i][2])) ++count2;
        if (reg3->match(x[i][0], x[i][1], x[i][2])) ++count3;
        if (reg4->match(x[i][0], x[i][1], x[i][2])) ++count4;
        if (reg5->match(x[i][0], x[i][1], x[i][2])) ++count5;
        if (reg6->match(x[i][0], x[i][1], x[i][2])) ++count6;
        if (reg7->match(x[i][0], x[i][1], x[i][2])) ++count7;
        if (reg8->match(x[i][0], x[i][1], x[i][2])) ++count8;
        if (reg9->match(x[i][0], x[i][1], x[i][2])) ++count9;
        if (reg10->match(x[i][0], x[i][1], x[i][2])) ++count10;
    }
    EXPECT_EQ(count1, 250);
    EXPECT_EQ(count2, 1132);
    EXPECT_EQ(count3, 4000);
    EXPECT_EQ(count4, 1000);
    EXPECT_EQ(count5, 907);
    EXPECT_EQ(count6, 4314);
    EXPECT_EQ(count7, 86);
    EXPECT_EQ(count8, 122);
    EXPECT_EQ(count9, count5);
    EXPECT_EQ(count10, count1);
}
} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (LAMMPS_NS::platform::mpi_vendor() == "Open MPI" && !Info::has_exceptions())
        std::cout << "Warning: using OpenMPI without exceptions. Death tests will be skipped\n";

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
