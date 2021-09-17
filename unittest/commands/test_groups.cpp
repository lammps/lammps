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

using LAMMPS_NS::utils::split_words;

namespace LAMMPS_NS {
using ::testing::ExitedWithCode;
using ::testing::MatchesRegex;
using ::testing::StrEq;

class GroupTest : public LAMMPSTest {
protected:
    Group *group;
    Domain *domain;

    void SetUp() override
    {
        testbinary = "GroupTest";
        LAMMPSTest::SetUp();
        group  = lmp->group;
        domain = lmp->domain;
    }

    void atomic_system()
    {
        BEGIN_HIDE_OUTPUT();
        command("units real");
        command("lattice sc 1.0 origin 0.125 0.125 0.125");
        command("region box block -2 2 -2 2 -2 2");
        command("create_box 8 box");
        command("create_atoms 1 box");
        command("mass * 1.0");
        command("region left block -2.0 -1.0 INF INF INF INF");
        command("region right block 0.5  2.0 INF INF INF INF");
        command("region top block INF INF -2.0 -1.0 INF INF");
        command("set region left type 2");
        command("set region right type 3");
        END_HIDE_OUTPUT();
    }

    void molecular_system()
    {
        BEGIN_HIDE_OUTPUT();
        command("fix props all property/atom mol rmass q");
        END_HIDE_OUTPUT();

        atomic_system();

        BEGIN_HIDE_OUTPUT();
        command("variable molid atom floor(id/4)+1");
        command("variable charge atom 2.0*sin(PI/32*id)");
        command("set atom * mol v_molid");
        command("set atom * charge v_charge");
        command("set type 1 mass 0.5");
        command("set type 2*4 mass 2.0");
        END_HIDE_OUTPUT();
    }
};

TEST_F(GroupTest, NoBox)
{
    ASSERT_EQ(group->ngroup, 1);
    TEST_FAILURE(".*ERROR: Group command before simulation box.*", command("group none empty"););
}

TEST_F(GroupTest, EmptyDelete)
{
    atomic_system();

    BEGIN_HIDE_OUTPUT();
    command("group new1 empty");
    command("group new2 empty");
    command("group new2 empty");
    command("group new3 empty");
    command("group new4 empty");
    command("group new5 empty");
    command("group new6 empty");
    command("fix 1 new2 nve");
    command("compute 1 new3 ke");
    command("dump 1 new4 atom 50 dump.melt");
    command("atom_modify first new5");
    END_HIDE_OUTPUT();
    ASSERT_EQ(group->ngroup, 7);
    TEST_FAILURE(".*ERROR: Illegal group command.*", command("group new3 xxx"););
    TEST_FAILURE(".*ERROR: Illegal group command.*", command("group new3 empty xxx"););
    TEST_FAILURE(".*ERROR: Group command requires atom attribute molecule.*",
                 command("group new2 include molecule"););

    BEGIN_HIDE_OUTPUT();
    group->assign("new1 delete");
    END_HIDE_OUTPUT();
    ASSERT_EQ(group->ngroup, 6);

    TEST_FAILURE(".*ERROR: Illegal group command.*", command("group new2 delete xxx"););
    TEST_FAILURE(".*ERROR: Cannot delete group all.*", command("group all delete"););
    TEST_FAILURE(".*ERROR: Could not find group delete.*", command("group new0 delete"););
    TEST_FAILURE(".*ERROR: Cannot delete group currently used by a fix.*",
                 command("group new2 delete"););
    TEST_FAILURE(".*ERROR: Cannot delete group currently used by a compute.*",
                 command("group new3 delete"););
    TEST_FAILURE(".*ERROR: Cannot delete group currently used by a dump.*",
                 command("group new4 delete"););
    TEST_FAILURE(".*ERROR: Cannot delete group currently used by atom_modify.*",
                 command("group new5 delete"););
}

TEST_F(GroupTest, RegionClear)
{
    atomic_system();

    BEGIN_HIDE_OUTPUT();
    command("group one region left");
    command("group two region right");
    command("group three empty");
    command("group four region left");
    command("group four region right");
    END_HIDE_OUTPUT();
    ASSERT_EQ(group->count(group->find("one")), 16);
    ASSERT_EQ(group->count(group->find("two")), 16);
    ASSERT_EQ(group->count(group->find("three")), 0);
    ASSERT_EQ(group->count(group->find("four")), 32);
    ASSERT_EQ(group->count(group->find("all")), lmp->atom->natoms);
    ASSERT_EQ(group->count_all(), lmp->atom->natoms);

    TEST_FAILURE(".*ERROR: Illegal group command.*", command("group three region left xxx"););
    TEST_FAILURE(".*ERROR: Group region ID does not exist.*", command("group four region dummy"););

    BEGIN_HIDE_OUTPUT();
    command("group one clear");
    command("group two clear");
    command("group three clear");
    command("group four clear");
    END_HIDE_OUTPUT();
    ASSERT_EQ(group->count(group->find("one")), 0);
    ASSERT_EQ(group->count(group->find("two")), 0);
    ASSERT_EQ(group->count(group->find("three")), 0);
    ASSERT_EQ(group->count(group->find("four")), 0);

    BEGIN_HIDE_OUTPUT();
    command("delete_atoms region box");
    END_HIDE_OUTPUT();
    ASSERT_EQ(group->count(group->find("all")), 0);
}

TEST_F(GroupTest, SelectRestart)
{
    atomic_system();

    int *flags = new int[lmp->atom->natoms];
    for (int i = 0; i < lmp->atom->natoms; ++i)
        flags[i] = i & 1;

    BEGIN_HIDE_OUTPUT();
    command("group one region left");
    command("group two region right");
    group->create("half", flags);
    group->find_or_create("three");
    group->find_or_create("one");
    command("group four union one two");
    command("group five subtract all half four");
    command("group top region top");
    command("group six intersect half top");
    END_HIDE_OUTPUT();
    ASSERT_EQ(group->count(group->find("one")), 16);
    ASSERT_EQ(group->count(group->find("two")), 16);
    ASSERT_EQ(group->count(group->find("three")), 0);
    ASSERT_EQ(group->count(group->find("half")), 32);
    ASSERT_EQ(group->count(group->find("four")), 32);
    ASSERT_EQ(group->count(group->find("five")), 16);
    ASSERT_EQ(group->count(group->find("six")), 8);
    ASSERT_EQ(group->count(group->find("half"), domain->find_region("top")), 8);
    ASSERT_DOUBLE_EQ(group->mass(group->find("half"), domain->find_region("top")), 8.0);

    BEGIN_HIDE_OUTPUT();
    command("write_restart group.restart");
    command("clear");
    command("read_restart group.restart");
    unlink("group.restart");
    END_HIDE_OUTPUT();
    group = lmp->group;
    ASSERT_EQ(group->count(group->find("one")), 16);
    ASSERT_EQ(group->count(group->find("two")), 16);
    ASSERT_EQ(group->count(group->find("three")), 0);
    ASSERT_EQ(group->count(group->find("half")), 32);
    ASSERT_EQ(group->count(group->find("four")), 32);
    ASSERT_EQ(group->count(group->find("five")), 16);
    ASSERT_DOUBLE_EQ(group->mass(group->find("six")), 8.0);

    BEGIN_HIDE_OUTPUT();
    command("group four clear");
    command("group five clear");
    command("group six clear");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*ERROR: Group ID does not exist.*", command("group four union one two xxx"););
    TEST_FAILURE(".*ERROR: Group ID does not exist.*",
                 command("group five subtract all half xxx"););
    TEST_FAILURE(".*ERROR: Group ID does not exist.*",
                 command("group five intersect half top xxx"););
}

TEST_F(GroupTest, Molecular)
{
    molecular_system();

    BEGIN_HIDE_OUTPUT();
    command("group one region left");
    command("group two region right");
    command("group half id 1:1000:2");
    command("group top region top");
    command("group three intersect half top");
    command("group three include molecule");
    END_HIDE_OUTPUT();
    ASSERT_EQ(group->count(group->find("one")), 16);
    ASSERT_EQ(group->count(group->find("two")), 16);
    ASSERT_EQ(group->count(group->find("three")), 15);
    ASSERT_DOUBLE_EQ(group->mass(group->find("half")), 40);
    ASSERT_DOUBLE_EQ(group->mass(group->find("half"), domain->find_region("top")), 10);
    ASSERT_NEAR(group->charge(group->find("top")), 0, 1.0e-14);
    ASSERT_DOUBLE_EQ(group->charge(group->find("right"), domain->find_region("top")), 0);

    TEST_FAILURE(".*ERROR: Illegal group command.*", command("group three include xxx"););
}

TEST_F(GroupTest, Dynamic)
{
    atomic_system();

    BEGIN_HIDE_OUTPUT();
    command("variable step atom id<=step");
    command("group half id 1:1000:2");
    command("group grow dynamic half var step every 1");
    END_HIDE_OUTPUT();
    ASSERT_EQ(group->count(group->find("grow")), 0);

    BEGIN_HIDE_OUTPUT();
    command("run 10 post no");
    END_HIDE_OUTPUT();
    ASSERT_EQ(group->count(group->find("grow")), 5);

    BEGIN_HIDE_OUTPUT();
    command("group grow dynamic half var step every 1");
    command("run 10 post no");
    END_HIDE_OUTPUT();
    ASSERT_EQ(group->count(group->find("grow")), 10);

    BEGIN_HIDE_OUTPUT();
    command("group grow static");
    command("run 10 post no");
    command("group part variable step");
    END_HIDE_OUTPUT();
    ASSERT_EQ(group->count(group->find("grow")), 10);
    ASSERT_EQ(group->count(group->find("part")), 30);

    BEGIN_HIDE_OUTPUT();
    command("group grow dynamic half var step every 1");
    command("run 10 post no");
    END_HIDE_OUTPUT();
    ASSERT_EQ(group->count(group->find("grow")), 20);
    TEST_FAILURE(".*ERROR: Cannot subtract groups using a dynamic group.*",
                 command("group chunk subtract half grow"););
    TEST_FAILURE(".*ERROR: Cannot union groups using a dynamic group.*",
                 command("group chunk union half grow"););
    TEST_FAILURE(".*ERROR: Cannot intersect groups using a dynamic group.*",
                 command("group chunk intersect half grow"););

    BEGIN_HIDE_OUTPUT();
    command("group grow delete");
    command("variable ramp equal step");
    END_HIDE_OUTPUT();
    ASSERT_EQ(group->ngroup, 3);

    TEST_FAILURE(".*ERROR: Group dynamic cannot reference itself.*",
                 command("group half dynamic half region top"););
    TEST_FAILURE(".*ERROR: Group dynamic parent group does not exist.*",
                 command("group half dynamic dummy region top"););

    TEST_FAILURE(".*ERROR: Variable for group is invalid style.*",
                 command("group ramp variable ramp"););
    TEST_FAILURE(".*ERROR: Variable name for group does not exist.*",
                 command("group ramp variable grow"););
}
} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (Info::get_mpi_vendor() == "Open MPI" && !LAMMPS_NS::Info::has_exceptions())
        std::cout << "Warning: using OpenMPI without exceptions. "
                     "Death tests will be skipped\n";

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
