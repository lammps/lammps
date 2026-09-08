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
#include "library.h"
#include "region.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <cstring>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

using LAMMPS_NS::utils::split_words;

namespace LAMMPS_NS {
using ::testing::ExitedWithCode;
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
    TEST_FAILURE(".*ERROR: Unknown group command keyword: xxx.*", command("group new3 xxx"););
    TEST_FAILURE(".*ERROR: Illegal group empty command.*", command("group new3 empty xxx"););
    TEST_FAILURE(".*ERROR: Group include molecule command requires atom attribute molecule.*",
                 command("group new2 include molecule"););

    BEGIN_HIDE_OUTPUT();
    group->assign("new1 delete");
    END_HIDE_OUTPUT();
    ASSERT_EQ(group->ngroup, 6);

    TEST_FAILURE(".*ERROR: Illegal group delete command: too many arguments.*",
                 command("group new2 delete xxx"););
    TEST_FAILURE(".*ERROR: Cannot delete group all.*", command("group all delete"););
    TEST_FAILURE(".*ERROR: Could not find group delete.*", command("group new0 delete"););
    TEST_FAILURE(".*ERROR: Cannot delete group new2 currently used by fix.*",
                 command("group new2 delete"););
    TEST_FAILURE(".*ERROR: Cannot delete group new3 currently used by compute.*",
                 command("group new3 delete"););
    TEST_FAILURE(".*ERROR: Cannot delete group new4 currently used by dump.*",
                 command("group new4 delete"););
    TEST_FAILURE(".*ERROR: Cannot delete group new5 currently used by atom_modify.*",
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

    TEST_FAILURE(".*ERROR: Illegal group region command.*",
                 command("group three region left xxx"););
    TEST_FAILURE(".*ERROR: Region dummy for group region does not exist.*",
                 command("group four region dummy"););

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
    ASSERT_EQ(group->count(group->find("half"), domain->get_region_by_id("top")), 8);
    ASSERT_DOUBLE_EQ(group->mass(group->find("half"), domain->get_region_by_id("top")), 8.0);

    BEGIN_HIDE_OUTPUT();
    command("write_restart group.restart");
    command("clear");
    command("read_restart group.restart");
    platform::unlink("group.restart");
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

    TEST_FAILURE(".*ERROR: Group ID xxx does not exist.*",
                 command("group four union one two xxx"););
    TEST_FAILURE(".*ERROR: Group ID xxx does not exist.*",
                 command("group five subtract all half xxx"););
    TEST_FAILURE(".*ERROR: Group ID xxx does not exist.*",
                 command("group five intersect half top xxx"););
    delete[] flags;
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
    ASSERT_DOUBLE_EQ(group->mass(group->find("half"), domain->get_region_by_id("top")), 10);
    ASSERT_NEAR(group->charge(group->find("top")), 0, 1.0e-14);
    ASSERT_NEAR(group->charge(group->find("three"), domain->get_region_by_id("top")), 0, 1.0e-14);

    TEST_FAILURE(".*ERROR: Unknown group include keyword xxx.*",
                 command("group three include xxx"););
}

TEST_F(GroupTest, Bitmap)
{
    atomic_system();

    BEGIN_HIDE_OUTPUT();
    command("group one region left");
    command("group two region right");
    command("group three empty");
    command("group four region left");
    command("group four region right");
    command("group six subtract four one");
    END_HIDE_OUTPUT();

    int bm_one   = group->get_bitmask_by_id(FLERR, "one", "unittest 1");
    int bm_two   = group->get_bitmask_by_id(FLERR, "two", "unittest 2");
    int bm_three = group->get_bitmask_by_id(FLERR, "three", "unittest 3");
    int bm_four  = group->get_bitmask_by_id(FLERR, "four", "unittest 4");
    int bm_six   = group->get_bitmask_by_id(FLERR, "six", "unittest 6");
    int nlocal   = lmp->atom->natoms;
    auto mask    = lmp->atom->mask;

    for (int i = 0; i < nlocal; ++i) {
        if ((mask[i] & bm_one) && (mask[i] & bm_two)) {
            EXPECT_NE((mask[i] & bm_four), 0);
        }
        if (mask[i] & bm_two) {
            EXPECT_NE((mask[i] & bm_six), 0);
        }
        EXPECT_EQ((mask[i] & bm_three), 0);
    }

    TEST_FAILURE(".*ERROR: Group ID five requested by unittest 5 does not exist.*",
                 group->get_bitmask_by_id(FLERR, "five", "unittest 5"););
}

TEST_F(GroupTest, DynamicAtomic)
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
    TEST_FAILURE(".*ERROR: Cannot subtract dynamic groups.*",
                 command("group chunk subtract half grow"););
    TEST_FAILURE(".*ERROR: Cannot union groups from a dynamic group.*",
                 command("group chunk union half grow"););
    TEST_FAILURE(".*ERROR: Cannot intersect groups using a dynamic group.*",
                 command("group chunk intersect half grow"););

    BEGIN_HIDE_OUTPUT();
    command("group grow delete");
    command("variable ramp equal step");
    END_HIDE_OUTPUT();
    ASSERT_EQ(group->ngroup, 3);

    BEGIN_HIDE_OUTPUT();
    command("region chunk block -1.1 1.1 -1.1 1.1 -1.1 0.1");
    command("group chunk dynamic all region chunk every 1");
    END_HIDE_OUTPUT();
    ASSERT_EQ(group->count(group->find("all")), 64);
    ASSERT_EQ(group->count(group->find("chunk")), 0);
    ASSERT_EQ(group->ngroup, 4);
    BEGIN_HIDE_OUTPUT();
    command("run 10 post no");
    END_HIDE_OUTPUT();
    ASSERT_EQ(group->count(group->find("chunk")), 4);
    BEGIN_HIDE_OUTPUT();
    command("group chunk delete");
    command("group chunk region chunk");
    command("group within dynamic chunk region chunk every 1 within 2.0");
    command("comm_modify cutoff 4.1");
    command("run 10 post no");
    END_HIDE_OUTPUT();
    ASSERT_EQ(group->count(group->find("chunk")), 4);
    ASSERT_EQ(group->count(group->find("within")), 52);

    TEST_FAILURE(".*ERROR: Group dynamic cannot reference itself.*",
                 command("group half dynamic half region top"););
    TEST_FAILURE(".*ERROR: Group dynamic parent group dummy does not exist.*",
                 command("group half dynamic dummy region top"););

    TEST_FAILURE(".*ERROR: Variable ramp for group is invalid style.*",
                 command("group ramp variable ramp"););
    TEST_FAILURE(".*ERROR: Variable name grow for group does not exist.*",
                 command("group ramp variable grow"););
}

TEST_F(GroupTest, DynamicMolecular)
{
    molecular_system();

    BEGIN_HIDE_OUTPUT();
    command("region chunk block -1.1 1.1 -1.1 1.1 -1.1 0.1");
    command("group chunk dynamic all region chunk every 10 include molecule");
    END_HIDE_OUTPUT();
    EXPECT_EQ(group->count(group->find("all")), 64);
    EXPECT_EQ(group->count(group->find("chunk")), 0);
    ASSERT_EQ(group->ngroup, 2);
    BEGIN_HIDE_OUTPUT();
    command("run 10 post no");
    END_HIDE_OUTPUT();
    EXPECT_EQ(group->count(group->find("chunk")), 8);

    BEGIN_HIDE_OUTPUT();
    command("group chunk delete");
    command("group chunk region chunk");
    command("group within dynamic chunk within 2.0 include molecule");
    command("group exclude dynamic chunk within 2.0 exclude chunk");
    command("comm_modify cutoff 4.1");
    command("run 10 post no");
    END_HIDE_OUTPUT();
    EXPECT_EQ(group->count(group->find("chunk")), 4);
    EXPECT_EQ(group->count(group->find("within")), 59);
    EXPECT_EQ(group->count(group->find("exclude")), 48);
}

static constexpr double EPSILON = 1.0e-13;

TEST_F(GroupTest, VariableFunctions)
{
    molecular_system();

    BEGIN_HIDE_OUTPUT();
    command("group one region left");
    command("group two region right");
    command("group three empty");
    command("group four region top");
    command("set atom 5 charge $(1.0+2.0*sin(PI/32*5))");
    command("set atom 12 charge $(-1.0+2.0*sin(PI/32*12))");
    command("set atom 25 charge $(2.0+2.0*sin(PI/32*25))");
    command("set atom 32 charge $(-2.0+2.0*sin(PI/32*32))");
    command("pair_style lj/cut 5.0");
    command("pair_coeff * * 0.4 3.0");
    command("pair_coeff 2 2 0.5 3.3");
    command("pair_coeff 3 3 0.2 3.5");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    int one   = group->find("one");
    int two   = group->find("two");
    int three = group->find("three");
    int four  = group->find("four");

    auto *right = domain->get_region_by_id("right");
    auto *left  = domain->get_region_by_id("left");
    auto *top   = domain->get_region_by_id("top");

    EXPECT_EQ(group->count_all(), 64);
    EXPECT_EQ(group->count(one), 16);
    EXPECT_EQ(group->count(two), 16);
    EXPECT_EQ(group->count(three), 0);
    EXPECT_EQ(group->count(four), 16);

    EXPECT_EQ(group->count(0, right), 16);
    EXPECT_EQ(group->count(one, top), 4);
    EXPECT_EQ(group->count(two, right), 16);
    EXPECT_EQ(group->count(two, left), 0);

    EXPECT_DOUBLE_EQ(group->mass(0), 80.0);
    EXPECT_DOUBLE_EQ(group->mass(two), 32.0);
    EXPECT_DOUBLE_EQ(group->mass(three), 0.0);

    EXPECT_DOUBLE_EQ(group->mass(0, right), 32.0);
    EXPECT_DOUBLE_EQ(group->mass(one, top), 8.0);
    EXPECT_DOUBLE_EQ(group->mass(two, right), 32.0);
    EXPECT_DOUBLE_EQ(group->mass(two, left), 0.0);

    EXPECT_NEAR(group->charge(0), 0.0, EPSILON);
    EXPECT_NEAR(group->charge(two), -3.0, EPSILON);
    EXPECT_NEAR(group->charge(four), 0.0, EPSILON);

    EXPECT_NEAR(group->charge(0, right), -3.0, EPSILON);
    EXPECT_NEAR(group->charge(one, top), 0.0, EPSILON);
    EXPECT_NEAR(group->charge(two, right), -3.0, EPSILON);
    EXPECT_NEAR(group->charge(two, left), 0.0, EPSILON);

    double bounds[6];
    group->bounds(0, bounds);
    EXPECT_DOUBLE_EQ(bounds[0], -1.875);
    EXPECT_DOUBLE_EQ(bounds[1], 1.125);
    EXPECT_DOUBLE_EQ(bounds[2], -1.875);
    EXPECT_DOUBLE_EQ(bounds[3], 1.125);
    EXPECT_DOUBLE_EQ(bounds[4], -1.875);
    EXPECT_DOUBLE_EQ(bounds[5], 1.125);

    group->bounds(one, bounds);
    EXPECT_DOUBLE_EQ(bounds[0], -1.875);
    EXPECT_DOUBLE_EQ(bounds[1], -1.875);
    EXPECT_DOUBLE_EQ(bounds[2], -1.875);
    EXPECT_DOUBLE_EQ(bounds[3], 1.125);
    EXPECT_DOUBLE_EQ(bounds[4], -1.875);
    EXPECT_DOUBLE_EQ(bounds[5], 1.125);

    group->bounds(0, bounds, right);
    EXPECT_DOUBLE_EQ(bounds[0], 1.125);
    EXPECT_DOUBLE_EQ(bounds[1], 1.125);
    EXPECT_DOUBLE_EQ(bounds[2], -1.875);
    EXPECT_DOUBLE_EQ(bounds[3], 1.125);
    EXPECT_DOUBLE_EQ(bounds[4], -1.875);
    EXPECT_DOUBLE_EQ(bounds[5], 1.125);

    group->bounds(one, bounds, top);
    EXPECT_DOUBLE_EQ(bounds[0], -1.875);
    EXPECT_DOUBLE_EQ(bounds[1], -1.875);
    EXPECT_DOUBLE_EQ(bounds[2], -1.875);
    EXPECT_DOUBLE_EQ(bounds[3], -1.875);
    EXPECT_DOUBLE_EQ(bounds[4], -1.875);
    EXPECT_DOUBLE_EQ(bounds[5], 1.125);

    double center[3];
    group->xcm(0, 80.0, center);
    EXPECT_DOUBLE_EQ(center[0], -0.375);
    EXPECT_DOUBLE_EQ(center[1], -0.375);
    EXPECT_DOUBLE_EQ(center[2], -0.375);

    group->xcm(two, 32.0, center);
    EXPECT_DOUBLE_EQ(center[0], 1.125);
    EXPECT_DOUBLE_EQ(center[1], -0.375);
    EXPECT_DOUBLE_EQ(center[2], -0.375);

    group->xcm(three, 0.0, center);
    EXPECT_DOUBLE_EQ(center[0], 0);
    EXPECT_DOUBLE_EQ(center[1], 0);
    EXPECT_DOUBLE_EQ(center[2], 0);

    group->xcm(one, 8.0, center, top);
    EXPECT_DOUBLE_EQ(center[0], -1.875);
    EXPECT_DOUBLE_EQ(center[1], -1.875);
    EXPECT_DOUBLE_EQ(center[2], -0.375);

    group->vcm(one, 80.0, center);
    EXPECT_DOUBLE_EQ(center[0], 0);
    EXPECT_DOUBLE_EQ(center[1], 0);
    EXPECT_DOUBLE_EQ(center[2], 0);

    group->vcm(one, 8.0, center, top);
    EXPECT_DOUBLE_EQ(center[0], 0);
    EXPECT_DOUBLE_EQ(center[1], 0);
    EXPECT_DOUBLE_EQ(center[2], 0);

    // the total force on a group of this equilibrated system is zero up to the
    // roundoff of the sum over its atoms.  the individual values depend on the
    // summation order and thus on the accelerator package in use, so only the
    // magnitude is checked here

    constexpr double FORCE_EPSILON = 1.0e-06;

    group->fcm(0, center);
    EXPECT_NEAR(center[0], 0, FORCE_EPSILON);
    EXPECT_NEAR(center[1], 0, FORCE_EPSILON);
    EXPECT_NEAR(center[2], 0, FORCE_EPSILON);

    group->fcm(two, center);
    EXPECT_NEAR(center[0], 0, FORCE_EPSILON);
    EXPECT_NEAR(center[1], 0, FORCE_EPSILON);
    EXPECT_NEAR(center[2], 0, FORCE_EPSILON);

    group->fcm(three, center);
    EXPECT_NEAR(center[0], 0, FORCE_EPSILON);
    EXPECT_NEAR(center[1], 0, FORCE_EPSILON);
    EXPECT_NEAR(center[2], 0, FORCE_EPSILON);

    group->fcm(one, center, top);
    EXPECT_NEAR(center[0], 0, FORCE_EPSILON);
    EXPECT_NEAR(center[1], 0, FORCE_EPSILON);
    EXPECT_NEAR(center[2], 0, FORCE_EPSILON);

    EXPECT_DOUBLE_EQ(group->ke(one), 0);
    EXPECT_DOUBLE_EQ(group->ke(one, top), 0);
}

// Group::inertia must include the moment of inertia of finite-size particles
// (issue #3710). Build small standalone systems whose inertia tensor can be
// computed analytically and check Group::inertia (the method behind the
// inertia()/omega() variable functions) directly.

TEST_F(GroupTest, InertiaEllipsoid)
{
    if (!info->has_style("atom", "ellipsoid")) GTEST_SKIP();

    // two axis-aligned ellipsoids, semi-axes a=1,b=2,c=3 (set shape uses
    // diameters), m1=1 at (0,0,0), m2=3 at (5,0,0); COM at x=3.75
    // -> tensor diag (10.4, 26.75, 22.75), off-diagonals 0

    BEGIN_HIDE_OUTPUT();
    command("units lj");
    command("atom_style ellipsoid");
    command("boundary f f f");
    command("region box block -20 20 -20 20 -20 20");
    command("create_box 1 box");
    command("create_atoms 1 single 0.0 0.0 0.0 units box");
    command("create_atoms 1 single 5.0 0.0 0.0 units box");
    command("set group all shape 2.0 4.0 6.0");
    command("set atom 1 mass 1.0");
    command("set atom 2 mass 3.0");
    command("pair_style zero 5.0");
    command("pair_coeff * *");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    const int all = group->find("all");
    const double masstotal = group->mass(all);
    double xcm[3], itensor[3][3];
    group->xcm(all, masstotal, xcm);
    group->inertia(all, xcm, itensor);
    group->inertia_extended(all, itensor);

    EXPECT_NEAR(itensor[0][0], 10.4, 1.0e-12);
    EXPECT_NEAR(itensor[1][1], 26.75, 1.0e-12);
    EXPECT_NEAR(itensor[2][2], 22.75, 1.0e-12);
    EXPECT_NEAR(itensor[0][1], 0.0, 1.0e-12);
    EXPECT_NEAR(itensor[1][2], 0.0, 1.0e-12);
    EXPECT_NEAR(itensor[0][2], 0.0, 1.0e-12);
}

TEST_F(GroupTest, InertiaSphere)
{
    if (!info->has_style("atom", "sphere")) GTEST_SKIP();

    // two finite spheres of radius 1, mass 1, at (0,0,0) and (4,0,0)
    // -> tensor diag (0.8, 8.8, 8.8), off-diagonals 0

    BEGIN_HIDE_OUTPUT();
    command("units lj");
    command("atom_style sphere");
    command("boundary f f f");
    command("region box block -20 20 -20 20 -20 20");
    command("create_box 1 box");
    command("create_atoms 1 single 0.0 0.0 0.0 units box");
    command("create_atoms 1 single 4.0 0.0 0.0 units box");
    command("set group all diameter 2.0");
    command("set group all mass 1.0");
    command("pair_style zero 5.0");
    command("pair_coeff * *");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    const int all = group->find("all");
    const double masstotal = group->mass(all);
    double xcm[3], itensor[3][3];
    group->xcm(all, masstotal, xcm);
    group->inertia(all, xcm, itensor);
    group->inertia_extended(all, itensor);

    EXPECT_NEAR(itensor[0][0], 0.8, 1.0e-12);
    EXPECT_NEAR(itensor[1][1], 8.8, 1.0e-12);
    EXPECT_NEAR(itensor[2][2], 8.8, 1.0e-12);
    EXPECT_NEAR(itensor[0][1], 0.0, 1.0e-12);
    EXPECT_NEAR(itensor[1][2], 0.0, 1.0e-12);
    EXPECT_NEAR(itensor[0][2], 0.0, 1.0e-12);
}

TEST_F(GroupTest, InertiaSuperellipsoid)
{
    // atom style ellipsoid/kk does not support the superellipsoid option
    if (lmp->suffix_enable) GTEST_SKIP() << "superellipsoid is not supported with an accelerator suffix";
    if (!info->has_style("atom", "ellipsoid")) GTEST_SKIP();

    // same configuration as InertiaEllipsoid, but atom_style "ellipsoid
    // superellipsoid"; default blockiness (2,2) reduces to an ellipsoid so
    // the analytic result is identical. Exercises the bonus_super branch.
    // Mass must be set before shape (stored inertia computed at set-shape).

    BEGIN_HIDE_OUTPUT();
    command("units lj");
    command("atom_style ellipsoid superellipsoid");
    command("boundary f f f");
    command("region box block -20 20 -20 20 -20 20");
    command("create_box 1 box");
    command("create_atoms 1 single 0.0 0.0 0.0 units box");
    command("create_atoms 1 single 5.0 0.0 0.0 units box");
    command("set atom 1 mass 1.0");
    command("set atom 2 mass 3.0");
    command("set group all shape 2.0 4.0 6.0");
    command("pair_style zero 5.0");
    command("pair_coeff * *");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    const int all = group->find("all");
    const double masstotal = group->mass(all);
    double xcm[3], itensor[3][3];
    group->xcm(all, masstotal, xcm);
    group->inertia(all, xcm, itensor);
    group->inertia_extended(all, itensor);

    EXPECT_NEAR(itensor[0][0], 10.4, 1.0e-12);
    EXPECT_NEAR(itensor[1][1], 26.75, 1.0e-12);
    EXPECT_NEAR(itensor[2][2], 22.75, 1.0e-12);
    EXPECT_NEAR(itensor[0][1], 0.0, 1.0e-12);
    EXPECT_NEAR(itensor[1][2], 0.0, 1.0e-12);
    EXPECT_NEAR(itensor[0][2], 0.0, 1.0e-12);
}

TEST_F(GroupTest, InertiaBody)
{
    if (!lammps_config_has_package("BODY")) GTEST_SKIP();

    // single body/nparticle at the origin with a known diagonal inertia
    // tensor (2,3,4); Group::inertia must return it unchanged

    const char *datafile = "group_inertia_body.data";
    FILE *fp = fopen(datafile, "w");
    ASSERT_NE(fp, nullptr);
    fputs("LAMMPS body nparticle test\n\n"
          "1 atoms\n1 bodies\n1 atom types\n"
          "-10 10 xlo xhi\n-10 10 ylo yhi\n-10 10 zlo zhi\n\n"
          "Atoms\n\n"
          "1 1 1 1.0 0.0 0.0 0.0 0 0 0\n\n"
          "Bodies\n\n"
          "1 1 9\n1\n2.0 3.0 4.0 0.0 0.0 0.0\n0.0 0.0 0.0\n",
          fp);
    fclose(fp);

    BEGIN_HIDE_OUTPUT();
    command("units lj");
    command("atom_style body nparticle 1 1");
    command("read_data " + std::string(datafile));
    command("pair_style zero 5.0");
    command("pair_coeff * *");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    const int all = group->find("all");
    const double masstotal = group->mass(all);
    double xcm[3], itensor[3][3];
    group->xcm(all, masstotal, xcm);
    group->inertia(all, xcm, itensor);
    group->inertia_extended(all, itensor);

    EXPECT_NEAR(itensor[0][0], 2.0, 1.0e-12);
    EXPECT_NEAR(itensor[1][1], 3.0, 1.0e-12);
    EXPECT_NEAR(itensor[2][2], 4.0, 1.0e-12);
    EXPECT_NEAR(itensor[0][1], 0.0, 1.0e-12);
    EXPECT_NEAR(itensor[1][2], 0.0, 1.0e-12);
    EXPECT_NEAR(itensor[0][2], 0.0, 1.0e-12);

    remove(datafile);
}

TEST_F(GroupTest, InertiaRegion)
{
    if (!info->has_style("atom", "ellipsoid")) GTEST_SKIP();

    // exercise the region overload: a region containing only the ellipsoid
    // at the origin, with cm at that origin -> only that particle's own spin
    // inertia contributes: 0.2*m*(b^2+c^2,a^2+c^2,a^2+b^2) = (2.6, 2.0, 1.0)

    BEGIN_HIDE_OUTPUT();
    command("units lj");
    command("atom_style ellipsoid");
    command("boundary f f f");
    command("region box block -20 20 -20 20 -20 20");
    command("create_box 1 box");
    command("create_atoms 1 single 0.0 0.0 0.0 units box");
    command("create_atoms 1 single 5.0 0.0 0.0 units box");
    command("set group all shape 2.0 4.0 6.0");
    command("set group all mass 1.0");
    command("region only block -1.0 1.0 -1.0 1.0 -1.0 1.0 units box");
    command("pair_style zero 5.0");
    command("pair_coeff * *");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    auto *region = domain->get_region_by_id("only");
    ASSERT_NE(region, nullptr);

    const int all = group->find("all");
    double cm[3] = {0.0, 0.0, 0.0};
    double itensor[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    group->inertia(all, cm, itensor, region);
    group->inertia_extended(all, itensor, region);

    EXPECT_NEAR(itensor[0][0], 2.6, 1.0e-12);
    EXPECT_NEAR(itensor[1][1], 2.0, 1.0e-12);
    EXPECT_NEAR(itensor[2][2], 1.0, 1.0e-12);
    EXPECT_NEAR(itensor[0][1], 0.0, 1.0e-12);
    EXPECT_NEAR(itensor[1][2], 0.0, 1.0e-12);
    EXPECT_NEAR(itensor[0][2], 0.0, 1.0e-12);
}

// Group::angmom_extended adds the intrinsic (spin) angular momentum of
// finite-size particles: OMEGA-type (sphere, line) via I_spin*omega,
// ANGMOM-type (ellipsoid, superellipsoid, triangle, body) via stored angmom.

TEST_F(GroupTest, AngmomSphere)
{
    if (!info->has_style("atom", "sphere")) GTEST_SKIP();

    // two finite spheres (r=1, m=1) at rest but spinning with omega=(0,0,5);
    // orbital L is zero (v=0), spin L = 0.4*m*r^2*omega per atom -> (0,0,4)

    BEGIN_HIDE_OUTPUT();
    command("units lj");
    command("atom_style sphere");
    command("boundary f f f");
    command("region box block -20 20 -20 20 -20 20");
    command("create_box 1 box");
    command("create_atoms 1 single 0.0 0.0 0.0 units box");
    command("create_atoms 1 single 2.0 0.0 0.0 units box");
    command("set group all diameter 2.0");
    command("set group all mass 1.0");
    command("set group all omega 0.0 0.0 5.0");
    command("pair_style zero 5.0");
    command("pair_coeff * *");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    const int all = group->find("all");
    const double masstotal = group->mass(all);
    double xcm[3], lmom[3];
    group->xcm(all, masstotal, xcm);
    group->angmom(all, xcm, lmom);
    group->angmom_extended(all, lmom);

    EXPECT_NEAR(lmom[0], 0.0, 1.0e-12);
    EXPECT_NEAR(lmom[1], 0.0, 1.0e-12);
    EXPECT_NEAR(lmom[2], 4.0, 1.0e-12);
}

TEST_F(GroupTest, AngmomEllipsoid)
{
    if (!info->has_style("atom", "ellipsoid")) GTEST_SKIP();

    // a single ellipsoid at rest with its angular momentum set to (1,2,3);
    // orbital L is zero, so the group angular momentum equals the stored
    // per-particle angular momentum

    BEGIN_HIDE_OUTPUT();
    command("units lj");
    command("atom_style ellipsoid");
    command("boundary f f f");
    command("region box block -20 20 -20 20 -20 20");
    command("create_box 1 box");
    command("create_atoms 1 single 0.0 0.0 0.0 units box");
    command("set group all shape 2.0 4.0 6.0");
    command("set group all mass 1.0");
    command("set group all angmom 1.0 2.0 3.0");
    command("pair_style zero 5.0");
    command("pair_coeff * *");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    const int all = group->find("all");
    const double masstotal = group->mass(all);
    double xcm[3], lmom[3];
    group->xcm(all, masstotal, xcm);
    group->angmom(all, xcm, lmom);
    group->angmom_extended(all, lmom);

    EXPECT_NEAR(lmom[0], 1.0, 1.0e-12);
    EXPECT_NEAR(lmom[1], 2.0, 1.0e-12);
    EXPECT_NEAR(lmom[2], 3.0, 1.0e-12);
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

    // finalize the KOKKOS package explicitly: otherwise Kokkos is torn down by
    // static destructors at program exit, leading to segfaults in some cases
    // same workaround as the force-style and FFT3d test drivers

    lammps_kokkos_finalize();

    MPI_Finalize();
    return rv;
}
