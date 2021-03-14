/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
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
#include "regrion.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstring>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

#if defined(OMPI_MAJOR_VERSION)
const bool have_openmpi = true;
#else
const bool have_openmpi = false;
#endif

using LAMMPS_NS::utils::split_words;

namespace LAMMPS_NS {
using ::testing::ExitedWithCode;
using ::testing::MatchesRegex;
using ::testing::StrEq;

#define TEST_FAILURE(errmsg, ...)                                 \
    if (Info::has_exceptions()) {                                 \
        ::testing::internal::CaptureStdout();                     \
        ASSERT_ANY_THROW({__VA_ARGS__});                          \
        auto mesg = ::testing::internal::GetCapturedStdout();     \
        if (verbose) std::cout << mesg;                           \
        ASSERT_THAT(mesg, MatchesRegex(errmsg));                  \
    } else {                                                      \
        if (!have_openmpi) {                                      \
            ::testing::internal::CaptureStdout();                 \
            ASSERT_DEATH({__VA_ARGS__}, "");                      \
            auto mesg = ::testing::internal::GetCapturedStdout(); \
            if (verbose) std::cout << mesg;                       \
            ASSERT_THAT(mesg, MatchesRegex(errmsg));              \
        }                                                         \
    }

class GroupTest : public ::testing::Test {
protected:
    LAMMPS *lmp;

    void SetUp() override
    {
        const char *args[] = {"GroupTest", "-log", "none", "-echo", "screen", "-nocite"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void TearDown() override
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        delete lmp;
        if (!verbose) ::testing::internal::GetCapturedStdout();
        std::cout.flush();
    }

    void atomic_system()
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp->input->one("units real");
        lmp->input->one("lattice sc 1.0 origin 0.125 0.125 0.125");
        lmp->input->one("region box block -2 2 -2 2 -2 2");
        lmp->input->one("create_box 8 box");
        lmp->input->one("create_atoms 1 box");
        lmp->input->one("mass * 1.0");
        lmp->input->one("region left block -2.0 -1.0 INF INF INF INF");
        lmp->input->one("region right block 0.5  2.0 INF INF INF INF");
        lmp->input->one("set region left type 2");
        lmp->input->one("set region right type 3");
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void molecular_system()
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp->input->one("atom_style full");
        atomic_system();
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }
};

TEST_F(GroupTest, NoBox)
{
    ASSERT_EQ(lmp->group->ngroup, 1);
    TEST_FAILURE(".*ERROR: Group command before simulation box.*",
                 lmp->input->one("group none empty"););
}

TEST_F(GroupTest, EmptyDelete)
{
    atomic_system();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("group new1 empty");
    lmp->input->one("group new2 empty");
    lmp->input->one("group new2 empty");
    lmp->input->one("group new3 empty");
    lmp->input->one("group new4 empty");
    lmp->input->one("group new5 empty");
    lmp->input->one("group new6 empty");
    lmp->input->one("fix 1 new2 nve");
    lmp->input->one("compute 1 new3 ke");
    lmp->input->one("dump 1 new4 atom 50 dump.melt");
    lmp->input->one("atom_modify first new5");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->group->ngroup, 7);
    TEST_FAILURE(".*ERROR: Illegal group command.*", lmp->input->one("group new3 empty xxx"););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->group->assign("new1 delete");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->group->ngroup, 6);

    TEST_FAILURE(".*ERROR: Illegal group command.*", lmp->input->one("group new2 delete xxx"););
    TEST_FAILURE(".*ERROR: Cannot delete group all.*", lmp->input->one("group all delete"););
    TEST_FAILURE(".*ERROR: Could not find group delete.*", lmp->input->one("group new0 delete"););
    TEST_FAILURE(".*ERROR: Cannot delete group currently used by a fix.*",
                 lmp->input->one("group new2 delete"););
    TEST_FAILURE(".*ERROR: Cannot delete group currently used by a compute.*",
                 lmp->input->one("group new3 delete"););
    TEST_FAILURE(".*ERROR: Cannot delete group currently used by a dump.*",
                 lmp->input->one("group new4 delete"););
    TEST_FAILURE(".*ERROR: Cannot delete group currently used by atom_modify.*",
                 lmp->input->one("group new5 delete"););
}

TEST_F(GroupTest, RegionClear)
{
    atomic_system();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("group one region left");
    lmp->input->one("group two region right");
    lmp->input->one("group three empty");
    lmp->input->one("group four region left");
    lmp->input->one("group four region right");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->group->count(lmp->group->find("one")), 16);
    ASSERT_EQ(lmp->group->count(lmp->group->find("two")), 16);
    ASSERT_EQ(lmp->group->count(lmp->group->find("three")), 0);
    ASSERT_EQ(lmp->group->count(lmp->group->find("four")), 32);
    ASSERT_EQ(lmp->group->count(lmp->group->find("all")), lmp->atom->natoms);
    ASSERT_EQ(lmp->group->count_all(), lmp->atom->natoms);

    TEST_FAILURE(".*ERROR: Illegal group command.*",
                 lmp->input->one("group three region left xxx"););
    TEST_FAILURE(".*ERROR: Group region ID does not exist.*",
                 lmp->input->one("group four region top"););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("group one clear");
    lmp->input->one("group two clear");
    lmp->input->one("group three clear");
    lmp->input->one("group four clear");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->group->count(lmp->group->find("one")), 0);
    ASSERT_EQ(lmp->group->count(lmp->group->find("two")), 0);
    ASSERT_EQ(lmp->group->count(lmp->group->find("three")), 0);
    ASSERT_EQ(lmp->group->count(lmp->group->find("four")), 0);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("delete_atoms region box");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->group->count(lmp->group->find("all")), 0);
}

TEST_F(GroupTest, SelectRestart)
{
    atomic_system();

    int *flags = new int[lmp->atom->natoms];
    for (int i = 0; i < lmp->atom->natoms; ++i)
        flags[i] = i & 1;

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("group one region left");
    lmp->input->one("group two region right");
    lmp->group->create("half", flags);
    lmp->group->find_or_create("three");
    lmp->group->find_or_create("one");
    lmp->input->one("group four union one two");
    lmp->input->one("group five subtract all half four");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->group->count(lmp->group->find("one")), 16);
    ASSERT_EQ(lmp->group->count(lmp->group->find("two")), 16);
    ASSERT_EQ(lmp->group->count(lmp->group->find("three")), 0);
    ASSERT_EQ(lmp->group->count(lmp->group->find("half")), 32);
    ASSERT_EQ(lmp->group->count(lmp->group->find("four")), 32);
    ASSERT_EQ(lmp->group->count(lmp->group->find("five")), 16);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("write_restart group.restart");
    lmp->input->one("clear");
    lmp->input->one("read_restart group.restart");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->group->count(lmp->group->find("one")), 16);
    ASSERT_EQ(lmp->group->count(lmp->group->find("two")), 16);
    ASSERT_EQ(lmp->group->count(lmp->group->find("three")), 0);
    ASSERT_EQ(lmp->group->count(lmp->group->find("half")), 32);
    ASSERT_EQ(lmp->group->count(lmp->group->find("four")), 32);
    ASSERT_EQ(lmp->group->count(lmp->group->find("five")), 16);
    ASSERT_EQ(lmp->group->count(lmp->group->find("half"), lmp->domain->find_region("left")), 8);

#if 0
    TEST_FAILURE(".*ERROR: Illegal group command.*",
                 lmp->input->one("group three region left xxx"););
    TEST_FAILURE(".*ERROR: Group region ID does not exist.*",
                 lmp->input->one("group four region top"););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("group one clear");
    lmp->input->one("group two clear");
    lmp->input->one("group three clear");
    lmp->input->one("group four clear");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->group->count(lmp->group->find("one")), 0);
    ASSERT_EQ(lmp->group->count(lmp->group->find("two")), 0);
    ASSERT_EQ(lmp->group->count(lmp->group->find("three")), 0);
    ASSERT_EQ(lmp->group->count(lmp->group->find("four")), 0);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("delete_atoms region box");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->group->count(lmp->group->find("all")), 0);
#endif
}
} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (have_openmpi && !LAMMPS_NS::Info::has_exceptions())
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
