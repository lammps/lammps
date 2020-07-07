/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom.h"
#include "fmt/format.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "output.h"
#include "update.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <mpi.h>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

using LAMMPS_NS::utils::split_words;

namespace LAMMPS_NS {
using ::testing::MatchesRegex;

#define GETIDX(i) lmp->atom->map(i)

#define TEST_FAILURE(...)                \
    if (Info::has_exceptions()) {        \
        ASSERT_ANY_THROW({__VA_ARGS__}); \
    } else {                             \
        ASSERT_DEATH({__VA_ARGS__}, ""); \
    }

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

class ResetMolIDsTest : public ::testing::Test {
protected:
    LAMMPS *lmp;

    void SetUp() override
    {
        const char *args[] = {"ResetMolIDsTest", "-log", "none", "-nocite", "-echo", "screen"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        Info *info = new Info(lmp);
        if (info->has_style("atom","full")) {
            lmp->input->one("variable input_dir index " STRINGIFY(TEST_INPUT_FOLDER));
            lmp->input->one("include ${input_dir}/in.fourmol");
        }
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void TearDown() override
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        delete lmp;
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }
};

TEST_F(ResetMolIDsTest, Plain)
{
    if (lmp->atom->natoms == 0) GTEST_SKIP();

    auto molid = lmp->atom->molecule;
    ASSERT_EQ(molid[GETIDX(1)], 1);
    ASSERT_EQ(molid[GETIDX(2)], 1);
    ASSERT_EQ(molid[GETIDX(3)], 1);
    ASSERT_EQ(molid[GETIDX(4)], 1);
    ASSERT_EQ(molid[GETIDX(5)], 1);
    ASSERT_EQ(molid[GETIDX(6)], 1);
    ASSERT_EQ(molid[GETIDX(7)], 1);
    ASSERT_EQ(molid[GETIDX(8)], 2);
    ASSERT_EQ(molid[GETIDX(9)], 2);
    ASSERT_EQ(molid[GETIDX(10)], 2);
    ASSERT_EQ(molid[GETIDX(11)], 2);
    ASSERT_EQ(molid[GETIDX(12)], 2);
    ASSERT_EQ(molid[GETIDX(13)], 2);
    ASSERT_EQ(molid[GETIDX(14)], 2);
    ASSERT_EQ(molid[GETIDX(15)], 2);
    ASSERT_EQ(molid[GETIDX(16)], 2);
    ASSERT_EQ(molid[GETIDX(17)], 2);
    ASSERT_EQ(molid[GETIDX(18)], 3);
    ASSERT_EQ(molid[GETIDX(19)], 3);
    ASSERT_EQ(molid[GETIDX(20)], 3);
    ASSERT_EQ(molid[GETIDX(21)], 4);
    ASSERT_EQ(molid[GETIDX(22)], 4);
    ASSERT_EQ(molid[GETIDX(23)], 4);
    ASSERT_EQ(molid[GETIDX(24)], 5);
    ASSERT_EQ(molid[GETIDX(25)], 5);
    ASSERT_EQ(molid[GETIDX(26)], 5);
    ASSERT_EQ(molid[GETIDX(27)], 6);
    ASSERT_EQ(molid[GETIDX(28)], 6);
    ASSERT_EQ(molid[GETIDX(29)], 6);

    // the original data file has two different molecule IDs
    // for two residues of the same molecule/fragment.
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("reset_mol_ids all");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_EQ(molid[GETIDX(1)], 1);
    ASSERT_EQ(molid[GETIDX(2)], 1);
    ASSERT_EQ(molid[GETIDX(3)], 1);
    ASSERT_EQ(molid[GETIDX(4)], 1);
    ASSERT_EQ(molid[GETIDX(5)], 1);
    ASSERT_EQ(molid[GETIDX(6)], 1);
    ASSERT_EQ(molid[GETIDX(7)], 1);
    ASSERT_EQ(molid[GETIDX(8)], 1);
    ASSERT_EQ(molid[GETIDX(9)], 1);
    ASSERT_EQ(molid[GETIDX(10)], 1);
    ASSERT_EQ(molid[GETIDX(11)], 1);
    ASSERT_EQ(molid[GETIDX(12)], 1);
    ASSERT_EQ(molid[GETIDX(13)], 1);
    ASSERT_EQ(molid[GETIDX(14)], 1);
    ASSERT_EQ(molid[GETIDX(15)], 1);
    ASSERT_EQ(molid[GETIDX(16)], 1);
    ASSERT_EQ(molid[GETIDX(17)], 1);
    ASSERT_EQ(molid[GETIDX(18)], 2);
    ASSERT_EQ(molid[GETIDX(19)], 2);
    ASSERT_EQ(molid[GETIDX(20)], 2);
    ASSERT_EQ(molid[GETIDX(21)], 3);
    ASSERT_EQ(molid[GETIDX(22)], 3);
    ASSERT_EQ(molid[GETIDX(23)], 3);
    ASSERT_EQ(molid[GETIDX(24)], 4);
    ASSERT_EQ(molid[GETIDX(25)], 4);
    ASSERT_EQ(molid[GETIDX(26)], 4);
    ASSERT_EQ(molid[GETIDX(27)], 5);
    ASSERT_EQ(molid[GETIDX(28)], 5);
    ASSERT_EQ(molid[GETIDX(29)], 5);
}

TEST_F(ResetMolIDsTest, DeletePlusAtomID)
{
    if (lmp->atom->natoms == 0) GTEST_SKIP();

    auto molid = lmp->atom->molecule;

    // delete two water molecules
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("group allwater molecule 3:6");
    lmp->input->one("group twowater molecule 4:6:2");
    lmp->input->one("delete_atoms group twowater compress no bond yes");
    lmp->input->one("reset_mol_ids all");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->atom->natoms, 23);
    ASSERT_EQ(lmp->atom->map_tag_max, 26);

    EXPECT_EQ(molid[GETIDX(1)], 1);
    EXPECT_EQ(molid[GETIDX(2)], 1);
    EXPECT_EQ(molid[GETIDX(3)], 1);
    EXPECT_EQ(molid[GETIDX(4)], 1);
    EXPECT_EQ(molid[GETIDX(5)], 1);
    EXPECT_EQ(molid[GETIDX(6)], 1);
    EXPECT_EQ(molid[GETIDX(7)], 1);
    EXPECT_EQ(molid[GETIDX(8)], 1);
    EXPECT_EQ(molid[GETIDX(9)], 1);
    EXPECT_EQ(molid[GETIDX(10)], 1);
    EXPECT_EQ(molid[GETIDX(11)], 1);
    EXPECT_EQ(molid[GETIDX(12)], 1);
    EXPECT_EQ(molid[GETIDX(13)], 1);
    EXPECT_EQ(molid[GETIDX(14)], 1);
    EXPECT_EQ(molid[GETIDX(15)], 1);
    EXPECT_EQ(molid[GETIDX(16)], 1);
    EXPECT_EQ(molid[GETIDX(17)], 1);
    EXPECT_EQ(molid[GETIDX(18)], 2);
    EXPECT_EQ(molid[GETIDX(19)], 2);
    EXPECT_EQ(molid[GETIDX(20)], 2);
    EXPECT_EQ(molid[GETIDX(24)], 3);
    EXPECT_EQ(molid[GETIDX(25)], 3);
    EXPECT_EQ(molid[GETIDX(26)], 3);

    // now also check and reset the atom ids

    EXPECT_GE(GETIDX(1), 0);
    EXPECT_GE(GETIDX(2), 0);
    EXPECT_GE(GETIDX(3), 0);
    EXPECT_GE(GETIDX(4), 0);
    EXPECT_GE(GETIDX(5), 0);
    EXPECT_GE(GETIDX(6), 0);
    EXPECT_GE(GETIDX(7), 0);
    EXPECT_GE(GETIDX(8), 0);
    EXPECT_GE(GETIDX(9), 0);
    EXPECT_GE(GETIDX(10), 0);
    EXPECT_GE(GETIDX(11), 0);
    EXPECT_GE(GETIDX(12), 0);
    EXPECT_GE(GETIDX(13), 0);
    EXPECT_GE(GETIDX(14), 0);
    EXPECT_GE(GETIDX(15), 0);
    EXPECT_GE(GETIDX(16), 0);
    EXPECT_GE(GETIDX(17), 0);
    EXPECT_GE(GETIDX(18), 0);
    EXPECT_GE(GETIDX(19), 0);
    EXPECT_GE(GETIDX(20), 0);
    EXPECT_EQ(GETIDX(21), -1);
    EXPECT_EQ(GETIDX(22), -1);
    EXPECT_EQ(GETIDX(23), -1);
    EXPECT_GE(GETIDX(24), 0);
    EXPECT_GE(GETIDX(25), 0);
    EXPECT_GE(GETIDX(26), 0);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("reset_ids");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_EQ(lmp->atom->map_tag_max, 23);
    for (int i = 1; i <= 23; ++i)
        EXPECT_GE(GETIDX(i), 0);
}

// TEST_FAILURE(lmp->input->one("reset_mol_ids all"););
// ASSERT_THAT(mesg, MatchesRegex(".*ERROR: Unknown command.*"));

} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

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
