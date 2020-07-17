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

#include "fmt/format.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "modify.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdlib>
#include <mpi.h>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

using LAMMPS_NS::utils::split_words;

namespace LAMMPS_NS {
using ::testing::ExitedWithCode;
using ::testing::MatchesRegex;
using ::testing::StrEq;

#define TEST_FAILURE(...)                \
    if (Info::has_exceptions()) {        \
        ASSERT_ANY_THROW({__VA_ARGS__}); \
    } else {                             \
        ASSERT_DEATH({__VA_ARGS__}, ""); \
    }

class KimCommandsTest : public ::testing::Test {
protected:
    LAMMPS *lmp;

    void SetUp() override
    {
        const char *args[] = {"KimCommandsTest", "-log", "none", "-echo", "screen", "-nocite"};
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
    }
};

TEST_F(KimCommandsTest, kim_init)
{
    if (!LAMMPS::is_installed_pkg("KIM")) GTEST_SKIP();

    ::testing::internal::CaptureStdout();
    lmp->input->one("kim_init LennardJones_Ar real");
    ::testing::internal::GetCapturedStdout();

    int ifix = lmp->modify->find_fix("KIM_MODEL_STORE");
    ASSERT_GE(ifix, 0);

    ::testing::internal::CaptureStdout();
    TEST_FAILURE(lmp->input->one("kim_init"););
    auto mesg = ::testing::internal::GetCapturedStdout();
    ASSERT_THAT(mesg, MatchesRegex(".*ERROR: Illegal kim_init command.*"));

    ::testing::internal::CaptureStdout();
    TEST_FAILURE(lmp->input->one("kim_init Unknown_Model real"););
    mesg = ::testing::internal::GetCapturedStdout();
    ASSERT_THAT(mesg, MatchesRegex(".*ERROR: KIM Model name not found.*"));
}

TEST_F(KimCommandsTest, kim_interactions_ar)
{
    if (!LAMMPS::is_installed_pkg("KIM")) GTEST_SKIP();

    ::testing::internal::CaptureStdout();
    lmp->input->one("kim_init LennardJones_Ar real");
    lmp->input->one("lattice fcc 4.4300");
    lmp->input->one("region box block 0 10 0 10 0 10");
    lmp->input->one("create_box 1 box");
    lmp->input->one("create_atoms 1 box");
    lmp->input->one("kim_interactions Ar");
    lmp->input->one("mass 1 39.95");
    ::testing::internal::GetCapturedStdout();

    int ifix = lmp->modify->find_fix("KIM_MODEL_STORE");
    ASSERT_GE(ifix, 0);
}

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
