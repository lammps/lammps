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
#include "output.h"
#include "update.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mpi.h>

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
        ASSERT_THAT(mesg, MatchesRegex(errmsg));                  \
    } else {                                                      \
        if (!have_openmpi) {                                      \
            ::testing::internal::CaptureStdout();                 \
            ASSERT_DEATH({__VA_ARGS__}, "");                      \
            auto mesg = ::testing::internal::GetCapturedStdout(); \
            ASSERT_THAT(mesg, MatchesRegex(errmsg));              \
        }                                                         \
    }

class SimpleCommandsTest : public ::testing::Test {
protected:
    LAMMPS *lmp;

    void SetUp() override
    {
        const char *args[] = {"SimpleCommandsTest", "-log", "none", "-echo", "screen", "-nocite"};
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

TEST_F(SimpleCommandsTest, UnknownCommand)
{
    TEST_FAILURE(".*ERROR: Unknown command.*", lmp->input->one("XXX one two"););
}

TEST_F(SimpleCommandsTest, Echo)
{
    ASSERT_EQ(lmp->input->echo_screen, 1);
    ASSERT_EQ(lmp->input->echo_log, 0);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("echo none");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->input->echo_screen, 0);
    ASSERT_EQ(lmp->input->echo_log, 0);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("echo both");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->input->echo_screen, 1);
    ASSERT_EQ(lmp->input->echo_log, 1);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("echo screen");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->input->echo_screen, 1);
    ASSERT_EQ(lmp->input->echo_log, 0);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("echo log");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->input->echo_screen, 0);
    ASSERT_EQ(lmp->input->echo_log, 1);

    TEST_FAILURE(".*ERROR: Illegal echo command.*", lmp->input->one("echo"););
    TEST_FAILURE(".*ERROR: Illegal echo command.*", lmp->input->one("echo xxx"););
}

TEST_F(SimpleCommandsTest, Log)
{
    ASSERT_EQ(lmp->logfile, nullptr);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("log simple_command_test.log");
    lmp->input->one("print 'test1'");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_NE(lmp->logfile, nullptr);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("log none");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->logfile, nullptr);

    std::string text;
    std::ifstream in;
    in.open("simple_command_test.log");
    in >> text;
    in.close();
    ASSERT_THAT(text, StrEq("test1"));

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("log simple_command_test.log append");
    lmp->input->one("print 'test2'");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_NE(lmp->logfile, nullptr);
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("log none");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->logfile, nullptr);

    in.open("simple_command_test.log");
    in >> text;
    ASSERT_THAT(text, StrEq("test1"));
    in >> text;
    ASSERT_THAT(text, StrEq("test2"));
    in.close();
    remove("simple_command_test.log");

    TEST_FAILURE(".*ERROR: Illegal log command.*", lmp->input->one("log"););
}

TEST_F(SimpleCommandsTest, Quit)
{
    ::testing::internal::CaptureStdout();
    lmp->input->one("echo none");
    ::testing::internal::GetCapturedStdout();
    TEST_FAILURE(".*ERROR: Expected integer .*", lmp->input->one("quit xxx"););

    // the following tests must be skipped with OpenMPI due to using threads
    if (have_openmpi) GTEST_SKIP();
    ASSERT_EXIT(lmp->input->one("quit"), ExitedWithCode(0), "");
    ASSERT_EXIT(lmp->input->one("quit 9"), ExitedWithCode(9), "");
}

TEST_F(SimpleCommandsTest, ResetTimestep)
{
    ASSERT_EQ(lmp->update->ntimestep, 0);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("reset_timestep 10");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->update->ntimestep, 10);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("reset_timestep 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->update->ntimestep, 0);

    TEST_FAILURE(".*ERROR: Timestep must be >= 0.*", lmp->input->one("reset_timestep -10"););
    TEST_FAILURE(".*ERROR: Illegal reset_timestep .*", lmp->input->one("reset_timestep"););
    TEST_FAILURE(".*ERROR: Illegal reset_timestep .*", lmp->input->one("reset_timestep 10 10"););
    TEST_FAILURE(".*ERROR: Expected integer .*", lmp->input->one("reset_timestep xxx"););
}

TEST_F(SimpleCommandsTest, Suffix)
{
    ASSERT_EQ(lmp->suffix_enable, 0);
    ASSERT_EQ(lmp->suffix, nullptr);
    ASSERT_EQ(lmp->suffix2, nullptr);

    TEST_FAILURE(".*ERROR: May only enable suffixes after defining one.*",
                 lmp->input->one("suffix on"););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("suffix one");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_THAT(lmp->suffix, StrEq("one"));

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("suffix hybrid two three");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_THAT(lmp->suffix, StrEq("two"));
    ASSERT_THAT(lmp->suffix2, StrEq("three"));

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("suffix four");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_THAT(lmp->suffix, StrEq("four"));
    ASSERT_EQ(lmp->suffix2, nullptr);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("suffix off");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->suffix_enable, 0);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("suffix on");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->suffix_enable, 1);

    TEST_FAILURE(".*ERROR: Illegal suffix command.*", lmp->input->one("suffix"););
    TEST_FAILURE(".*ERROR: Illegal suffix command.*", lmp->input->one("suffix hybrid"););
    TEST_FAILURE(".*ERROR: Illegal suffix command.*", lmp->input->one("suffix hybrid one"););
}

TEST_F(SimpleCommandsTest, Thermo)
{
    ASSERT_EQ(lmp->output->thermo_every, 0);
    ASSERT_EQ(lmp->output->var_thermo, nullptr);
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("thermo 2");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->output->thermo_every, 2);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("variable step equal logfreq(10,3,10)");
    lmp->input->one("thermo v_step");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_THAT(lmp->output->var_thermo, StrEq("step"));

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("thermo 10");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->output->thermo_every, 10);
    ASSERT_EQ(lmp->output->var_thermo, nullptr);

    TEST_FAILURE(".*ERROR: Illegal thermo command.*", lmp->input->one("thermo"););
    TEST_FAILURE(".*ERROR: Illegal thermo command.*", lmp->input->one("thermo -1"););
    TEST_FAILURE(".*ERROR: Expected integer.*", lmp->input->one("thermo xxx"););
}

TEST_F(SimpleCommandsTest, TimeStep)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("timestep 1");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->update->dt, 1.0);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("timestep 0.1");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->update->dt, 0.1);

    // zero timestep is legal and works (atoms don't move)
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("timestep 0.0");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->update->dt, 0.0);

    // negative timestep also creates a viable MD.
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("timestep -0.1");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->update->dt, -0.1);

    TEST_FAILURE(".*ERROR: Illegal timestep command.*", lmp->input->one("timestep"););
    TEST_FAILURE(".*ERROR: Expected floating point.*", lmp->input->one("timestep xxx"););
}

TEST_F(SimpleCommandsTest, Units)
{
    const char *names[] = {"lj", "real", "metal", "si", "cgs", "electron", "micro", "nano"};
    const double dt[]   = {0.005, 1.0, 0.001, 1.0e-8, 1.0e-8, 0.001, 2.0, 0.00045};
    std::size_t num     = sizeof(names) / sizeof(const char *);
    ASSERT_EQ(num, sizeof(dt) / sizeof(double));

    ASSERT_THAT(lmp->update->unit_style, StrEq("lj"));
    for (std::size_t i = 0; i < num; ++i) {
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp->input->one(fmt::format("units {}", names[i]));
        if (!verbose) ::testing::internal::GetCapturedStdout();
        ASSERT_THAT(lmp->update->unit_style, StrEq(names[i]));
        ASSERT_EQ(lmp->update->dt, dt[i]);
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_THAT(lmp->update->unit_style, StrEq("lj"));

    TEST_FAILURE(".*ERROR: Illegal units command.*", lmp->input->one("units unknown"););
}

TEST_F(SimpleCommandsTest, Shell)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("shell putenv TEST_VARIABLE=simpletest");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    char * test_var = getenv("TEST_VARIABLE");
    ASSERT_NE(test_var, nullptr);
    ASSERT_THAT(test_var, StrEq("simpletest"));

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("shell putenv TEST_VARIABLE=simpletest");
    lmp->input->one("shell putenv TEST_VARIABLE2=simpletest2 OTHER_VARIABLE=2");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    char * test_var2 = getenv("TEST_VARIABLE2");
    char * other_var = getenv("OTHER_VARIABLE");

    ASSERT_NE(test_var2, nullptr);
    ASSERT_THAT(test_var2, StrEq("simpletest2"));

    ASSERT_NE(other_var, nullptr);
    ASSERT_THAT(other_var, StrEq("2"));
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
