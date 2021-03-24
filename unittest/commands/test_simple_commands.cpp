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

#include "citeme.h"
#include "force.h"
#include "info.h"
#include "input.h"
#include "output.h"
#include "update.h"
#include "utils.h"

#include "fmt/format.h"
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

    void command(const std::string &cmd) { lmp->input->one(cmd); }
};

TEST_F(SimpleCommandsTest, UnknownCommand)
{
    TEST_FAILURE(".*ERROR: Unknown command.*", command("XXX one two"););
}

TEST_F(SimpleCommandsTest, Echo)
{
    ASSERT_EQ(lmp->input->echo_screen, 1);
    ASSERT_EQ(lmp->input->echo_log, 0);

    if (!verbose) ::testing::internal::CaptureStdout();
    command("echo none");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->input->echo_screen, 0);
    ASSERT_EQ(lmp->input->echo_log, 0);

    if (!verbose) ::testing::internal::CaptureStdout();
    command("echo both");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->input->echo_screen, 1);
    ASSERT_EQ(lmp->input->echo_log, 1);

    if (!verbose) ::testing::internal::CaptureStdout();
    command("echo screen");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->input->echo_screen, 1);
    ASSERT_EQ(lmp->input->echo_log, 0);

    if (!verbose) ::testing::internal::CaptureStdout();
    command("echo log");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->input->echo_screen, 0);
    ASSERT_EQ(lmp->input->echo_log, 1);

    TEST_FAILURE(".*ERROR: Illegal echo command.*", command("echo"););
    TEST_FAILURE(".*ERROR: Illegal echo command.*", command("echo xxx"););
}

TEST_F(SimpleCommandsTest, Log)
{
    ASSERT_EQ(lmp->logfile, nullptr);

    if (!verbose) ::testing::internal::CaptureStdout();
    command("log simple_command_test.log");
    command("print 'test1'");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_NE(lmp->logfile, nullptr);

    if (!verbose) ::testing::internal::CaptureStdout();
    command("log none");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->logfile, nullptr);

    std::string text;
    std::ifstream in;
    in.open("simple_command_test.log");
    in >> text;
    in.close();
    ASSERT_THAT(text, StrEq("test1"));

    if (!verbose) ::testing::internal::CaptureStdout();
    command("log simple_command_test.log append");
    command("print 'test2'");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_NE(lmp->logfile, nullptr);
    if (!verbose) ::testing::internal::CaptureStdout();
    command("log none");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->logfile, nullptr);

    in.open("simple_command_test.log");
    in >> text;
    ASSERT_THAT(text, StrEq("test1"));
    in >> text;
    ASSERT_THAT(text, StrEq("test2"));
    in.close();
    remove("simple_command_test.log");

    TEST_FAILURE(".*ERROR: Illegal log command.*", command("log"););
}

TEST_F(SimpleCommandsTest, Newton)
{
    // default setting is "on" for both
    ASSERT_EQ(lmp->force->newton_pair, 1);
    ASSERT_EQ(lmp->force->newton_bond, 1);
    if (!verbose) ::testing::internal::CaptureStdout();
    command("newton off");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->force->newton_pair, 0);
    ASSERT_EQ(lmp->force->newton_bond, 0);
    if (!verbose) ::testing::internal::CaptureStdout();
    command("newton on off");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->force->newton_pair, 1);
    ASSERT_EQ(lmp->force->newton_bond, 0);
    if (!verbose) ::testing::internal::CaptureStdout();
    command("newton off on");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->force->newton_pair, 0);
    ASSERT_EQ(lmp->force->newton_bond, 1);
    if (!verbose) ::testing::internal::CaptureStdout();
    command("newton on");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->force->newton_pair, 1);
    ASSERT_EQ(lmp->force->newton_bond, 1);
}

TEST_F(SimpleCommandsTest, Partition)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("echo none");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    TEST_FAILURE(".*ERROR: Illegal partition command .*", command("partition xxx 1 echo none"););
    TEST_FAILURE(".*ERROR: Numeric index 2 is out of bounds.*",
                 command("partition yes 2 echo none"););

    ::testing::internal::CaptureStdout();
    command("partition yes 1 print 'test'");
    auto text = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << text;
    ASSERT_THAT(text, StrEq("test\n"));

    ::testing::internal::CaptureStdout();
    command("partition no 1 print 'test'");
    text = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << text;
    ASSERT_THAT(text, StrEq(""));
}

TEST_F(SimpleCommandsTest, Quit)
{
    ::testing::internal::CaptureStdout();
    command("echo none");
    ::testing::internal::GetCapturedStdout();
    TEST_FAILURE(".*ERROR: Expected integer .*", command("quit xxx"););

    // the following tests must be skipped with OpenMPI due to using threads
    if (have_openmpi) GTEST_SKIP();
    ASSERT_EXIT(command("quit"), ExitedWithCode(0), "");
    ASSERT_EXIT(command("quit 9"), ExitedWithCode(9), "");
}

TEST_F(SimpleCommandsTest, ResetTimestep)
{
    ASSERT_EQ(lmp->update->ntimestep, 0);

    if (!verbose) ::testing::internal::CaptureStdout();
    command("reset_timestep 10");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->update->ntimestep, 10);

    if (!verbose) ::testing::internal::CaptureStdout();
    command("reset_timestep 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->update->ntimestep, 0);

    TEST_FAILURE(".*ERROR: Timestep must be >= 0.*", command("reset_timestep -10"););
    TEST_FAILURE(".*ERROR: Illegal reset_timestep .*", command("reset_timestep"););
    TEST_FAILURE(".*ERROR: Illegal reset_timestep .*", command("reset_timestep 10 10"););
    TEST_FAILURE(".*ERROR: Expected integer .*", command("reset_timestep xxx"););
}

TEST_F(SimpleCommandsTest, Suffix)
{
    ASSERT_EQ(lmp->suffix_enable, 0);
    ASSERT_EQ(lmp->suffix, nullptr);
    ASSERT_EQ(lmp->suffix2, nullptr);

    TEST_FAILURE(".*ERROR: May only enable suffixes after defining one.*", command("suffix on"););

    if (!verbose) ::testing::internal::CaptureStdout();
    command("suffix one");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_THAT(lmp->suffix, StrEq("one"));

    if (!verbose) ::testing::internal::CaptureStdout();
    command("suffix hybrid two three");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_THAT(lmp->suffix, StrEq("two"));
    ASSERT_THAT(lmp->suffix2, StrEq("three"));

    if (!verbose) ::testing::internal::CaptureStdout();
    command("suffix four");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_THAT(lmp->suffix, StrEq("four"));
    ASSERT_EQ(lmp->suffix2, nullptr);

    if (!verbose) ::testing::internal::CaptureStdout();
    command("suffix off");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->suffix_enable, 0);

    if (!verbose) ::testing::internal::CaptureStdout();
    command("suffix on");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->suffix_enable, 1);

    TEST_FAILURE(".*ERROR: Illegal suffix command.*", command("suffix"););
    TEST_FAILURE(".*ERROR: Illegal suffix command.*", command("suffix hybrid"););
    TEST_FAILURE(".*ERROR: Illegal suffix command.*", command("suffix hybrid one"););
}

TEST_F(SimpleCommandsTest, Thermo)
{
    ASSERT_EQ(lmp->output->thermo_every, 0);
    ASSERT_EQ(lmp->output->var_thermo, nullptr);
    if (!verbose) ::testing::internal::CaptureStdout();
    command("thermo 2");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->output->thermo_every, 2);

    if (!verbose) ::testing::internal::CaptureStdout();
    command("variable step equal logfreq(10,3,10)");
    command("thermo v_step");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_THAT(lmp->output->var_thermo, StrEq("step"));

    if (!verbose) ::testing::internal::CaptureStdout();
    command("thermo 10");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->output->thermo_every, 10);
    ASSERT_EQ(lmp->output->var_thermo, nullptr);

    TEST_FAILURE(".*ERROR: Illegal thermo command.*", command("thermo"););
    TEST_FAILURE(".*ERROR: Illegal thermo command.*", command("thermo -1"););
    TEST_FAILURE(".*ERROR: Expected integer.*", command("thermo xxx"););
}

TEST_F(SimpleCommandsTest, TimeStep)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("timestep 1");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->update->dt, 1.0);

    if (!verbose) ::testing::internal::CaptureStdout();
    command("timestep 0.1");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->update->dt, 0.1);

    // zero timestep is legal and works (atoms don't move)
    if (!verbose) ::testing::internal::CaptureStdout();
    command("timestep 0.0");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->update->dt, 0.0);

    // negative timestep also creates a viable MD.
    if (!verbose) ::testing::internal::CaptureStdout();
    command("timestep -0.1");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->update->dt, -0.1);

    TEST_FAILURE(".*ERROR: Illegal timestep command.*", command("timestep"););
    TEST_FAILURE(".*ERROR: Expected floating point.*", command("timestep xxx"););
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
        command(fmt::format("units {}", names[i]));
        if (!verbose) ::testing::internal::GetCapturedStdout();
        ASSERT_THAT(lmp->update->unit_style, StrEq(names[i]));
        ASSERT_EQ(lmp->update->dt, dt[i]);
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    command("clear");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_THAT(lmp->update->unit_style, StrEq("lj"));

    TEST_FAILURE(".*ERROR: Illegal units command.*", command("units unknown"););
}

TEST_F(SimpleCommandsTest, Shell)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("shell putenv TEST_VARIABLE=simpletest");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    char *test_var = getenv("TEST_VARIABLE");
    ASSERT_NE(test_var, nullptr);
    ASSERT_THAT(test_var, StrEq("simpletest"));

    if (!verbose) ::testing::internal::CaptureStdout();
    command("shell putenv TEST_VARIABLE=simpletest");
    command("shell putenv TEST_VARIABLE2=simpletest2 OTHER_VARIABLE=2");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    char *test_var2 = getenv("TEST_VARIABLE2");
    char *other_var = getenv("OTHER_VARIABLE");

    ASSERT_NE(test_var2, nullptr);
    ASSERT_THAT(test_var2, StrEq("simpletest2"));

    ASSERT_NE(other_var, nullptr);
    ASSERT_THAT(other_var, StrEq("2"));
}

TEST_F(SimpleCommandsTest, CiteMe)
{
    ASSERT_EQ(lmp->citeme, nullptr);

    lmp->citeme = new LAMMPS_NS::CiteMe(lmp, CiteMe::TERSE, CiteMe::TERSE, nullptr);

    ::testing::internal::CaptureStdout();
    lmp->citeme->add("test citation one:\n 1\n");
    lmp->citeme->add("test citation two:\n 2\n");
    lmp->citeme->add("test citation one:\n 1\n");
    lmp->citeme->flush();
    std::string text = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << text;

    // find the two unique citations, but not the third
    ASSERT_THAT(text, MatchesRegex(".*one.*two.*"));
    ASSERT_THAT(text, Not(MatchesRegex(".*one.*two.*one.*")));

    ::testing::internal::CaptureStdout();
    lmp->citeme->add("test citation one:\n 0\n");
    lmp->citeme->add("test citation two:\n 2\n");
    lmp->citeme->add("test citation three:\n 3\n");
    lmp->citeme->flush();

    text = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << text;

    // find the forth (only differs in long citation) and sixth added citation
    ASSERT_THAT(text, MatchesRegex(".*one.*three.*"));
    ASSERT_THAT(text, Not(MatchesRegex(".*two.*")));

    ::testing::internal::CaptureStdout();
    lmp->citeme->add("test citation one:\n 1\n");
    lmp->citeme->add("test citation two:\n 2\n");
    lmp->citeme->add("test citation one:\n 0\n");
    lmp->citeme->add("test citation two:\n 2\n");
    lmp->citeme->add("test citation three:\n 3\n");
    lmp->citeme->flush();

    text = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << text;

    // no new citation. no CITE-CITE-CITE- lines
    ASSERT_THAT(text, Not(MatchesRegex(".*CITE-CITE-CITE-CITE.*")));
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
