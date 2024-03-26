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

#include "citeme.h"
#include "comm.h"
#include "force.h"
#include "info.h"
#include "input.h"
#include "output.h"
#include "update.h"
#include "utils.h"
#include "variable.h"

#include "../testing/core.h"
#include "../testing/utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mpi.h>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {
using ::testing::ContainsRegex;
using ::testing::ExitedWithCode;
using ::testing::StrEq;

class SimpleCommandsTest : public LAMMPSTest {};

TEST_F(SimpleCommandsTest, UnknownCommand)
{
    TEST_FAILURE(".*ERROR: Unknown command.*", command("XXX one two"););
}

TEST_F(SimpleCommandsTest, Echo)
{
    ASSERT_EQ(lmp->input->echo_screen, 1);
    ASSERT_EQ(lmp->input->echo_log, 0);

    BEGIN_HIDE_OUTPUT();
    command("echo none");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->input->echo_screen, 0);
    ASSERT_EQ(lmp->input->echo_log, 0);

    BEGIN_HIDE_OUTPUT();
    command("echo both");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->input->echo_screen, 1);
    ASSERT_EQ(lmp->input->echo_log, 1);

    BEGIN_HIDE_OUTPUT();
    command("echo screen");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->input->echo_screen, 1);
    ASSERT_EQ(lmp->input->echo_log, 0);

    BEGIN_HIDE_OUTPUT();
    command("echo log");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->input->echo_screen, 0);
    ASSERT_EQ(lmp->input->echo_log, 1);

    TEST_FAILURE(".*ERROR: Illegal echo command.*", command("echo"););
    TEST_FAILURE(".*ERROR: Unknown echo keyword: xxx.*", command("echo xxx"););
}

TEST_F(SimpleCommandsTest, Log)
{
    ASSERT_EQ(lmp->logfile, nullptr);

    BEGIN_HIDE_OUTPUT();
    command("log simple_command_test.log");
    command("print 'test1'");
    END_HIDE_OUTPUT();
    ASSERT_NE(lmp->logfile, nullptr);

    BEGIN_HIDE_OUTPUT();
    command("log none");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->logfile, nullptr);

    std::string text;
    std::ifstream in;
    in.open("simple_command_test.log");
    in >> text;
    in.close();
    ASSERT_THAT(text, StrEq("test1"));

    BEGIN_HIDE_OUTPUT();
    command("log simple_command_test.log append");
    command("print 'test2'");
    END_HIDE_OUTPUT();
    ASSERT_NE(lmp->logfile, nullptr);
    BEGIN_HIDE_OUTPUT();
    command("log none");
    END_HIDE_OUTPUT();
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
    BEGIN_HIDE_OUTPUT();
    command("newton off");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->force->newton_pair, 0);
    ASSERT_EQ(lmp->force->newton_bond, 0);
    BEGIN_HIDE_OUTPUT();
    command("newton on off");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->force->newton_pair, 1);
    ASSERT_EQ(lmp->force->newton_bond, 0);
    BEGIN_HIDE_OUTPUT();
    command("newton off on");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->force->newton_pair, 0);
    ASSERT_EQ(lmp->force->newton_bond, 1);
    BEGIN_HIDE_OUTPUT();
    command("newton on");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->force->newton_pair, 1);
    ASSERT_EQ(lmp->force->newton_bond, 1);
}

TEST_F(SimpleCommandsTest, Partition)
{
    BEGIN_HIDE_OUTPUT();
    command("echo none");
    END_HIDE_OUTPUT();
    TEST_FAILURE(".*ERROR: Expected boolean parameter instead of 'xxx'.*",
                 command("partition xxx 1 echo none"););
    TEST_FAILURE(".*ERROR: Numeric index 2 is out of bounds.*",
                 command("partition yes 2 echo none"););

    BEGIN_CAPTURE_OUTPUT();
    command("partition yes 1 print 'test'");
    auto text = END_CAPTURE_OUTPUT();
    ASSERT_THAT(text, StrEq("test\n"));

    BEGIN_CAPTURE_OUTPUT();
    command("partition no 1 print 'test'");
    text = END_CAPTURE_OUTPUT();
    ASSERT_THAT(text, StrEq(""));
}

TEST_F(SimpleCommandsTest, Processors)
{
    // default setting is "*" for all dimensions
    ASSERT_EQ(lmp->comm->user_procgrid[0], 0);
    ASSERT_EQ(lmp->comm->user_procgrid[1], 0);
    ASSERT_EQ(lmp->comm->user_procgrid[2], 0);

    BEGIN_HIDE_OUTPUT();
    command("processors 1 1 1");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->comm->user_procgrid[0], 1);
    ASSERT_EQ(lmp->comm->user_procgrid[1], 1);
    ASSERT_EQ(lmp->comm->user_procgrid[2], 1);

    BEGIN_HIDE_OUTPUT();
    command("processors * 1 *");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->comm->user_procgrid[0], 0);
    ASSERT_EQ(lmp->comm->user_procgrid[1], 1);
    ASSERT_EQ(lmp->comm->user_procgrid[2], 0);

    BEGIN_HIDE_OUTPUT();
    command("processors 0 0 0");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->comm->user_procgrid[0], 0);
    ASSERT_EQ(lmp->comm->user_procgrid[1], 0);
    ASSERT_EQ(lmp->comm->user_procgrid[2], 0);

    TEST_FAILURE(".*ERROR: Illegal processors command .*", command("processors -1 0 0"););
    TEST_FAILURE(".*ERROR: Specified processors != physical processors.*",
                 command("processors 100 100 100"););
}

TEST_F(SimpleCommandsTest, Quit)
{
    BEGIN_HIDE_OUTPUT();
    command("echo none");
    END_HIDE_OUTPUT();
    TEST_FAILURE(".*ERROR: Expected integer .*", command("quit xxx"););

    // the following tests must be skipped with OpenMPI or MPICH 4.1 and later due to using threads
    if (platform::mpi_vendor() == "Open MPI") GTEST_SKIP();
#if defined(MPICH_NUMVERSION)
    if (MPICH_NUMVERSION >= 40100000) GTEST_SKIP();
#endif
    ASSERT_EXIT(command("quit"), ExitedWithCode(0), "");
    ASSERT_EXIT(command("quit 9"), ExitedWithCode(9), "");
}

TEST_F(SimpleCommandsTest, ResetTimestep)
{
    ASSERT_EQ(lmp->update->ntimestep, 0);
    ASSERT_EQ(lmp->update->atimestep, 0);
    ASSERT_DOUBLE_EQ(lmp->update->atime, 0.0);

    BEGIN_HIDE_OUTPUT();
    command("reset_timestep 10");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->update->ntimestep, 10);
    ASSERT_EQ(lmp->update->atimestep, 10);
    ASSERT_DOUBLE_EQ(lmp->update->atime, lmp->update->dt * 10);

    BEGIN_HIDE_OUTPUT();
    command("reset_timestep 0");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->update->ntimestep, 0);
    ASSERT_EQ(lmp->update->atimestep, 0);
    ASSERT_DOUBLE_EQ(lmp->update->atime, 0.0);

    BEGIN_HIDE_OUTPUT();
    command("reset_timestep 10 time 100.0");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->update->ntimestep, 10);
    ASSERT_EQ(lmp->update->atimestep, 10);
    ASSERT_DOUBLE_EQ(lmp->update->atime, 100.0);

    TEST_FAILURE(".*ERROR: Timestep must be >= 0.*", command("reset_timestep -10"););
    TEST_FAILURE(".*ERROR: Illegal reset_timestep .*", command("reset_timestep"););
    TEST_FAILURE(".*ERROR: Unknown reset_timestep option 10.*", command("reset_timestep 10 10"););
    TEST_FAILURE(".*ERROR: Illegal reset_timestep .*", command("reset_timestep 10 time"););
    TEST_FAILURE(".*ERROR: Expected floating .*", command("reset_timestep 10 time xxx"););
    TEST_FAILURE(".*ERROR: Expected integer .*", command("reset_timestep xxx"););
}

TEST_F(SimpleCommandsTest, Suffix)
{
    ASSERT_EQ(lmp->suffix_enable, 0);
    ASSERT_EQ(lmp->suffix, nullptr);
    ASSERT_EQ(lmp->suffix2, nullptr);

    TEST_FAILURE(".*ERROR: May only enable suffixes after defining one.*", command("suffix on"););
    TEST_FAILURE(".*ERROR: May only enable suffixes after defining one.*", command("suffix yes"););
    TEST_FAILURE(".*ERROR: May only enable suffixes after defining one.*", command("suffix true"););

    BEGIN_HIDE_OUTPUT();
    command("suffix one");
    command("suffix yes");
    END_HIDE_OUTPUT();
    ASSERT_THAT(lmp->suffix, StrEq("one"));
    ASSERT_EQ(lmp->suffix_enable, 1);
    ASSERT_THAT(utils::strip_style_suffix("one/four", lmp), StrEq("one/four"));
    ASSERT_THAT(utils::strip_style_suffix("four/one", lmp), StrEq("four"));

    BEGIN_HIDE_OUTPUT();
    command("suffix hybrid two three");
    END_HIDE_OUTPUT();
    ASSERT_THAT(lmp->suffix, StrEq("two"));
    ASSERT_THAT(lmp->suffix2, StrEq("three"));
    ASSERT_THAT(utils::strip_style_suffix("one/four", lmp), StrEq("one/four"));
    ASSERT_THAT(utils::strip_style_suffix("one/two", lmp), StrEq("one"));
    ASSERT_THAT(utils::strip_style_suffix("one/three", lmp), StrEq("one"));
    ASSERT_THAT(utils::strip_style_suffix("four/one", lmp), StrEq("four/one"));

    BEGIN_HIDE_OUTPUT();
    command("suffix four");
    END_HIDE_OUTPUT();
    ASSERT_THAT(lmp->suffix, StrEq("four"));
    ASSERT_EQ(lmp->suffix2, nullptr);

    BEGIN_HIDE_OUTPUT();
    command("suffix off");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->suffix_enable, 0);
    ASSERT_THAT(utils::strip_style_suffix("one/four", lmp), StrEq("one/four"));

    BEGIN_HIDE_OUTPUT();
    command("suffix yes");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->suffix_enable, 1);
    ASSERT_THAT(utils::strip_style_suffix("one/four", lmp), StrEq("one"));

    BEGIN_HIDE_OUTPUT();
    command("suffix no");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->suffix_enable, 0);

    BEGIN_HIDE_OUTPUT();
    command("suffix true");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->suffix_enable, 1);

    BEGIN_HIDE_OUTPUT();
    command("suffix false");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->suffix_enable, 0);

    BEGIN_HIDE_OUTPUT();
    command("suffix on");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->suffix_enable, 1);

    TEST_FAILURE(".*ERROR: Illegal suffix command.*", command("suffix"););
    TEST_FAILURE(".*ERROR: Illegal suffix command.*", command("suffix hybrid"););
    TEST_FAILURE(".*ERROR: Illegal suffix command.*", command("suffix hybrid one"););
}

TEST_F(SimpleCommandsTest, Thermo)
{
    ASSERT_EQ(lmp->output->thermo_every, 0);
    ASSERT_EQ(lmp->output->var_thermo, nullptr);
    BEGIN_HIDE_OUTPUT();
    command("thermo 2");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->output->thermo_every, 2);

    BEGIN_HIDE_OUTPUT();
    command("variable step equal logfreq(10,3,10)");
    command("thermo v_step");
    END_HIDE_OUTPUT();
    ASSERT_THAT(lmp->output->var_thermo, StrEq("step"));

    BEGIN_HIDE_OUTPUT();
    command("thermo 10");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->output->thermo_every, 10);
    ASSERT_EQ(lmp->output->var_thermo, nullptr);

    TEST_FAILURE(".*ERROR: Illegal thermo command.*", command("thermo"););
    TEST_FAILURE(".*ERROR: Illegal thermo output frequency.*", command("thermo -1"););
    TEST_FAILURE(".*ERROR: Expected integer.*", command("thermo xxx"););
}

TEST_F(SimpleCommandsTest, TimeStep)
{
    BEGIN_HIDE_OUTPUT();
    command("timestep 1");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->update->dt, 1.0);

    BEGIN_HIDE_OUTPUT();
    command("timestep 0.1");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->update->dt, 0.1);

    // zero timestep is legal and works (atoms don't move)
    BEGIN_HIDE_OUTPUT();
    command("timestep 0.0");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->update->dt, 0.0);

    // negative timestep also creates a viable MD.
    BEGIN_HIDE_OUTPUT();
    command("timestep -0.1");
    END_HIDE_OUTPUT();
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
        BEGIN_HIDE_OUTPUT();
        command(fmt::format("units {}", names[i]));
        END_HIDE_OUTPUT();
        ASSERT_THAT(lmp->update->unit_style, StrEq(names[i]));
        ASSERT_EQ(lmp->update->dt, dt[i]);
    }

    BEGIN_HIDE_OUTPUT();
    command("clear");
    END_HIDE_OUTPUT();
    ASSERT_THAT(lmp->update->unit_style, StrEq("lj"));

    TEST_FAILURE(".*ERROR: Illegal units command.*", command("units unknown"););
}

#if defined(LMP_PLUGIN)
TEST_F(SimpleCommandsTest, Plugin)
{
    const char *bindir = getenv("LAMMPS_PLUGIN_BIN_DIR");
    const char *config = getenv("CMAKE_CONFIG_TYPE");
    if (!bindir) GTEST_SKIP();
    std::string loadfmt = platform::path_join("plugin load ", bindir);
    if (config) loadfmt = platform::path_join(loadfmt, config);
    loadfmt = platform::path_join(loadfmt, "{}plugin.so");
    ::testing::internal::CaptureStdout();
    lmp->input->one(fmt::format(loadfmt, "hello"));
    auto text = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << text;
    ASSERT_THAT(text, ContainsRegex(".*Loading plugin: Hello world command.*"));

    ::testing::internal::CaptureStdout();
    lmp->input->one(fmt::format(loadfmt, "xxx"));
    text = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << text;
    ASSERT_THAT(text, ContainsRegex(".*Open of file .*xxx.* failed.*"));

    ::testing::internal::CaptureStdout();
    lmp->input->one(fmt::format(loadfmt, "nve2"));
    text = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << text;
    ASSERT_THAT(text, ContainsRegex(".*Loading plugin: NVE2 variant fix style.*"));
    ::testing::internal::CaptureStdout();
    lmp->input->one("plugin list");
    text = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << text;
    ASSERT_THAT(text, ContainsRegex(".*1: command style plugin hello"
                                    ".*2: fix style plugin nve2.*"));

    ::testing::internal::CaptureStdout();
    lmp->input->one(fmt::format(loadfmt, "hello"));
    text = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << text;
    ASSERT_THAT(text, ContainsRegex(".*Ignoring load of command style hello: "
                                    "must unload existing hello plugin.*"));

    ::testing::internal::CaptureStdout();
    lmp->input->one("plugin unload command hello");
    text = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << text;
    ASSERT_THAT(text, ContainsRegex(".*Unloading command style hello.*"));

    ::testing::internal::CaptureStdout();
    lmp->input->one("plugin unload pair nve2");
    text = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << text;
    ASSERT_THAT(text, ContainsRegex(".*Ignoring unload of pair style nve2: not from a plugin.*"));

    ::testing::internal::CaptureStdout();
    lmp->input->one("plugin unload fix nve2");
    text = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << text;
    ASSERT_THAT(text, ContainsRegex(".*Unloading fix style nve2.*"));

    ::testing::internal::CaptureStdout();
    lmp->input->one("plugin unload fix nve");
    text = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << text;
    ASSERT_THAT(text, ContainsRegex(".*Ignoring unload of fix style nve: not from a plugin.*"));

    ::testing::internal::CaptureStdout();
    lmp->input->one("plugin list");
    text = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << text;
    ASSERT_THAT(text, ContainsRegex(".*Currently loaded plugins.*"));
}
#endif

TEST_F(SimpleCommandsTest, Shell)
{
    BEGIN_HIDE_OUTPUT();
    command("shell putenv TEST_VARIABLE=simpletest");
    END_HIDE_OUTPUT();

    const char *test_var = getenv("TEST_VARIABLE");
    ASSERT_NE(test_var, nullptr);
    ASSERT_THAT(test_var, StrEq("simpletest"));

    BEGIN_HIDE_OUTPUT();
    command("shell putenv TEST_VARIABLE");
    command("shell putenv TEST_VARIABLE2=simpletest OTHER_VARIABLE=2");
    END_HIDE_OUTPUT();

    test_var = getenv("TEST_VARIABLE2");
    ASSERT_NE(test_var, nullptr);
    ASSERT_THAT(test_var, StrEq("simpletest"));

    test_var = getenv("OTHER_VARIABLE");
    ASSERT_NE(test_var, nullptr);
    ASSERT_THAT(test_var, StrEq("2"));

    test_var = getenv("TEST_VARIABLE");
    ASSERT_NE(test_var, nullptr);
#if defined(_WIN32)
    // we cannot create empty environment variables on Windows so platform::putenv() sets their
    // value to "1"
    ASSERT_THAT(test_var, StrEq("1"));
#else
    ASSERT_THAT(test_var, StrEq(""));
#endif
}

TEST_F(SimpleCommandsTest, CiteMe)
{
    ASSERT_EQ(lmp->citeme, nullptr);

    lmp->citeme = new LAMMPS_NS::CiteMe(lmp, CiteMe::TERSE, CiteMe::TERSE, nullptr);

    BEGIN_CAPTURE_OUTPUT();
    lmp->citeme->add("test citation one:\n 1\n");
    lmp->citeme->add("test citation two:\n 2\n");
    lmp->citeme->add("test citation one:\n 1\n");
    lmp->citeme->flush();
    std::string text = END_CAPTURE_OUTPUT();

    // find the two unique citations, but not the third
    ASSERT_THAT(text, ContainsRegex("test citation one.\n.*test citation two.*"));
    ASSERT_THAT(text, Not(ContainsRegex(
                          "test citation one.\n.*test citation two.*\n.*test citation one.*")));

    BEGIN_CAPTURE_OUTPUT();
    lmp->citeme->add("test citation one:\n 0\n");
    lmp->citeme->add("test citation two:\n 2\n");
    lmp->citeme->add("test citation three:\n 3\n");
    lmp->citeme->flush();

    text = END_CAPTURE_OUTPUT();

    // find the forth (only differs in long citation) and sixth added citation
    ASSERT_THAT(text, ContainsRegex("test citation one.*\n.*test citation three.*"));
    ASSERT_THAT(text, Not(ContainsRegex("test_citation two.*\n")));

    BEGIN_CAPTURE_OUTPUT();
    lmp->citeme->add("test citation one:\n 1\n");
    lmp->citeme->add("test citation two:\n 2\n");
    lmp->citeme->add("test citation one:\n 0\n");
    lmp->citeme->add("test citation two:\n 2\n");
    lmp->citeme->add("test citation three:\n 3\n");
    lmp->citeme->flush();

    text = END_CAPTURE_OUTPUT();

    // no new citation. no CITE-CITE-CITE- lines
    ASSERT_THAT(text, Not(ContainsRegex(".*CITE-CITE-CITE-CITE.*")));
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
