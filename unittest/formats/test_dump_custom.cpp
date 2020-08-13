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

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "fmt/format.h"
#include "utils.h"
#include "../testing/core.h"
#include "../testing/systems/melt.h"
#include "../testing/utils.h"

char * BINARY2TXT_BINARY = nullptr;

class DumpCustomTest : public MeltTest {
};

TEST_F(DumpCustomTest, run0)
{
    auto dump_file = "dump_custom_run0.melt";

    if (!verbose) ::testing::internal::CaptureStdout();
    command(fmt::format("dump id all custom 1 {} id type x y vx fx", dump_file));
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();


    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 41);
    ASSERT_STREQ(lines[4].c_str(), "ITEM: BOX BOUNDS pp pp pp");
    ASSERT_EQ(utils::split_words(lines[5]).size(), 2);
    ASSERT_STREQ(lines[8].c_str(), "ITEM: ATOMS id type x y vx fx");
    ASSERT_EQ(utils::split_words(lines[9]).size(), 6);
    delete_file(dump_file);
}

TEST_F(DumpCustomTest, triclinic_run0)
{
    auto dump_file = "dump_custom_tri_run0.melt";
    if (!verbose) ::testing::internal::CaptureStdout();

    command("change_box all triclinic");
    command(fmt::format("dump id all custom 1 {} id type x y vx fx", dump_file));
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS(dump_file);

    auto lines = read_lines(dump_file);
    ASSERT_STREQ(lines[4].c_str(), "ITEM: BOX BOUNDS xy xz yz pp pp pp");
    ASSERT_EQ(utils::split_words(lines[5]).size(), 3);

    ASSERT_EQ(lines.size(), 41);
    delete_file(dump_file);
}


int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = utils::split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }

    BINARY2TXT_BINARY = getenv("BINARY2TXT_BINARY");

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
