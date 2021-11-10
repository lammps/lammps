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

#include "../testing/core.h"
#include "../testing/systems/melt.h"
#include "../testing/utils.h"
#include "fmt/format.h"
#include "output.h"
#include "thermo.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <string>

using ::testing::Eq;

char *BINARY2TXT_BINARY = nullptr;
bool verbose            = false;

class DumpLocalTest : public MeltTest {
    std::string dump_style = "local";

public:
    void enable_triclinic()
    {
        BEGIN_HIDE_OUTPUT();
        command("change_box all triclinic");
        END_HIDE_OUTPUT();
    }

    void generate_dump(std::string dump_file, std::string dump_options,
                       std::string dump_modify_options, int ntimesteps)
    {
        BEGIN_HIDE_OUTPUT();
        command(fmt::format("dump id all {} 1 {} {}", dump_style, dump_file, dump_options));

        if (!dump_modify_options.empty()) {
            command(fmt::format("dump_modify id {}", dump_modify_options));
        }

        command(fmt::format("run {} post no", ntimesteps));
        END_HIDE_OUTPUT();
    }

    void continue_dump(int ntimesteps)
    {
        BEGIN_HIDE_OUTPUT();
        command(fmt::format("run {} pre no post no", ntimesteps));
        END_HIDE_OUTPUT();
    }

    void SetUp() override
    {
        MeltTest::SetUp();

        BEGIN_HIDE_OUTPUT();
        command("compute comp all pair/local dist eng");
        END_HIDE_OUTPUT();
    }
};

TEST_F(DumpLocalTest, run0)
{
    auto dump_file = "dump_local_run0.melt";
    generate_dump(dump_file, "index c_comp[1]", "", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 873);

    ASSERT_THAT(lines[0], Eq("ITEM: TIMESTEP"));
    ASSERT_EQ(std::stoi(lines[1]), 0);

    ASSERT_THAT(lines[2], Eq("ITEM: NUMBER OF ENTRIES"));
    ASSERT_EQ(std::stoi(lines[3]), 864);

    ASSERT_THAT(lines[4], Eq("ITEM: BOX BOUNDS pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[5]).size(), 2);
    ASSERT_EQ(utils::split_words(lines[6]).size(), 2);
    ASSERT_EQ(utils::split_words(lines[7]).size(), 2);
    ASSERT_THAT(lines[8], Eq("ITEM: ENTRIES index c_comp[1] "));
    ASSERT_EQ(utils::split_words(lines[9]).size(), 2);
    ASSERT_THAT(lines[9], Eq("1 1.18765 "));
    delete_file(dump_file);
}

TEST_F(DumpLocalTest, label_run0)
{
    auto dump_file = "dump_local_label_run0.melt";
    generate_dump(dump_file, "index c_comp[1]", "label ELEMENTS", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_THAT(lines[2], Eq("ITEM: NUMBER OF ELEMENTS"));
    ASSERT_THAT(lines[8], Eq("ITEM: ELEMENTS index c_comp[1] "));
    delete_file(dump_file);
}

TEST_F(DumpLocalTest, format_line_run0)
{
    auto dump_file = "dump_local_format_line_run0.melt";
    generate_dump(dump_file, "index c_comp[1]", "format line \"%d %20.8g\"", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);

    ASSERT_EQ(utils::split_words(lines[9]).size(), 2);
    ASSERT_THAT(lines[9], Eq("1            1.1876539 "));
    delete_file(dump_file);
}

TEST_F(DumpLocalTest, format_int_run0)
{
    auto dump_file = "dump_local_format_int_run0.melt";
    generate_dump(dump_file, "index c_comp[1]", "format int \"%20d\"", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);

    ASSERT_EQ(utils::split_words(lines[9]).size(), 2);
    ASSERT_THAT(lines[9], Eq("                   1 1.18765 "));
    delete_file(dump_file);
}

TEST_F(DumpLocalTest, format_float_run0)
{
    auto dump_file = "dump_local_format_float_run0.melt";
    generate_dump(dump_file, "index c_comp[1]", "format float \"%20.5g\"", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);

    ASSERT_EQ(utils::split_words(lines[9]).size(), 2);
    ASSERT_THAT(lines[9], Eq("1               1.1877 "));
    delete_file(dump_file);
}

TEST_F(DumpLocalTest, format_column_run0)
{
    auto dump_file = "dump_local_format_column_run0.melt";
    generate_dump(dump_file, "index c_comp[1]", "format 1 \"%20d\"", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);

    ASSERT_EQ(utils::split_words(lines[9]).size(), 2);
    ASSERT_THAT(lines[9], Eq("                   1 1.18765 "));
    delete_file(dump_file);
}

TEST_F(DumpLocalTest, no_buffer_run0)
{
    auto dump_file = "dump_local_format_line_run0.melt";
    generate_dump(dump_file, "index c_comp[1]", "buffer no", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 873);

    ASSERT_THAT(lines[0], Eq("ITEM: TIMESTEP"));
    ASSERT_EQ(std::stoi(lines[1]), 0);

    ASSERT_THAT(lines[2], Eq("ITEM: NUMBER OF ENTRIES"));
    ASSERT_EQ(std::stoi(lines[3]), 864);

    ASSERT_THAT(lines[4], Eq("ITEM: BOX BOUNDS pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[5]).size(), 2);
    ASSERT_EQ(utils::split_words(lines[6]).size(), 2);
    ASSERT_EQ(utils::split_words(lines[7]).size(), 2);
    ASSERT_THAT(lines[8], Eq("ITEM: ENTRIES index c_comp[1] "));
    ASSERT_EQ(utils::split_words(lines[9]).size(), 2);
    ASSERT_THAT(lines[9], Eq("1 1.18765 "));
    delete_file(dump_file);
}

TEST_F(DumpLocalTest, with_units_run0)
{
    auto dump_file = "dump_with_units_run0.melt";
    generate_dump(dump_file, "index c_comp[1]", "units yes", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 875);

    ASSERT_THAT(lines[0], Eq("ITEM: UNITS"));
    ASSERT_THAT(lines[1], Eq("lj"));

    ASSERT_THAT(lines[2], Eq("ITEM: TIMESTEP"));
    ASSERT_EQ(std::stoi(lines[3]), 0);

    ASSERT_THAT(lines[4], Eq("ITEM: NUMBER OF ENTRIES"));
    ASSERT_EQ(std::stoi(lines[5]), 864);
}

TEST_F(DumpLocalTest, with_time_run0)
{
    auto dump_file = "dump_with_time_run0.melt";
    generate_dump(dump_file, "index c_comp[1]", "time yes", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 875);

    ASSERT_THAT(lines[0], Eq("ITEM: TIME"));
    ASSERT_THAT(std::stof(lines[1]), 0.0);

    ASSERT_THAT(lines[2], Eq("ITEM: TIMESTEP"));
    ASSERT_EQ(std::stoi(lines[3]), 0);

    ASSERT_THAT(lines[4], Eq("ITEM: NUMBER OF ENTRIES"));
    ASSERT_EQ(std::stoi(lines[5]), 864);
}

TEST_F(DumpLocalTest, triclinic_run0)
{
    auto dump_file = "dump_local_triclinic_run0.melt";
    enable_triclinic();
    generate_dump(dump_file, "index c_comp[1]", "", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);

    ASSERT_THAT(lines[4], Eq("ITEM: BOX BOUNDS xy xz yz pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[5]).size(), 3);
    ASSERT_EQ(utils::split_words(lines[6]).size(), 3);
    ASSERT_EQ(utils::split_words(lines[7]).size(), 3);
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
