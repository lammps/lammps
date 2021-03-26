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

#include "../testing/core.h"
#include "../testing/systems/melt.h"
#include "../testing/utils.h"
#include "fmt/format.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using ::testing::Eq;

bool verbose = false;

class DumpCfgTest : public MeltTest {
    std::string dump_style = "cfg";

public:
    void generate_dump(std::string dump_file, std::string fields, std::string dump_modify_options,
                       int ntimesteps)
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        command(fmt::format("dump id all {} 1 {} {}", dump_style, dump_file, fields));

        if (!dump_modify_options.empty()) {
            command(fmt::format("dump_modify id {}", dump_modify_options));
        }

        command(fmt::format("run {}", ntimesteps));
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }
};

TEST_F(DumpCfgTest, invalid_options)
{
    TEST_FAILURE(".*Dump cfg arguments must start with 'mass type xs ys zs'.*",
                 command("dump id all cfg 1 dump.cfg id type proc procp1 mass x y z"););
}

TEST_F(DumpCfgTest, require_multifile)
{
    auto dump_file = "dump.melt.cfg_run.cfg";
    auto fields =
        "mass type xs ys zs id proc procp1 x y z ix iy iz xu yu zu xsu ysu zsu vx vy vz fx fy fz";

    if (!verbose) ::testing::internal::CaptureStdout();
    command(fmt::format("dump id all cfg 1 {} {}", dump_file, fields));
    if (!verbose) ::testing::internal::GetCapturedStdout();

    TEST_FAILURE(".*Dump cfg requires one snapshot per file.*", command("run 0"););
}

TEST_F(DumpCfgTest, run0)
{
    auto dump_file = "dump_cfg_run*.melt.cfg";
    auto fields    = "mass type xs ys zs id proc procp1 x y z ix iy iz vx vy vz fx fy fz";

    generate_dump(dump_file, fields, "", 0);

    ASSERT_FILE_EXISTS("dump_cfg_run0.melt.cfg");
    auto lines = read_lines("dump_cfg_run0.melt.cfg");
    ASSERT_EQ(lines.size(), 124);
    ASSERT_THAT(lines[0], Eq("Number of particles = 32"));
    delete_file("dump_cfg_run0.melt.cfg");
}

TEST_F(DumpCfgTest, unwrap_run0)
{
    auto dump_file = "dump_cfg_unwrap_run*.melt.cfg";
    auto fields    = "mass type xsu ysu zsu id proc procp1 x y z ix iy iz vx vy vz fx fy fz";

    generate_dump(dump_file, fields, "", 0);

    ASSERT_FILE_EXISTS("dump_cfg_unwrap_run0.melt.cfg");
    auto lines = read_lines("dump_cfg_unwrap_run0.melt.cfg");
    ASSERT_EQ(lines.size(), 124);
    ASSERT_THAT(lines[0], Eq("Number of particles = 32"));
    delete_file("dump_cfg_unwrap_run0.melt.cfg");
}

TEST_F(DumpCfgTest, no_buffer_run0)
{
    auto dump_file = "dump_cfg_no_buffer_run*.melt.cfg";
    auto fields    = "mass type xsu ysu zsu id proc procp1 x y z ix iy iz vx vy vz fx fy fz";

    generate_dump(dump_file, fields, "buffer no", 0);

    ASSERT_FILE_EXISTS("dump_cfg_no_buffer_run0.melt.cfg");
    auto lines = read_lines("dump_cfg_no_buffer_run0.melt.cfg");
    ASSERT_EQ(lines.size(), 124);
    ASSERT_THAT(lines[0], Eq("Number of particles = 32"));
    delete_file("dump_cfg_no_buffer_run0.melt.cfg");
}

TEST_F(DumpCfgTest, no_unwrap_no_buffer_run0)
{
    auto dump_file = "dump_cfg_no_unwrap_no_buffer_run*.melt.cfg";
    auto fields    = "mass type xs ys zs id proc procp1 x y z ix iy iz vx vy vz fx fy fz";

    generate_dump(dump_file, fields, "buffer no", 0);

    ASSERT_FILE_EXISTS("dump_cfg_no_unwrap_no_buffer_run0.melt.cfg");
    auto lines = read_lines("dump_cfg_no_unwrap_no_buffer_run0.melt.cfg");
    ASSERT_EQ(lines.size(), 124);
    ASSERT_THAT(lines[0], Eq("Number of particles = 32"));
    delete_file("dump_cfg_no_unwrap_no_buffer_run0.melt.cfg");
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

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
