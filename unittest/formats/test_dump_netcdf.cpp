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
#include "library.h"
#include "output.h"
#include "thermo.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <string>

using ::testing::Eq;

char *NCDUMP_BINARY = nullptr;
bool verbose        = false;

class DumpNetCDFTest : public MeltTest {
    std::string dump_style = "netcdf";

public:
    void set_style(const std::string &new_style) { dump_style = new_style; }

    void enable_triclinic()
    {
        BEGIN_HIDE_OUTPUT();
        command("change_box all triclinic");
        END_HIDE_OUTPUT();
    }

    std::string dump_filename(std::string ident)
    {
        return fmt::format("dump_{}_{}.nc", dump_style, ident);
    }

    void generate_dump(std::string dump_file, std::string fields, std::string dump_modify_options,
                       int ntimesteps)
    {
        BEGIN_HIDE_OUTPUT();
        command(fmt::format("dump id all {} 1 {} {}", dump_style, dump_file, fields));

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

    void close_dump()
    {
        BEGIN_HIDE_OUTPUT();
        command("undump id");
        END_HIDE_OUTPUT();
    }

    std::string convert_binary_to_text(std::string binary_file)
    {
        BEGIN_HIDE_OUTPUT();
        std::string cmdline = fmt::format("{0} {1} > {1}.txt", NCDUMP_BINARY, binary_file);
        system(cmdline.c_str());
        END_HIDE_OUTPUT();
        return fmt::format("{}.txt", binary_file);
    }
};

TEST_F(DumpNetCDFTest, run0_plain)
{
    if (!lammps_has_style(lmp, "dump", "netcdf")) GTEST_SKIP();
    auto dump_file = dump_filename("run0");
    auto fields    = "id type proc procp1 mass x y z ix iy iz xu yu zu vx vy vz fx fy fz";
    set_style("netcdf");
    generate_dump(dump_file, fields, "", 0);

    ASSERT_FILE_EXISTS(dump_file);
    if (NCDUMP_BINARY) {
        auto converted_file = convert_binary_to_text(dump_file);
        auto lines          = read_lines(converted_file);
        auto words          = utils::split_words(lines[0]);
        ASSERT_EQ(lines.size(), 233);
        ASSERT_THAT(words[0], Eq("netcdf"));
        ASSERT_THAT(words[1]+".nc", Eq(dump_file));
        words = utils::split_words(lines[3]);
        ASSERT_THAT(words[0], Eq("atom"));
        ASSERT_THAT(words[2], Eq("32"));
        delete_file(converted_file);
    }
    delete_file(dump_file);
}

TEST_F(DumpNetCDFTest, run0_mpi)
{
    if (!lammps_has_style(lmp, "dump", "netcdf/mpiio")) GTEST_SKIP();
    auto dump_file = dump_filename("run0");
    auto fields    = "id type proc procp1 mass x y z ix iy iz xu yu zu vx vy vz fx fy fz";
    set_style("netcdf/mpiio");
    generate_dump(dump_file, fields, "", 0);

    ASSERT_FILE_EXISTS(dump_file);
    if (NCDUMP_BINARY) {
        auto converted_file = convert_binary_to_text(dump_file);
        auto lines          = read_lines(converted_file);
        auto words          = utils::split_words(lines[0]);
        ASSERT_EQ(lines.size(), 234);
        ASSERT_THAT(words[0], Eq("netcdf"));
        ASSERT_THAT(words[1]+".nc", Eq(dump_file));
        words = utils::split_words(lines[3]);
        ASSERT_THAT(words[0], Eq("atom"));
        ASSERT_THAT(words[2], Eq("32"));
        delete_file(converted_file);
    }
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

    NCDUMP_BINARY = getenv("NCDUMP_BINARY");

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
