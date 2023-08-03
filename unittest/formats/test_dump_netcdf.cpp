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

#include "../testing/core.h"
#include "../testing/systems/melt.h"
#include "../testing/utils.h"
#include "fmt/format.h"
#include "library.h"
#include "output.h"
#include "thermo.h"
#include "utils.h"
#include "version.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <string>

using ::testing::Eq;

char *NCDUMP_EXECUTABLE = nullptr;
bool verbose            = false;

namespace LAMMPS_NS {
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
        std::string cmdline = fmt::format("{0} {1} > {1}.txt", NCDUMP_EXECUTABLE, binary_file);
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
    if (NCDUMP_EXECUTABLE) {
        auto converted_file = convert_binary_to_text(dump_file);
        auto lines          = read_lines(converted_file);
        auto header         = utils::split_words(lines[0]);
        ASSERT_EQ(lines.size(), 233);
        ASSERT_THAT(header[0], Eq("netcdf"));
        ASSERT_THAT(header[1] + ".nc", Eq(dump_file));

        // check dimensions section
        auto section = std::find(lines.begin(), lines.end(), "dimensions:");
        for (auto line = ++section; line < lines.end(); ++line) {
            auto words = utils::split_words(*line);
            if ((words.size() < 1) || (words[0] == "variables:")) break;
            if (words[0] == "atom") ASSERT_THAT(words[2], Eq("32"));
            if (words[0] == "label") ASSERT_THAT(words[2], Eq("10"));
            if (words[0] == "Voigt") ASSERT_THAT(words[2], Eq("6"));
            if (words[0] == "spatial") ASSERT_THAT(words[2], Eq("3"));
        }

        // check variables section
        section = std::find(lines.begin(), lines.end(), "variables:");
        for (auto line = ++section; line < lines.end(); ++line) {
            auto words = utils::split_words(*line);
            if ((words.size() < 2) || (words[0] == "data:")) break;
            if (words[0] == "time:units") ASSERT_THAT(words[2], Eq("lj"));
            if (words[0] == "time:scale_factor") ASSERT_THAT(words[2], Eq("0.005f"));
            if (words[0] == "cell_origin:units") ASSERT_THAT(words[2], Eq("lj"));
            if (words[0] == "cell_angles:units") ASSERT_THAT(words[2], Eq("degree"));
            if (words[1] == "id(frame,") ASSERT_THAT(words[2], Eq("atom)"));
            if (words[1] == "type(frame,") ASSERT_THAT(words[2], Eq("atom)"));
            if (words[1] == "proc(frame,") ASSERT_THAT(words[2], Eq("atom)"));
            if (words[1] == "procp1(frame,") ASSERT_THAT(words[2], Eq("atom)"));
            if (words[1] == "mass(frame,") ASSERT_THAT(words[2], Eq("atom)"));
            if (words[1] == "ix(frame,") ASSERT_THAT(words[2], Eq("atom)"));
            if (words[1] == "iy(frame,") ASSERT_THAT(words[2], Eq("atom)"));
            if (words[1] == "iz(frame,") ASSERT_THAT(words[2], Eq("atom)"));
            if (words[0] == ":Conventions") ASSERT_THAT(words[2], Eq("AMBER"));
            if (words[0] == ":ConventionVersion") ASSERT_THAT(words[2], Eq("1.0"));
            if (words[0] == ":program") ASSERT_THAT(words[2], Eq("LAMMPS"));
            if (words[0] == ":programVersion") ASSERT_THAT(words[2], Eq(LAMMPS_VERSION));
        }

        // check data section
        section = std::find(lines.begin(), lines.end(), "data:");
        for (auto line = ++section; line < lines.end(); ++line) {
            auto words = utils::split_words(*line);
            if (words.size() > 0) {
                if (words[0] == "spatial") ASSERT_THAT(words[2], Eq("xyz"));
                if (words[0] == "cell_spatial") ASSERT_THAT(words[2], Eq("abc"));
                if (words[0] == "cell_origin") {
                    ++line;
                    words = utils::split_words(*line);
                    ASSERT_THAT(words[0], Eq("0,"));
                    ASSERT_THAT(words[1], Eq("0,"));
                    ASSERT_THAT(words[2], Eq("0"));
                }
                if (words[0] == "cell_lengths") {
                    ++line;
                    words = utils::split_words(*line);
                    ASSERT_THAT(words[0], Eq("3.359192,"));
                    ASSERT_THAT(words[1], Eq("3.359192,"));
                    ASSERT_THAT(words[2], Eq("3.359192"));
                }
                if (words[0] == "cell_angles") {
                    ++line;
                    words = utils::split_words(*line);
                    ASSERT_THAT(words[0], Eq("90,"));
                    ASSERT_THAT(words[1], Eq("90,"));
                    ASSERT_THAT(words[2], Eq("90"));
                }
                if (words[0] == "id") {
                    ++line;
                    words = utils::split_words(*line);
                    ASSERT_THAT(words[0], Eq("1,"));
                    ASSERT_THAT(words[1], Eq("2,"));
                    ASSERT_THAT(words[2], Eq("3,"));
                    ASSERT_THAT(words[3], Eq("4,"));
                    ASSERT_THAT(words[4], Eq("5,"));
                    ASSERT_THAT(words[5], Eq("6,"));
                    ASSERT_THAT(words[6], Eq("7,"));
                    ASSERT_THAT(words[7], Eq("8,"));
                    ASSERT_THAT(words[8], Eq("9,"));
                    ASSERT_THAT(words[9], Eq("10,"));
                }
                if (words[0] == "mass") {
                    ++line;
                    words = utils::split_words(*line);
                    ASSERT_THAT(words[0], Eq("1,"));
                    ASSERT_THAT(words[1], Eq("1,"));
                    ASSERT_THAT(words[2], Eq("1,"));
                    ASSERT_THAT(words[3], Eq("1,"));
                    ASSERT_THAT(words[4], Eq("1,"));
                    ASSERT_THAT(words[5], Eq("1,"));
                    ASSERT_THAT(words[6], Eq("1,"));
                    ASSERT_THAT(words[7], Eq("1,"));
                    ASSERT_THAT(words[8], Eq("1,"));
                    ASSERT_THAT(words[9], Eq("1,"));
                }
                if (words[0] == "coordinates") {
                    ++line;
                    words = utils::split_words(*line);
                    ASSERT_THAT(words[0], Eq("0,"));
                    ASSERT_THAT(words[1], Eq("0,"));
                    ASSERT_THAT(words[2], Eq("0,"));
                    ++line;
                    words = utils::split_words(*line);
                    ASSERT_THAT(words[0], Eq("0.8397981,"));
                    ASSERT_THAT(words[1], Eq("0.8397981,"));
                    ASSERT_THAT(words[2], Eq("0,"));
                    ++line;
                    words = utils::split_words(*line);
                    ASSERT_THAT(words[0], Eq("0.8397981,"));
                    ASSERT_THAT(words[1], Eq("0,"));
                    ASSERT_THAT(words[2], Eq("0.8397981,"));
                    ++line;
                    words = utils::split_words(*line);
                    ASSERT_THAT(words[0], Eq("0,"));
                    ASSERT_THAT(words[1], Eq("0.8397981,"));
                    ASSERT_THAT(words[2], Eq("0.8397981,"));
                }
                if (words[0] == "ix") {
                    ++line;
                    words = utils::split_words(*line);
                    ASSERT_THAT(words[0], Eq("0,"));
                    ASSERT_THAT(words[1], Eq("0,"));
                    ASSERT_THAT(words[2], Eq("0,"));
                    ASSERT_THAT(words[3], Eq("0,"));
                    ASSERT_THAT(words[4], Eq("0,"));
                    ASSERT_THAT(words[5], Eq("0,"));
                    ASSERT_THAT(words[6], Eq("0,"));
                    ASSERT_THAT(words[7], Eq("0,"));
                    ASSERT_THAT(words[8], Eq("0,"));
                    ASSERT_THAT(words[9], Eq("0,"));
                }
            }
        }
        delete_file(converted_file);
    }
    delete_file(dump_file);
}

TEST_F(DumpNetCDFTest, run0_mpi)
{
    if (!lammps_has_style(lmp, "dump", "netcdf/mpiio")) GTEST_SKIP();
    auto dump_file = dump_filename("mpi0");
    auto fields    = "id type proc procp1 mass x y z ix iy iz xu yu zu vx vy vz fx fy fz";
    set_style("netcdf/mpiio");
    generate_dump(dump_file, fields, "", 0);

    ASSERT_FILE_EXISTS(dump_file);
    if (NCDUMP_EXECUTABLE) {
        auto converted_file = convert_binary_to_text(dump_file);
        auto lines          = read_lines(converted_file);
        auto header         = utils::split_words(lines[0]);
        ASSERT_EQ(lines.size(), 234);
        ASSERT_THAT(header[0], Eq("netcdf"));
        ASSERT_THAT(header[1] + ".nc", Eq(dump_file));

        // check dimensions section
        auto section = std::find(lines.begin(), lines.end(), "dimensions:");
        for (auto line = ++section; line < lines.end(); ++line) {
            auto words = utils::split_words(*line);
            if ((words.size() < 1) || (words[0] == "variables:")) break;
            if (words[0] == "atom") ASSERT_THAT(words[2], Eq("32"));
            if (words[0] == "label") ASSERT_THAT(words[2], Eq("10"));
            if (words[0] == "Voigt") ASSERT_THAT(words[2], Eq("6"));
            if (words[0] == "spatial") ASSERT_THAT(words[2], Eq("3"));
        }

        // check variables section
        section = std::find(lines.begin(), lines.end(), "variables:");
        for (auto line = ++section; line < lines.end(); ++line) {
            auto words = utils::split_words(*line);
            if ((words.size() < 2) || (words[0] == "data:")) break;
            if (words[0] == "time:units") ASSERT_THAT(words[2], Eq("lj"));
            if (words[0] == "time:scale_factor") ASSERT_THAT(words[2], Eq("0.005f"));
            if (words[0] == "cell_origin:units") ASSERT_THAT(words[2], Eq("lj"));
            if (words[0] == "cell_angles:units") ASSERT_THAT(words[2], Eq("degree"));
            if (words[1] == "id(frame,") ASSERT_THAT(words[2], Eq("atom)"));
            if (words[1] == "type(frame,") ASSERT_THAT(words[2], Eq("atom)"));
            if (words[1] == "proc(frame,") ASSERT_THAT(words[2], Eq("atom)"));
            if (words[1] == "procp1(frame,") ASSERT_THAT(words[2], Eq("atom)"));
            if (words[1] == "mass(frame,") ASSERT_THAT(words[2], Eq("atom)"));
            if (words[1] == "ix(frame,") ASSERT_THAT(words[2], Eq("atom)"));
            if (words[1] == "iy(frame,") ASSERT_THAT(words[2], Eq("atom)"));
            if (words[1] == "iz(frame,") ASSERT_THAT(words[2], Eq("atom)"));
            if (words[0] == ":Conventions") ASSERT_THAT(words[2], Eq("AMBER"));
            if (words[0] == ":ConventionVersion") ASSERT_THAT(words[2], Eq("1.0"));
            if (words[0] == ":program") ASSERT_THAT(words[2], Eq("LAMMPS"));
            if (words[0] == ":programVersion") ASSERT_THAT(words[2], Eq(LAMMPS_VERSION));
        }

        // check data section
        section = std::find(lines.begin(), lines.end(), "data:");
        for (auto line = ++section; line < lines.end(); ++line) {
            auto words = utils::split_words(*line);
            if (words.size() > 0) {
                if (words[0] == "spatial") ASSERT_THAT(words[2], Eq("xyz"));
                if (words[0] == "cell_spatial") ASSERT_THAT(words[2], Eq("abc"));
                if (words[0] == "cell_origin") {
                    ++line;
                    words = utils::split_words(*line);
                    ASSERT_THAT(words[0], Eq("0,"));
                    ASSERT_THAT(words[1], Eq("0,"));
                    ASSERT_THAT(words[2], Eq("0"));
                }
                if (words[0] == "cell_lengths") {
                    ++line;
                    words = utils::split_words(*line);
                    ASSERT_THAT(words[0], Eq("3.359192,"));
                    ASSERT_THAT(words[1], Eq("3.359192,"));
                    ASSERT_THAT(words[2], Eq("3.359192"));
                }
                if (words[0] == "cell_angles") {
                    ++line;
                    words = utils::split_words(*line);
                    ASSERT_THAT(words[0], Eq("90,"));
                    ASSERT_THAT(words[1], Eq("90,"));
                    ASSERT_THAT(words[2], Eq("90"));
                }
                if (words[0] == "id") {
                    ++line;
                    words = utils::split_words(*line);
                    ASSERT_THAT(words[0], Eq("1,"));
                    ASSERT_THAT(words[1], Eq("2,"));
                    ASSERT_THAT(words[2], Eq("3,"));
                    ASSERT_THAT(words[3], Eq("4,"));
                    ASSERT_THAT(words[4], Eq("5,"));
                    ASSERT_THAT(words[5], Eq("6,"));
                    ASSERT_THAT(words[6], Eq("7,"));
                    ASSERT_THAT(words[7], Eq("8,"));
                    ASSERT_THAT(words[8], Eq("9,"));
                    ASSERT_THAT(words[9], Eq("10,"));
                }
                if (words[0] == "mass") {
                    ++line;
                    words = utils::split_words(*line);
                    ASSERT_THAT(words[0], Eq("1,"));
                    ASSERT_THAT(words[1], Eq("1,"));
                    ASSERT_THAT(words[2], Eq("1,"));
                    ASSERT_THAT(words[3], Eq("1,"));
                    ASSERT_THAT(words[4], Eq("1,"));
                    ASSERT_THAT(words[5], Eq("1,"));
                    ASSERT_THAT(words[6], Eq("1,"));
                    ASSERT_THAT(words[7], Eq("1,"));
                    ASSERT_THAT(words[8], Eq("1,"));
                    ASSERT_THAT(words[9], Eq("1,"));
                }
                if (words[0] == "coordinates") {
                    ++line;
                    words = utils::split_words(*line);
                    ASSERT_THAT(words[0], Eq("0,"));
                    ASSERT_THAT(words[1], Eq("0,"));
                    ASSERT_THAT(words[2], Eq("0,"));
                    ++line;
                    words = utils::split_words(*line);
                    ASSERT_THAT(words[0], Eq("0.8397981,"));
                    ASSERT_THAT(words[1], Eq("0.8397981,"));
                    ASSERT_THAT(words[2], Eq("0,"));
                    ++line;
                    words = utils::split_words(*line);
                    ASSERT_THAT(words[0], Eq("0.8397981,"));
                    ASSERT_THAT(words[1], Eq("0,"));
                    ASSERT_THAT(words[2], Eq("0.8397981,"));
                    ++line;
                    words = utils::split_words(*line);
                    ASSERT_THAT(words[0], Eq("0,"));
                    ASSERT_THAT(words[1], Eq("0.8397981,"));
                    ASSERT_THAT(words[2], Eq("0.8397981,"));
                }
                if (words[0] == "ix") {
                    ++line;
                    words = utils::split_words(*line);
                    ASSERT_THAT(words[0], Eq("0,"));
                    ASSERT_THAT(words[1], Eq("0,"));
                    ASSERT_THAT(words[2], Eq("0,"));
                    ASSERT_THAT(words[3], Eq("0,"));
                    ASSERT_THAT(words[4], Eq("0,"));
                    ASSERT_THAT(words[5], Eq("0,"));
                    ASSERT_THAT(words[6], Eq("0,"));
                    ASSERT_THAT(words[7], Eq("0,"));
                    ASSERT_THAT(words[8], Eq("0,"));
                    ASSERT_THAT(words[9], Eq("0,"));
                }
            }
        }
        delete_file(converted_file);
    }
    delete_file(dump_file);
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

    NCDUMP_EXECUTABLE = getenv("NCDUMP_EXECUTABLE");

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
