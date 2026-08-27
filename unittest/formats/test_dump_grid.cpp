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

// Tests for the "grid" and "grid/vtk" dump styles.  The system is the melt
// example with a 2x2x2 grid on top of it, like examples/grid/in.grid.3d, so
// that each grid cell holds the same number of atoms.

#include "../testing/core.h"
#include "../testing/systems/melt.h"
#include "../testing/utils.h"
#include "fmt/format.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using ::testing::Eq;
using ::testing::HasSubstr;
using ::testing::Not;
using ::testing::StartsWith;

bool verbose = false;

namespace LAMMPS_NS {

// read a whole file, dropping any '\r' so that assertions with '\n' literals
// also hold on Windows, where dump grid writes through a text mode FILE

static std::string slurp(const std::string &filename)
{
    std::ifstream in(filename, std::ios::binary);
    std::ostringstream buf;
    buf << in.rdbuf();
    auto rv = buf.str();
    rv.erase(std::remove(rv.begin(), rv.end(), '\r'), rv.end());
    return rv;
}

class DumpGridTest : public MeltTest {
protected:
    void InitSystem() override
    {
        MeltTest::InitSystem();

        // the 32 atoms of the melt system are spread over 8 grid cells, so
        // every cell contains 4 of them

        HIDE_OUTPUT([&] {
            command("fix ave all ave/grid 1 1 1 2 2 2 vx vy vz");
            command("compute prop all property/grid 2 2 2 id ix iy iz");
        });
    }

public:
    void enable_triclinic()
    {
        BEGIN_HIDE_OUTPUT();
        command("change_box all triclinic");
        END_HIDE_OUTPUT();
    }

    void generate_dump(const std::string &style, const std::string &dump_file,
                       const std::string &fields, const std::string &dump_modify_options,
                       int ntimesteps)
    {
        BEGIN_HIDE_OUTPUT();
        command(fmt::format("dump id all {} 1 {} {}", style, dump_file, fields));

        if (!dump_modify_options.empty()) {
            command(fmt::format("dump_modify id {}", dump_modify_options));
        }

        command(fmt::format("run {} post no", ntimesteps));
        END_HIDE_OUTPUT();
    }
};

//-------------------------------------------------------------------------------------------------
// dump grid
//-------------------------------------------------------------------------------------------------

TEST_F(DumpGridTest, run0)
{
    const std::string dump_file = "dump_grid_run0.melt";
    generate_dump("grid", dump_file, "c_prop:grid:data[*] f_ave:grid:count f_ave:grid:data[*]",
                  "sort 1", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 19);
    ASSERT_THAT(lines[0], Eq("ITEM: TIMESTEP"));
    ASSERT_THAT(lines[1], Eq("0"));
    ASSERT_THAT(lines[2], Eq("ITEM: BOX BOUNDS pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[3]).size(), 2);
    ASSERT_THAT(lines[6], Eq("ITEM: DIMENSION"));
    ASSERT_THAT(lines[7], Eq("3"));
    ASSERT_THAT(lines[8], Eq("ITEM: GRID SIZE nx ny nz"));
    ASSERT_THAT(lines[9], Eq("2 2 2"));
    ASSERT_THAT(lines[10], Eq("ITEM: GRID CELLS c_prop:grid:data[1] c_prop:grid:data[2] "
                              "c_prop:grid:data[3] c_prop:grid:data[4] f_ave:grid:count "
                              "f_ave:grid:data[1] f_ave:grid:data[2] f_ave:grid:data[3]"));

    // one line per grid cell, sorted by cell id.  the first four columns are
    // the cell id and its x, y, z index, the fifth is the atom count

    ASSERT_EQ(utils::split_words(lines[11]).size(), 8);
    ASSERT_THAT(lines[11], StartsWith("1 1 1 1 4 "));
    ASSERT_THAT(lines[12], StartsWith("2 2 1 1 4 "));
    ASSERT_THAT(lines[13], StartsWith("3 1 2 1 4 "));
    ASSERT_THAT(lines[18], StartsWith("8 2 2 2 4 "));

    delete_file(dump_file);
}

TEST_F(DumpGridTest, with_units_and_time_run0)
{
    const std::string dump_file = "dump_grid_units.melt";
    generate_dump("grid", dump_file, "f_ave:grid:count", "units yes time yes sort 1", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 23);
    ASSERT_THAT(lines[0], Eq("ITEM: UNITS"));
    ASSERT_THAT(lines[1], Eq("lj"));
    ASSERT_THAT(lines[2], Eq("ITEM: TIME"));
    ASSERT_THAT(lines[4], Eq("ITEM: TIMESTEP"));

    // the cells are written in id order, but the id itself is only a column
    // when it is requested as an attribute

    ASSERT_THAT(lines[14], Eq("ITEM: GRID CELLS f_ave:grid:count"));
    ASSERT_THAT(lines[15], Eq("4"));

    delete_file(dump_file);
}

TEST_F(DumpGridTest, triclinic_run0)
{
    const std::string dump_file = "dump_grid_tri.melt";
    enable_triclinic();
    generate_dump("grid", dump_file, "f_ave:grid:count", "sort 1", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 19);
    ASSERT_THAT(lines[2], Eq("ITEM: BOX BOUNDS xy xz yz pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[3]).size(), 3);

    delete_file(dump_file);
}

TEST_F(DumpGridTest, multi_file_run1)
{
    generate_dump("grid", "dump_grid_multi_*.melt", "f_ave:grid:count", "sort 1", 1);

    for (const auto &step : {"0", "1"}) {
        auto dump_file = fmt::format("dump_grid_multi_{}.melt", step);
        ASSERT_FILE_EXISTS(dump_file);
        ASSERT_EQ(count_lines(dump_file), 19);
        delete_file(dump_file);
    }
}

TEST_F(DumpGridTest, per_processor_file_run0)
{
    // the "%" becomes the id of the writing processor

    generate_dump("grid", "dump_grid_p%.melt", "f_ave:grid:count", "sort 1", 0);

    const std::string dump_file = "dump_grid_p0.melt";
    ASSERT_FILE_EXISTS(dump_file);
    ASSERT_EQ(count_lines(dump_file), 19);
    delete_file(dump_file);
}

TEST_F(DumpGridTest, no_attributes)
{
    TEST_FAILURE(".*No dump grid arguments specified.*",
                 command("dump id all grid 1 dump_grid_invalid.melt"););
}

TEST_F(DumpGridTest, invalid_attribute)
{
    TEST_FAILURE(".*Invalid attribute xxx in dump grid command.*",
                 command("dump id all grid 1 dump_grid_invalid.melt xxx"););
}

TEST_F(DumpGridTest, mismatched_grid_sizes)
{
    BEGIN_HIDE_OUTPUT();
    command("compute other all property/grid 3 3 3 id ix");
    command("dump id all grid 1 dump_grid_invalid.melt c_prop:grid:data[1] c_other:grid:data[1]");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*Dump grid field grid sizes do not match.*", command("run 0 post no"););
}

//-------------------------------------------------------------------------------------------------
// dump grid/vtk
//-------------------------------------------------------------------------------------------------

TEST_F(DumpGridTest, vtk_rectilinear_run0)
{
    const std::string dump_file = "dump_grid_rect.0.vtr";
    generate_dump("grid/vtk", "dump_grid_rect.*.vtr", "f_ave:grid:count", "", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto text = slurp(dump_file);
    EXPECT_THAT(text, HasSubstr(R"(<VTKFile type="RectilinearGrid" version="0.1")"));

    // the extent counts the points bounding the cells, so it is one more than
    // the number of cells in each direction

    EXPECT_THAT(text, HasSubstr(R"(<RectilinearGrid WholeExtent="0 2 0 2 0 2">)"));
    EXPECT_THAT(text, HasSubstr(R"(<Piece Extent="0 2 0 2 0 2">)"));

    // a single field becomes a scalar cell array, one value per grid cell

    EXPECT_THAT(text, HasSubstr(R"(<DataArray type="Float32" Name="Scalar" format="ascii">)"));
    EXPECT_THAT(text, HasSubstr("4 4 4 4 4 4 4 4 \n"));
    EXPECT_THAT(text, HasSubstr(R"(<DataArray type="Float32" Name="x" format="ascii">)"));

    delete_file(dump_file);
}

TEST_F(DumpGridTest, vtk_legacy_run0)
{
    const std::string dump_file = "dump_grid_legacy.0.vtk";
    generate_dump("grid/vtk", "dump_grid_legacy.*.vtk", "f_ave:grid:data[*]", "", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_GT(lines.size(), 4);
    ASSERT_THAT(lines[0], Eq("# vtk DataFile Version 5.1"));
    ASSERT_THAT(lines[2], Eq("ASCII"));

    auto text = slurp(dump_file);
    EXPECT_THAT(text, HasSubstr("DATASET RECTILINEAR_GRID\nDIMENSIONS 3 3 3\n"));
    EXPECT_THAT(text, HasSubstr("X_COORDINATES 3 float\n0 1.679"));

    // three fields become a vector cell array with three components

    EXPECT_THAT(text, HasSubstr("CELL_DATA 8\nFIELD FieldData 1\nVector 3 8 float\n"));

    delete_file(dump_file);
}

TEST_F(DumpGridTest, vtk_image_run0)
{
    const std::string dump_file = "dump_grid_image.0.vti";
    generate_dump("grid/vtk", "dump_grid_image.*.vti", "f_ave:grid:count", "", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto text = slurp(dump_file);
    EXPECT_THAT(text, HasSubstr(R"(<VTKFile type="ImageData" version="0.1")"));

    // all cells have the same size, so origin and spacing describe the grid

    EXPECT_THAT(text, HasSubstr(R"(<ImageData WholeExtent="0 2 0 2 0 2" Origin="0 0 0" )"
                                R"(Spacing="1.6795961913825073 1.6795961913825073 1.6795961913825073">)"));
    EXPECT_THAT(text, HasSubstr(R"(<DataArray type="Float32" Name="Scalar" format="ascii">)"));
    EXPECT_THAT(text, Not(HasSubstr("<Coordinates>")));

    delete_file(dump_file);
}

TEST_F(DumpGridTest, vtk_unknown_extension_run0)
{
    const std::string dump_file = "dump_grid_unknown.0.melt";
    generate_dump("grid/vtk", "dump_grid_unknown.*.melt", "f_ave:grid:count", "", 0);

    // an unknown file name extension keeps the XML rectilinear grid default

    ASSERT_FILE_EXISTS(dump_file);
    EXPECT_THAT(slurp(dump_file), HasSubstr(R"(<VTKFile type="RectilinearGrid")"));

    delete_file(dump_file);
}

TEST_F(DumpGridTest, vtk_binary_run0)
{
    const std::string xml_file    = "dump_grid_bin.0.vtr";
    const std::string legacy_file = "dump_grid_bin.0.vtk";

    BEGIN_HIDE_OUTPUT();
    command("dump xml all grid/vtk 1 dump_grid_bin.*.vtr f_ave:grid:count");
    command("dump_modify xml binary yes");
    command("dump legacy all grid/vtk 1 dump_grid_bin.*.vtk f_ave:grid:count");
    command("dump_modify legacy binary yes");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    ASSERT_FILE_EXISTS(xml_file);
    auto xml = slurp(xml_file);
    EXPECT_THAT(xml, HasSubstr(R"(<DataArray type="Float32" Name="Scalar" format="binary">)"));
    EXPECT_THAT(xml, Not(HasSubstr(R"(format="ascii")")));

    ASSERT_FILE_EXISTS(legacy_file);
    EXPECT_THAT(slurp(legacy_file), HasSubstr("BINARY\nDATASET RECTILINEAR_GRID\n"));

    delete_file(xml_file);
    delete_file(legacy_file);
}

TEST_F(DumpGridTest, vtk_double_precision_run0)
{
    const std::string xml_file    = "dump_grid_double.0.vtr";
    const std::string legacy_file = "dump_grid_double.0.vtk";

    // dump_modify double yes switches the grid coordinates and the cell data
    // to double precision

    BEGIN_HIDE_OUTPUT();
    command("dump xml all grid/vtk 1 dump_grid_double.*.vtr f_ave:grid:count");
    command("dump_modify xml double yes");
    command("dump legacy all grid/vtk 1 dump_grid_double.*.vtk f_ave:grid:data[*]");
    command("dump_modify legacy double yes");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    ASSERT_FILE_EXISTS(xml_file);
    auto xml = slurp(xml_file);
    EXPECT_THAT(xml, HasSubstr(R"(<DataArray type="Float64" Name="Scalar" format="ascii">)"));
    EXPECT_THAT(xml, HasSubstr(R"(<DataArray type="Float64" Name="x" format="ascii">)"));

    ASSERT_FILE_EXISTS(legacy_file);
    auto legacy = slurp(legacy_file);
    EXPECT_THAT(legacy, HasSubstr("X_COORDINATES 3 double\n0 1.6795961913825073"));
    EXPECT_THAT(legacy, HasSubstr("Vector 3 8 double\n"));

    delete_file(xml_file);
    delete_file(legacy_file);
}

TEST_F(DumpGridTest, vtk_needs_one_file_per_snapshot)
{
    BEGIN_HIDE_OUTPUT();
    command("dump id all grid/vtk 1 dump_grid_invalid.vtr f_ave:grid:count");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*Dump grid/vtk requires one snapshot per file.*", command("run 0 post no"););

    // the file is created before the dump style can complain about its name

    delete_file("dump_grid_invalid.vtr");
}

TEST_F(DumpGridTest, vtk_needs_sorting)
{
    BEGIN_HIDE_OUTPUT();
    command("dump id all grid/vtk 1 dump_grid_invalid.*.vtr f_ave:grid:count");
    command("dump_modify id sort off");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*Dump grid/vtk requires sorting on IDs.*", command("run 0 post no"););
}

TEST_F(DumpGridTest, vtk_rejects_triclinic)
{
    enable_triclinic();
    BEGIN_HIDE_OUTPUT();
    command("dump id all grid/vtk 1 dump_grid_invalid.*.vtr f_ave:grid:count");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*Dump grid/vtk does not support triclinic simulation boxes.*",
                 command("run 0 post no"););
}

TEST_F(DumpGridTest, vtk_invalid_filename)
{
    TEST_FAILURE(".*Invalid dump grid/vtk filename.*",
                 command("dump id all grid/vtk 1 dump_grid_invalid_%.*.vtr f_ave:grid:count"););
}

TEST_F(DumpGridTest, vtk_invalid_number_of_fields)
{
    TEST_FAILURE(".*Dump grid/vtk requires one or three fields.*",
                 command("dump id all grid/vtk 1 dump_grid_invalid.*.vtr f_ave:grid:count "
                         "f_ave:grid:data[1]"););
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
