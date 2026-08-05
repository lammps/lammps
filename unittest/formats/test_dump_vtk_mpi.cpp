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

// Tests for the parts of the "vtk" dump style that need more than one MPI
// process: the per processor piece files, the parallel summary file that
// refers to them, and collecting a complete snapshot in a single file.

#include "../testing/core.h"
#include "../testing/systems/melt.h"
#include "../testing/utils.h"
#include "comm.h"
#include "fmt/format.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

#include "../testing/test_mpi_main.h"

using ::testing::HasSubstr;
using ::testing::Not;

namespace LAMMPS_NS {

static std::string slurp(const std::string &filename)
{
    std::ifstream in(filename, std::ios::binary);
    std::ostringstream buf;
    buf << in.rdbuf();
    return buf.str();
}

// read the value of an integer XML attribute, e.g. NumberOfPoints="8"

static int xml_attribute(const std::string &text, const std::string &name)
{
    const auto pos = text.find(name + "=\"");
    if (pos == std::string::npos) return -1;
    return atoi(text.c_str() + pos + name.size() + 2);
}

class DumpVTKParallelTest : public MeltTest {
public:
    int me() const { return lmp->comm->me; }
    int nprocs() const { return lmp->comm->nprocs; }

    void generate_dump(const std::string &dump_file, const std::string &dump_modify_options)
    {
        BEGIN_HIDE_OUTPUT();
        command(fmt::format("dump id all vtk 1 {} id type", dump_file));

        if (!dump_modify_options.empty()) {
            command(fmt::format("dump_modify id {}", dump_modify_options));
        }

        command("run 0 post no");
        END_HIDE_OUTPUT();

        // the piece files are written by their own processor, so wait for all
        // of them before looking at any file that is not our own

        MPI_Barrier(MPI_COMM_WORLD);
    }
};

TEST_F(DumpVTKParallelTest, polydata_pieces)
{
    generate_dump("dump_vtk_par_%.vtp", "");

    // every processor writes its own piece and the pieces together hold the
    // complete system

    const auto piece_file = fmt::format("dump_vtk_par__{}.vtp", me());
    ASSERT_FILE_EXISTS(piece_file);

    int mypoints = xml_attribute(slurp(piece_file), "NumberOfPoints");
    ASSERT_GE(mypoints, 0);
    int npoints = 0;
    MPI_Allreduce(&mypoints, &npoints, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    ASSERT_EQ(npoints, 32);

    // the summary file lists one piece per processor

    if (me() == 0) {
        const std::string summary_file = "dump_vtk_par_.pvtp";
        ASSERT_FILE_EXISTS(summary_file);
        auto summary = slurp(summary_file);
        EXPECT_THAT(summary, HasSubstr(R"(<VTKFile type="PPolyData" version="1.0")"));
        for (int i = 0; i < nprocs(); ++i)
            EXPECT_THAT(summary,
                        HasSubstr(fmt::format(R"(<Piece Source="dump_vtk_par__{}.vtp"/>)", i)));
        EXPECT_THAT(summary, Not(HasSubstr(fmt::format(R"(dump_vtk_par__{}.vtp)", nprocs()))));

        // the box is written once, by the first processor

        ASSERT_FILE_EXISTS("dump_vtk_par__boundingBox.vtr");
        delete_file(summary_file);
        delete_file("dump_vtk_par__boundingBox.vtr");
    }
    delete_file(piece_file);
}

TEST_F(DumpVTKParallelTest, unstructured_pieces)
{
    generate_dump("dump_vtk_pgrid_%.vtu", "");

    const auto piece_file = fmt::format("dump_vtk_pgrid__{}.vtu", me());
    ASSERT_FILE_EXISTS(piece_file);
    EXPECT_THAT(slurp(piece_file), HasSubstr(R"(<VTKFile type="UnstructuredGrid")"));

    if (me() == 0) {
        const std::string summary_file = "dump_vtk_pgrid_.pvtu";
        ASSERT_FILE_EXISTS(summary_file);
        auto summary = slurp(summary_file);
        EXPECT_THAT(summary, HasSubstr(R"(<VTKFile type="PUnstructuredGrid" version="1.0")"));
        for (int i = 0; i < nprocs(); ++i)
            EXPECT_THAT(summary,
                        HasSubstr(fmt::format(R"(<Piece Source="dump_vtk_pgrid__{}.vtu"/>)", i)));
        delete_file(summary_file);
        delete_file("dump_vtk_pgrid__boundingBox.vtr");
    }
    delete_file(piece_file);
}

TEST_F(DumpVTKParallelTest, fewer_files_than_processors)
{
    if (nprocs() < 2) GTEST_SKIP();

    // with "nfile 2" only two processors write, each collecting the data of
    // its cluster of processors

    generate_dump("dump_vtk_nfile_%.vtp", "nfile 2");

    if (me() == 0) {
        int npoints = 0;
        for (int i = 0; i < 2; ++i) {
            const auto piece_file = fmt::format("dump_vtk_nfile__{}.vtp", i);
            ASSERT_FILE_EXISTS(piece_file);
            npoints += xml_attribute(slurp(piece_file), "NumberOfPoints");
            delete_file(piece_file);
        }
        ASSERT_EQ(npoints, 32);

        // no processor beyond the two file writers has written anything

        for (int i = 2; i < nprocs(); ++i)
            ASSERT_FILE_NOT_EXISTS(fmt::format("dump_vtk_nfile__{}.vtp", i));

        const std::string summary_file = "dump_vtk_nfile_.pvtp";
        ASSERT_FILE_EXISTS(summary_file);
        auto summary = slurp(summary_file);
        EXPECT_THAT(summary, HasSubstr(R"(<Piece Source="dump_vtk_nfile__0.vtp"/>)"));
        EXPECT_THAT(summary, HasSubstr(R"(<Piece Source="dump_vtk_nfile__1.vtp"/>)"));
        EXPECT_THAT(summary, Not(HasSubstr(R"(dump_vtk_nfile__2.vtp)")));
        delete_file(summary_file);
        delete_file("dump_vtk_nfile__boundingBox.vtr");
    }
}

TEST_F(DumpVTKParallelTest, legacy_collects_all_atoms)
{
    // the legacy format has no per processor output: the "%" becomes 0 and
    // the first processor writes all atoms of the system into that one file

    generate_dump("dump_vtk_legacy_%.vtk", "sort id");

    const std::string dump_file = "dump_vtk_legacy_0.vtk";
    const std::string box_file  = "dump_vtk_legacy_0_boundingBox.vtk";
    ASSERT_FILE_EXISTS(dump_file);
    ASSERT_FILE_EXISTS(box_file);
    for (int i = 1; i < nprocs(); ++i)
        ASSERT_FILE_NOT_EXISTS(fmt::format("dump_vtk_legacy_{}.vtk", i));

    // sorting is what puts the atoms of the different processors in order

    auto text = slurp(dump_file);
    EXPECT_THAT(text, HasSubstr("POINTS 32 float\n"));
    EXPECT_THAT(text, HasSubstr("id 1 32 int\n1 2 3 4 5 6 7 8 9\n"));

    // all processors read the same file, so it may only be removed once they
    // are done with it

    MPI_Barrier(MPI_COMM_WORLD);
    if (me() == 0) {
        delete_file(dump_file);
        delete_file(box_file);
    }
}
} // namespace LAMMPS_NS
