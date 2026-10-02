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

// Tests for compute voronoi/atom, derived from the checks of the inputs in
// examples/voronoi on lattices small enough to run in a fraction of a
// second.  Where the examples compare sums against the box volume, these
// tests also check the individual cells against the known Voronoi
// polyhedra of the lattices.  The same executable is registered as a
// serial test and as a test with 4 MPI ranks: the tessellation is built
// per sub-domain from the owned plus ghost atoms, and the edge histogram
// and the occupation counts are summed across ranks.

#include "../testing/core.h"
#include "../testing/test_mpi_main.h"

#include "atom.h"
#include "comm.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "library.h"
#include "lmptype.h"
#include "variable.h"
#include "fmt/format.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <functional>
#include <map>
#include <mpi.h>
#include <string>
#include <utility>
#include <vector>

namespace LAMMPS_NS {

// one row of the local array: the atom IDs on either side of a Voronoi
// cell face (0 for a face with the boundary of the domain) and its area
struct Face {
    tagint itag, jtag;
    double area;
};

class VoronoiTest : public LAMMPSTest {
protected:
    // define the compute with the ID "v1" plus global sums of its first two
    // columns
    void define_compute(const std::string &group, const std::string &args)
    {
        HIDE_OUTPUT([&] {
            command(fmt::format("compute v1 {} voronoi/atom {}", group, args));
            command("compute vol all reduce sum c_v1[1]");
            command("compute faces all reduce sum c_v1[2]");
        });
    }

    // invoke the compute through a "run 0" that prints the sums (and any
    // further thermo keywords given) as thermo output.  the thermo output
    // re-invokes the compute on every call, also when the timestep did not
    // advance since the last call
    void run_thermo(const std::string &thermo = "")
    {
        HIDE_OUTPUT([&] {
            command("thermo_style custom step c_vol c_faces " + thermo);
        });
        evaluate();
    }

    void setup_compute(const std::string &group, const std::string &args,
                       const std::string &thermo = "")
    {
        define_compute(group, args);
        run_thermo(thermo);
    }

    void evaluate()
    {
        HIDE_OUTPUT([&] {
            command("run 0 post no");
        });
    }

    void remove_compute()
    {
        HIDE_OUTPUT([&] {
            command("uncompute v1");
            command("uncompute vol");
            command("uncompute faces");
        });
    }

    // check that a command raises the given error.  the TEST_FAILURE macro
    // matches the message in the captured screen output, which only MPI
    // rank 0 prints, so the message is matched in the exception instead,
    // which all ranks receive.  the check must also not return early on a
    // mismatch, since the ranks would then diverge in the collective calls
    // that follow
    void expect_error(const std::string &errmsg, const std::function<void()> &f)
    {
#if defined(LAMMPS_SKIP_DEATH_TESTS)
        (void)errmsg;
        (void)f;
#else
        std::string mesg;
        BEGIN_HIDE_OUTPUT();
        try {
            f();
        } catch (LAMMPSException &e) {
            mesg = e.what();
        }
        END_HIDE_OUTPUT();
        EXPECT_THAT(mesg, ContainsRegex(errmsg));
#endif
    }

    double scalar(const std::string &id)
    {
        auto *ptr = (double *)lammps_extract_compute(lmp, id.c_str(), LMP_STYLE_GLOBAL,
                                                     LMP_TYPE_SCALAR);
        return ptr ? *ptr : NAN;
    }

    double formula(const std::string &expression)
    {
        return lmp->input->variable->compute_equal(expression);
    }

    // the global vector of the compute, i.e. the histogram of the number of
    // edges per face, summed over all MPI ranks
    double *histogram()
    {
        return (double *)lammps_extract_compute(lmp, "v1", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR);
    }

    double **peratom(const std::string &id = "v1")
    {
        return (double **)lammps_extract_compute(lmp, id.c_str(), LMP_STYLE_ATOM, LMP_TYPE_ARRAY);
    }

    // the rows of the local array collected from all MPI ranks
    std::vector<Face> faces()
    {
        auto *nrows = (int *)lammps_extract_compute(lmp, "v1", LMP_STYLE_LOCAL, LMP_SIZE_ROWS);
        auto **rows = (double **)lammps_extract_compute(lmp, "v1", LMP_STYLE_LOCAL, LMP_TYPE_ARRAY);
        std::vector<double> mine;
        if (nrows && rows) {
            for (int i = 0; i < *nrows; ++i) {
                mine.push_back(rows[i][0]);
                mine.push_back(rows[i][1]);
                mine.push_back(rows[i][2]);
            }
        }

        const int nprocs = lmp->comm->nprocs;
        std::vector<int> counts(nprocs), offsets(nprocs);
        int mycount = (int)mine.size();
        MPI_Allgather(&mycount, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
        int total = 0;
        for (int i = 0; i < nprocs; ++i) {
            offsets[i] = total;
            total += counts[i];
        }
        // a rank without rows must still pass a valid buffer
        std::vector<double> all(total + 1);
        mine.push_back(0.0);
        MPI_Allgatherv(mine.data(), mycount, MPI_DOUBLE, all.data(), counts.data(), offsets.data(),
                       MPI_DOUBLE, MPI_COMM_WORLD);

        std::vector<Face> list;
        for (int i = 0; i < total; i += 3)
            list.push_back({(tagint)all[i], (tagint)all[i + 1], all[i + 2]});
        return list;
    }

    // expect every owned atom (of the given type, if non-zero) to have a
    // cell with the given volume and number of faces
    void expect_cells(double volume, int nfaces, int itype = 0, const std::string &id = "v1")
    {
        auto **cell = peratom(id);
        ASSERT_NE(cell, nullptr);
        const int *type   = lmp->atom->type;
        const tagint *tag = lmp->atom->tag;
        for (int i = 0; i < lmp->atom->nlocal; ++i) {
            if (itype && (type[i] != itype)) continue;
            EXPECT_NEAR(cell[i][0], volume, 1.0e-10) << "volume of atom " << tag[i];
            EXPECT_DOUBLE_EQ(cell[i][1], nfaces) << "faces of atom " << tag[i];
        }
    }

    // expect the given number of faces per number of edges, and none for
    // any other number of edges.  the last bin of the histogram collects
    // the faces with more edges than the maximum
    void expect_histogram(int maxedge, const std::map<int, double> &expected)
    {
        auto *histo = histogram();
        ASSERT_NE(histo, nullptr);
        for (int i = 0; i <= maxedge; ++i) {
            const int nedges = i + 1;
            const double count = expected.count(nedges) ? expected.at(nedges) : 0.0;
            EXPECT_DOUBLE_EQ(histo[i], count) << "faces with " << nedges << " edges";
        }
    }
};

// 6x6x6 unit cells of an fcc lattice with lattice constant 1 (864 atoms).
// the corner atoms of the unit cells are type 2 and form a simple cubic
// lattice, the face centered atoms are type 1 (the Cu3Au structure).  the
// Voronoi cell of every atom is a rhombic dodecahedron of volume 1/4 with
// 12 rhombic faces of area sqrt(2)/8 and edges of length sqrt(3)/4.  an
// atom map is required by the occupation keyword

class ComputeVoronoiTest : public VoronoiTest {
protected:
    static constexpr int natoms      = 864;
    static constexpr double boxvol   = 216.0;
    static constexpr double cellvol  = 0.25;
    static constexpr double facearea = 0.17677669529663687;    // sqrt(2)/8

    void SetUp() override
    {
        testbinary = "ComputeVoronoiTest";
        LAMMPSTest::SetUp();
        if (!Info::has_package("VORONOI")) GTEST_SKIP();
        HIDE_OUTPUT([&] {
            command("units metal");
            command("atom_style atomic");
            command("atom_modify map array");
            command("boundary p p p");
            command("lattice fcc 1.0 origin 0.25 0.25 0.25");
            command("region box block 0 6 0 6 0 6");
            command("create_box 2 box");
            command("create_atoms 1 box basis 1 2");
            command("mass * 1.0");
            command("pair_style zero 2.0");
            command("pair_coeff * *");
            command("neighbor 0.5 bin");
            command("group type1 type 1");
            command("group type2 type 2");
        });
        ASSERT_EQ(lmp->atom->natoms, natoms);
    }
};

TEST_F(ComputeVoronoiTest, cells)
{
    setup_compute("all", "");
    EXPECT_NEAR(scalar("vol"), boxvol, 1.0e-8);
    EXPECT_DOUBLE_EQ(scalar("faces"), 12.0 * natoms);
    expect_cells(cellvol, 12);
}

TEST_F(ComputeVoronoiTest, only_group)
{
    // the tessellation of the type 2 atoms alone is the simple cubic
    // lattice with cubic cells; the type 1 atoms get no cell at all
    define_compute("type2", "only_group");
    HIDE_OUTPUT([&] {
        command("compute vol1 type1 reduce sum c_v1[1]");
        command("compute vol2 type2 reduce sum c_v1[1]");
    });
    run_thermo("c_vol1 c_vol2");
    EXPECT_NEAR(scalar("vol2"), boxvol, 1.0e-8);
    EXPECT_DOUBLE_EQ(scalar("vol1"), 0.0);
    EXPECT_DOUBLE_EQ(scalar("faces"), 6.0 * natoms / 4);
    expect_cells(1.0, 6, 2);
    expect_cells(0.0, 0, 1);
}

TEST_F(ComputeVoronoiTest, group_without_only_group)
{
    // without only_group the tessellation includes all atoms and only the
    // output is restricted to the group
    setup_compute("type2", "");
    EXPECT_NEAR(scalar("vol"), boxvol / 4, 1.0e-8);
    EXPECT_DOUBLE_EQ(scalar("faces"), 12.0 * natoms / 4);
    expect_cells(cellvol, 12, 2);
    expect_cells(0.0, 0, 1);
}

TEST_F(ComputeVoronoiTest, radius)
{
    // radical Voronoi tessellation: the face between two atoms with radii
    // r1 and r2 at distance d is moved from the midpoint towards the
    // smaller atom by (r1^2 - r2^2)/(2 d).  with r = 0.3 for the type 2
    // atoms and r = 0.1 for the type 1 atoms this scales the dodecahedron
    // of a type 2 atom by s = 1 + (r2^2 - r1^2)/d^2 = 1.16, which moves its
    // 4-fold vertices to 0.58 from the atom, past the faces with the type 2
    // atoms at distance 1 along the cube axes; those cut off 6 pyramids
    // with a square base of volume (4/3) h^3 each, with h = 0.58 - 0.5.
    // all type 1 cells are equivalent and share the remaining volume
    const double s    = 1.16;
    const double vol2 = cellvol * s * s * s - 6.0 * (4.0 / 3.0) * std::pow(0.5 * s - 0.5, 3);
    const double vol1 = (4.0 * cellvol - vol2) / 3.0;

    HIDE_OUTPUT([&] {
        command("variable r atom (type==1)*0.1+(type==2)*0.3");
    });
    setup_compute("all", "radius v_r");
    EXPECT_NEAR(scalar("vol"), boxvol, 1.0e-8);
    expect_cells(vol2, 18, 2);
    expect_cells(vol1, 12, 1);
}

TEST_F(ComputeVoronoiTest, edge_histogram_fcc)
{
    // all 12 faces of a rhombic dodecahedron have 4 edges
    setup_compute("all", "edge_histo 8", "c_v1[4]");
    EXPECT_DOUBLE_EQ(scalar("faces"), 12.0 * natoms);
    expect_histogram(8, {{4, 12.0 * natoms}});
}

TEST_F(ComputeVoronoiTest, edge_histogram_group)
{
    // the histogram only includes the faces of the atoms in the group
    setup_compute("type2", "edge_histo 8", "c_v1[4]");
    expect_histogram(8, {{4, 12.0 * natoms / 4}});
}

TEST_F(ComputeVoronoiTest, edge_histogram_simple_cubic)
{
    // cubes: 6 faces with 4 edges each
    setup_compute("type2", "only_group edge_histo 8", "c_v1[4]");
    expect_histogram(8, {{4, 6.0 * natoms / 4}});
}

TEST_F(ComputeVoronoiTest, edge_histogram_perturbed)
{
    // small random displacements split the 4-fold vertices of the
    // dodecahedra into very short edges and tiny faces.  without a
    // threshold every face lands in some bin of the histogram; with an
    // edge length threshold the short edges are not counted, the tiny
    // faces are dropped and the histogram is that of the ideal lattice
    HIDE_OUTPUT([&] {
        command("displace_atoms all random 0.01 0.01 0.01 31423");
    });
    setup_compute("all", "edge_histo 8", "c_v1[4]");
    EXPECT_NEAR(scalar("vol"), boxvol, 1.0e-8);
    auto *histo = histogram();
    ASSERT_NE(histo, nullptr);
    double total = 0.0;
    for (int i = 0; i <= 8; ++i) total += histo[i];
    EXPECT_DOUBLE_EQ(total, scalar("faces"));

    remove_compute();
    setup_compute("all", "edge_histo 8 edge_threshold 0.1", "c_v1[4]");
    expect_histogram(8, {{4, 12.0 * natoms}});
}

TEST_F(ComputeVoronoiTest, edge_histogram_bcc)
{
    // the cell of a bcc lattice is a truncated octahedron of volume 1/2
    // with 6 square and 8 hexagonal faces
    HIDE_OUTPUT([&] {
        command("delete_atoms group all");
        command("lattice bcc 1.0 origin 0.25 0.25 0.25");
        command("create_atoms 1 box");
    });
    ASSERT_EQ(lmp->atom->natoms, natoms / 2);
    setup_compute("all", "edge_histo 8", "c_v1[4] c_v1[6]");
    EXPECT_NEAR(scalar("vol"), boxvol, 1.0e-8);
    expect_cells(0.5, 14);
    expect_histogram(8, {{4, 6.0 * natoms / 2}, {6, 8.0 * natoms / 2}});
}

TEST_F(ComputeVoronoiTest, face_threshold)
{
    // faces below the area threshold are not counted and not listed
    setup_compute("all", "face_threshold 0.1 neighbors yes");
    expect_cells(cellvol, 12);
    EXPECT_EQ(faces().size(), (std::size_t)(12 * natoms));

    remove_compute();
    setup_compute("all", "face_threshold 0.2 neighbors yes");
    expect_cells(cellvol, 0);
    EXPECT_EQ(faces().size(), (std::size_t)0);
}

TEST_F(ComputeVoronoiTest, surface)
{
    // the surface area of a rhombic dodecahedron with edge length e is
    // 8 sqrt(2) e^2.  a type 1 atom shares 4 of its 12 faces with type 2
    // atoms, a type 2 atom shares none with other type 2 atoms
    const double area = 8.0 * std::sqrt(2.0) * 3.0 / 16.0;
    define_compute("all", "surface all");
    HIDE_OUTPUT([&] {
        command("compute v2 all voronoi/atom surface type2");
        command("compute stot all reduce sum c_v1[3]");
        command("compute sgrp all reduce sum c_v2[3]");
    });
    run_thermo("c_stot c_sgrp");
    EXPECT_NEAR(scalar("stot"), area * natoms, 1.0e-8);
    EXPECT_NEAR(scalar("sgrp"), 4.0 * facearea * 3 * natoms / 4, 1.0e-8);

    auto **all = peratom("v1");
    auto **grp = peratom("v2");
    ASSERT_NE(all, nullptr);
    ASSERT_NE(grp, nullptr);
    const int *type = lmp->atom->type;
    for (int i = 0; i < lmp->atom->nlocal; ++i) {
        EXPECT_NEAR(all[i][2], area, 1.0e-10);
        if (type[i] == 1)
            EXPECT_NEAR(grp[i][2], 4.0 * facearea, 1.0e-10);
        else
            EXPECT_DOUBLE_EQ(grp[i][2], 0.0);
    }
}

TEST_F(ComputeVoronoiTest, neighbors)
{
    setup_compute("all", "neighbors yes");
    auto list = faces();
    ASSERT_EQ(list.size(), (std::size_t)(12 * natoms));

    // every face is between two atoms and has the same area
    std::map<std::pair<tagint, tagint>, double> area;
    for (const auto &face : list) {
        EXPECT_GT(face.itag, 0);
        EXPECT_GT(face.jtag, 0);
        EXPECT_NEAR(face.area, facearea, 1.0e-10);
        area[{face.itag, face.jtag}] = face.area;
    }
    ASSERT_EQ(area.size(), list.size());

    // and is listed once from either side with the same area
    for (const auto &[pair, value] : area) {
        auto other = area.find({pair.second, pair.first});
        ASSERT_NE(other, area.end()) << "face " << pair.first << "-" << pair.second;
        EXPECT_NEAR(value, other->second, 1.0e-10);
    }
}

TEST_F(ComputeVoronoiTest, occupation)
{
    // the occupation keyword builds the cells once and then reports for
    // every atom how many atoms are in the cell it originally occupied
    // (column 1) and how many atoms are in the cell it is currently in
    // (column 2).  the source atom is the type 2 atom at (1.25,1.25,1.25),
    // its destination the type 1 nearest neighbor at (1.75,1.75,1.25)
    HIDE_OUTPUT([&] {
        command("variable ids atom id");
        command("region src sphere 1.25 1.25 1.25 0.1 units box");
        command("region dst sphere 1.75 1.75 1.25 0.1 units box");
        command("group src region src");
        command("group dst region dst");
        command("compute srcid src reduce max v_ids");
        command("compute dstid dst reduce max v_ids");
    });
    setup_compute("all", "occupation", "c_srcid c_dstid");
    ASSERT_DOUBLE_EQ(formula("count(src)"), 1.0);
    ASSERT_DOUBLE_EQ(formula("count(dst)"), 1.0);
    const auto src = (tagint)scalar("srcid");
    const auto dst = (tagint)scalar("dstid");
    ASSERT_NE(src, dst);

    // one atom per cell
    EXPECT_DOUBLE_EQ(scalar("vol"), natoms);
    EXPECT_DOUBLE_EQ(scalar("faces"), natoms);
    expect_cells(1.0, 1);

    // move the source atom next to the destination atom: a vacancy at the
    // source site and two atoms sharing the destination cell
    HIDE_OUTPUT([&] {
        command(fmt::format("set atom {} x 1.95 y 1.95 z 1.45", src));
    });
    evaluate();
    EXPECT_DOUBLE_EQ(scalar("vol"), natoms);
    EXPECT_DOUBLE_EQ(scalar("faces"), natoms + 2);
    EXPECT_DOUBLE_EQ(formula(fmt::format("C_v1[{}][1]", src)), 0.0);
    EXPECT_DOUBLE_EQ(formula(fmt::format("C_v1[{}][2]", src)), 2.0);
    EXPECT_DOUBLE_EQ(formula(fmt::format("C_v1[{}][1]", dst)), 2.0);
    EXPECT_DOUBLE_EQ(formula(fmt::format("C_v1[{}][2]", dst)), 2.0);

    auto **cell      = peratom();
    const tagint *tag = lmp->atom->tag;
    for (int i = 0; i < lmp->atom->nlocal; ++i) {
        if ((tag[i] == src) || (tag[i] == dst)) continue;
        EXPECT_DOUBLE_EQ(cell[i][0], 1.0) << "atom " << tag[i];
        EXPECT_DOUBLE_EQ(cell[i][1], 1.0) << "atom " << tag[i];
    }

    // and back to the initial state
    HIDE_OUTPUT([&] {
        command(fmt::format("set atom {} x 1.25 y 1.25 z 1.25", src));
    });
    evaluate();
    EXPECT_DOUBLE_EQ(scalar("vol"), natoms);
    EXPECT_DOUBLE_EQ(scalar("faces"), natoms);
    expect_cells(1.0, 1);
}

TEST_F(ComputeVoronoiTest, triclinic)
{
    // a sheared lattice is still a lattice, so all cells remain equal
    HIDE_OUTPUT([&] {
        command("change_box all triclinic");
        command("change_box all xy final 2.0 remap units box");
    });
    setup_compute("all", "");
    EXPECT_NEAR(scalar("vol"), boxvol, 1.0e-8);
    auto **cell = peratom();
    ASSERT_NE(cell, nullptr);
    for (int i = 0; i < lmp->atom->nlocal; ++i)
        EXPECT_NEAR(cell[i][0], cellvol, 1.0e-10) << "atom " << lmp->atom->tag[i];
}

TEST_F(ComputeVoronoiTest, repeated_evaluation)
{
    // a second "run 0" without moving any atoms must give the same
    // histogram, also when nothing else invokes the per-atom data first
    HIDE_OUTPUT([&] {
        command("compute v1 all voronoi/atom edge_histo 8");
        command("thermo_style custom step c_v1[4]");
    });
    evaluate();
    expect_histogram(8, {{4, 12.0 * natoms}});
    evaluate();
    expect_histogram(8, {{4, 12.0 * natoms}});
}

TEST_F(ComputeVoronoiTest, histogram_invoked_twice)
{
    // fix ave/time invokes the histogram at the end of the step and the
    // thermo output invokes it again on the same step.  the second call
    // must not sum the counts across the MPI ranks a second time
    HIDE_OUTPUT([&] {
        command("compute v1 all voronoi/atom edge_histo 8");
        command("fix av all ave/time 1 1 1 c_v1[4]");
        command("thermo_style custom step c_v1[4] f_av");
        command("run 2 post no");
    });
    expect_histogram(8, {{4, 12.0 * natoms}});
    EXPECT_DOUBLE_EQ(formula("f_av"), 12.0 * natoms);
}

TEST_F(ComputeVoronoiTest, histogram_before_peratom)
{
    // the thermo output invokes the histogram first, and then the per-atom
    // data through the reduction.  the tessellation repeated for the
    // per-atom data must not replace the summed histogram with the counts
    // of the local rank
    HIDE_OUTPUT([&] {
        command("compute v1 all voronoi/atom edge_histo 8");
        command("compute vol all reduce sum c_v1[1]");
        command("thermo_style custom step c_v1[4] c_vol");
        command("run 1 post no");
    });
    expect_histogram(8, {{4, 12.0 * natoms}});
    EXPECT_NEAR(scalar("vol"), boxvol, 1.0e-8);
}

TEST_F(ComputeVoronoiTest, unsupported_settings)
{
    expect_error("ERROR: Illegal compute voronoi/atom command", [&] {
        command("compute v1 all voronoi/atom bogus");
    });
    expect_error("ERROR: Illegal compute voronoi/atom command", [&] {
        command("compute v1 all voronoi/atom radius r");
    });
    expect_error("ERROR: Illegal compute voronoi/atom command", [&] {
        command("compute v1 all voronoi/atom edge_histo");
    });
    expect_error("ERROR: Illegal compute voronoi/atom command .occupation and .surface or edges", [&] {
        command("compute v1 all voronoi/atom occupation edge_histo 8");
    });
    expect_error("ERROR: Could not find compute/voronoi surface group ID", [&] {
        command("compute v1 all voronoi/atom surface nosuchgroup");
    });

    // the radius variable is only looked up when the compute is invoked
    HIDE_OUTPUT([&] {
        command("compute v1 all voronoi/atom radius v_nosuchvariable");
        command("compute vol all reduce sum c_v1[1]");
        command("thermo_style custom step c_vol");
    });
    expect_error("ERROR: Variable name for voronoi radius does not exist", [&] {
        command("run 0 post no");
    });
}

TEST_F(ComputeVoronoiTest, occupation_requires_map)
{
    HIDE_OUTPUT([&] {
        command("clear");
        command("units metal");
        command("atom_style atomic");
        command("region box block 0 2 0 2 0 2");
        command("create_box 1 box");
    });
    expect_error("ERROR: Compute voronoi/atom occupation requires an atom map", [&] {
        command("compute v1 all voronoi/atom occupation");
    });
}

// 4x4 unit cells of a 2d hexagonal lattice with nearest neighbor distance
// 1 (32 atoms) in a box of height 1.  the cell of every atom is a
// hexagonal prism of volume sqrt(3)/2 with 6 rectangular side faces of
// area 1/sqrt(3) and two hexagonal faces with the boundaries in z

class ComputeVoronoi2dTest : public VoronoiTest {
protected:
    static constexpr int natoms      = 32;
    static constexpr double cellvol  = 0.86602540378443865;    // sqrt(3)/2
    static constexpr double sidearea = 0.57735026918962576;    // 1/sqrt(3)

    void SetUp() override
    {
        testbinary = "ComputeVoronoi2dTest";
        LAMMPSTest::SetUp();
        if (!Info::has_package("VORONOI")) GTEST_SKIP();
        HIDE_OUTPUT([&] {
            command("dimension 2");
            command("units metal");
            command("atom_style atomic");
            command("atom_modify map array");
            command("boundary p p p");
            command("lattice hex 1.0 origin 0.25 0.25 0.0");
            command("region box block 0 4 0 4 -0.5 0.5");
            command("region atoms block 0 4 0 4 0.0 0.0");
            command("create_box 1 box");
            command("create_atoms 1 region atoms");
            command("mass * 1.0");
            command("pair_style zero 2.0");
            command("pair_coeff * *");
            command("neighbor 0.5 bin");
        });
        ASSERT_EQ(lmp->atom->natoms, natoms);
    }
};

TEST_F(ComputeVoronoi2dTest, periodic)
{
    setup_compute("all", "neighbors yes edge_histo 6", "c_v1[4] c_v1[6]");
    const double boxvol = formula("lx*ly*lz");
    EXPECT_NEAR(scalar("vol"), boxvol, 1.0e-8);
    expect_cells(cellvol, 8);
    expect_histogram(6, {{4, 6.0 * natoms}, {6, 2.0 * natoms}});

    auto list = faces();
    ASSERT_EQ(list.size(), (std::size_t)(8 * natoms));
    int nboundary = 0;
    double total  = 0.0;
    for (const auto &face : list) {
        if (face.jtag == 0) {
            ++nboundary;
            EXPECT_NEAR(face.area, cellvol, 1.0e-10);
        } else {
            EXPECT_NEAR(face.area, sidearea, 1.0e-10);
        }
        total += face.area;
    }
    EXPECT_EQ(nboundary, 2 * natoms);
    EXPECT_NEAR(total, (6.0 * sidearea + 2.0 * cellvol) * natoms, 1.0e-8);
}

TEST_F(ComputeVoronoi2dTest, finite)
{
    // with non-periodic boundaries in x and y the cells are clipped at
    // the box boundary, so they still fill the box.  the histogram and
    // the total face area are the reference values from the logs of
    // examples/voronoi/in.voronoi.data, which sets up the same system
    HIDE_OUTPUT([&] {
        command("change_box all boundary f f p");
    });
    define_compute("all", "neighbors yes edge_histo 6");
    HIDE_OUTPUT([&] {
        command("compute area all reduce sum c_v1[3] inputs local");
    });
    run_thermo("c_v1[4] c_v1[5] c_v1[6] c_area");
    EXPECT_NEAR(scalar("vol"), formula("lx*ly*lz"), 1.0e-8);
    expect_histogram(6, {{4, 186.0}, {5, 12.0}, {6, 36.0}});
    EXPECT_DOUBLE_EQ(scalar("faces"), 186.0 + 12.0 + 36.0);
    EXPECT_NEAR(scalar("area"), 171.39013, 1.0e-5);

    auto **cell = peratom();
    ASSERT_NE(cell, nullptr);
    for (int i = 0; i < lmp->atom->nlocal; ++i)
        EXPECT_GT(cell[i][0], 0.0) << "atom " << lmp->atom->tag[i];

    auto list = faces();
    ASSERT_EQ(list.size(), (std::size_t)(186 + 12 + 36));
    double total = 0.0;
    for (const auto &face : list) total += face.area;
    EXPECT_NEAR(total, scalar("area"), 1.0e-10);
}

} // namespace LAMMPS_NS
