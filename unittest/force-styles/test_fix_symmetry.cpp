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

#include "atom.h"
#include "domain.h"
#include "input.h"
#include "lammps.h"
#include "platform.h"

#include "../../src/EXTRA-FIX/symmetry_group.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

bool verbose = false;

using LAMMPS_NS::Atom;
using LAMMPS_NS::tagint;
using ::testing::Not;

// Symmetry-data file content for a tests that need a symfile on disk.
// Returns the absolute path to a newly written file in the working dir.
static std::string write_symfile(const std::string &name, const std::string &content)
{
    std::string path = name;
    std::ofstream out(path);
    out << content;
    out.close();
    return path;
}

// a space group operator s' = R s + t in fractional coordinates, as used by
// the test helpers below to check the orbit relations independently of the fix
struct TestOp {
    double R[3][3];
    double t[3];
};

static const TestOp OP_IDENTITY = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, {0, 0, 0}};
static const TestOp OP_INVERSION = {{{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}, {0, 0, 0}};

class FixSymmetryTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "FixSymmetryTest";
        LAMMPSTest::SetUp();
    }

    // 8-atom P-1 system in an orthogonal [0, L]^3 box. Atoms 1..4 are
    // asymmetric representatives; atoms 5..8 are their lattice-inversion
    // images (Cartesian: x_image = L - x_asym componentwise).
    void setup_pinv_system(double L = 8.0)
    {
        // HIDE_OUTPUT lambda form is exception-safe -- if a command throws,
        // the captured stdout is released before the exception propagates.
        // Plain BEGIN/END_HIDE_OUTPUT would leak the capture and break the
        // next test with "Only one stdout capturer can exist at a time".
        HIDE_OUTPUT([&] {
            command("units lj");
            command("atom_style atomic");
            command("boundary p p p");
            command("atom_modify map array");
            command(fmt::format("region box block 0 {} 0 {} 0 {}", L, L, L));
            command("create_box 1 box");
            command("create_atoms 1 single 1.2 1.5 1.8");
            command("create_atoms 1 single 3.0 1.2 0.9");
            command("create_atoms 1 single 1.5 3.0 1.1");
            command("create_atoms 1 single 0.9 1.1 3.2");
            command(fmt::format("create_atoms 1 single {} {} {}", L - 1.2, L - 1.5, L - 1.8));
            command(fmt::format("create_atoms 1 single {} {} {}", L - 3.0, L - 1.2, L - 0.9));
            command(fmt::format("create_atoms 1 single {} {} {}", L - 1.5, L - 3.0, L - 1.1));
            command(fmt::format("create_atoms 1 single {} {} {}", L - 0.9, L - 1.1, L - 3.2));
            command("mass 1 1.0");
            command("pair_style lj/cut 2.5");
            command("pair_coeff * * 1.0 1.0");
            command("velocity all create 0.7 12345 mom yes rot yes");
            command("fix 1 all nve");
        });
    }

    // Returns the Cartesian position of the atom with global tag t on rank 0
    // (the only rank in this serial test runner). Caller must ensure the
    // map is enabled; in this fixture it is via setup_pinv_system().
    void get_pos(tagint t, double xyz[3])
    {
        int li = lmp->atom->map(t);
        ASSERT_GE(li, 0);
        ASSERT_LT(li, lmp->atom->nlocal);
        xyz[0] = lmp->atom->x[li][0];
        xyz[1] = lmp->atom->x[li][1];
        xyz[2] = lmp->atom->x[li][2];
    }

    // largest deviation (in fractional coordinates, modulo lattice translations)
    // of the atoms tags[k] from the image R_k s + t_k of atom tags[0]
    double orbit_residual(const std::vector<tagint> &tags, const std::vector<TestOp> &ops)
    {
        double x0[3], s0[3];
        get_pos(tags[0], x0);
        lmp->domain->x2lamda(x0, s0);
        double maxres = 0.0;
        for (std::size_t k = 1; k < tags.size(); ++k) {
            double xk[3], sk[3];
            get_pos(tags[k], xk);
            lmp->domain->x2lamda(xk, sk);
            for (int i = 0; i < 3; ++i) {
                double d = sk[i] - ops[k].t[i];
                for (int j = 0; j < 3; ++j)
                    d -= ops[k].R[i][j] * s0[j];
                maxres = std::max(maxres, std::fabs(d - std::round(d)));
            }
        }
        return maxres;
    }

    // per-atom positions and velocities keyed by atom ID
    using Snapshot = std::map<tagint, std::array<double, 6>>;

    Snapshot take_snapshot()
    {
        Snapshot snap;
        const auto *atom = lmp->atom;
        for (int i = 0; i < atom->nlocal; ++i)
            snap[atom->tag[i]] = {atom->x[i][0], atom->x[i][1], atom->x[i][2],
                                  atom->v[i][0], atom->v[i][1], atom->v[i][2]};
        return snap;
    }

    void restore_snapshot(const Snapshot &snap)
    {
        auto *atom = lmp->atom;
        for (int i = 0; i < atom->nlocal; ++i) {
            const auto &d = snap.at(atom->tag[i]);
            for (int k = 0; k < 3; ++k) {
                atom->x[i][k] = d[k];
                atom->v[i][k] = d[k + 3];
            }
        }
    }

    // largest minimum image distance component between two snapshots
    double max_position_difference(const Snapshot &a, const Snapshot &b)
    {
        const double prd[3] = {lmp->domain->xprd, lmp->domain->yprd, lmp->domain->zprd};
        double maxdiff = 0.0;
        for (const auto &kv : a) {
            const auto &other = b.at(kv.first);
            for (int k = 0; k < 3; ++k) {
                double d = (kv.second[k] - other[k]) / prd[k];
                maxdiff = std::max(maxdiff, std::fabs(d - std::round(d)) * prd[k]);
            }
        }
        return maxdiff;
    }
};

// operators of frequently used groups in coordinate triplet notation
static const char *JSON_PINV = R"({ "group_name": "P-1", "ops": [ "x,y,z", "-x,-y,-z" ] })";
static const char *JSON_PM = R"({ "group_name": "Pm", "ops": [ "x,y,z", "x,-y,z" ] })";

// --------------------------------------------------------------------------
//   Operator parser: coordinate triplet notation to matrix and translation
// --------------------------------------------------------------------------
TEST(SymmetryOperatorParser, CoordinateTriplets)
{
    struct Case {
        const char *str;
        double R[3][3];
        double t[3];
    };
    const Case cases[] = {
        {"x,y,z", {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, {0, 0, 0}},
        {"-x,-y,-z", {{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}, {0, 0, 0}},
        {"-x+1/2, y, -z", {{-1, 0, 0}, {0, 1, 0}, {0, 0, -1}}, {0.5, 0, 0}},
        {"y,-x,z+1/4", {{0, 1, 0}, {-1, 0, 0}, {0, 0, 1}}, {0, 0, 0.25}},
        {"-y,x-y,z+1/3", {{0, -1, 0}, {1, -1, 0}, {0, 0, 1}}, {0, 0, 1.0 / 3.0}},
        {"1/2+X,1/2-Y,Z", {{1, 0, 0}, {0, -1, 0}, {0, 0, 1}}, {0.5, 0.5, 0}},
        {"x-y,x,z-0.5", {{1, -1, 0}, {1, 0, 0}, {0, 0, 1}}, {0, 0, -0.5}},
        {"2*x-y,y,z", {{2, -1, 0}, {0, 1, 0}, {0, 0, 1}}, {0, 0, 0}},
    };
    for (const auto &c : cases) {
        auto op = LAMMPS_NS::SymmetryGroup::parse_xyz(c.str);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j)
                EXPECT_DOUBLE_EQ(op.R[i][j], c.R[i][j]) << c.str << " R[" << i << "][" << j << "]";
            EXPECT_DOUBLE_EQ(op.t[i], c.t[i]) << c.str << " t[" << i << "]";
        }
    }

    const char *bad[] = {"x,y", "x,y,z,x", "x,,z", "x,y,a", "x,y,z+", "x,y,1/0", "x,y,1/2z",
                         "xy,y,z", "x,y,z*", "x,--y,z", "x,y,1/2/3"};
    for (const auto *str : bad)
        EXPECT_THROW(LAMMPS_NS::SymmetryGroup::parse_xyz(str), std::invalid_argument) << str;
}

// --------------------------------------------------------------------------
//   P1 (identity only) does not constrain anything
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, P1IsNoOp)
{
    std::string p1 = write_symfile("p1.json", R"({ "group_name": "P1", "ops": [ "x,y,z" ] })");

    setup_pinv_system(8.0);
    auto output = CAPTURE_OUTPUT([&] {
        command(fmt::format("fix sym all symmetry {}", p1));
        command("run 0 post no");
    });
    EXPECT_THAT(output, ContainsRegex("0 atoms in 0 orbits"));
    const auto start = take_snapshot();
    HIDE_OUTPUT([&] { command("run 20 post no"); });
    const auto with_fix = take_snapshot();

    restore_snapshot(start);
    HIDE_OUTPUT([&] {
        command("unfix sym");
        command("run 20 post no");
    });
    EXPECT_EQ(max_position_difference(with_fix, take_snapshot()), 0.0);

    LAMMPS_NS::platform::unlink(p1);
}

// --------------------------------------------------------------------------
//   P-1 inversion: the sum of the positions of each atom and its image
//   is a lattice vector throughout the run
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, PinvPreservesSymmetryUnderNVE)
{
    std::string sym = write_symfile("pinv.json", JSON_PINV);

    setup_pinv_system(8.0);
    auto output = CAPTURE_OUTPUT([&] {
        command(fmt::format("fix sym all symmetry {}", sym));
        command("run 50 post no");
    });
    EXPECT_THAT(output, ContainsRegex("Fix symmetry sym: group P-1 with 2 operators, 8 atoms in 4 "
                                      "orbits .0 on special positions. with 12 independent"));

    for (tagint a = 1; a <= 4; ++a)
        EXPECT_LT(orbit_residual({a, a + 4}, {OP_IDENTITY, OP_INVERSION}), 1.0e-12);

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   Validation independent of the implementation: starting from a symmetric
//   configuration with symmetric velocities, the exact dynamics stays
//   symmetric, so the constrained trajectory must match plain NVE (without
//   the fix) to round-off over a short run. Any error in the projection
//   (wrong normalization, missing contribution, wrong image refold) would
//   show up as a difference far above round-off within a few steps.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, EquivalentToPlainNVE)
{
    std::string sym = write_symfile("pinv_equiv.json", JSON_PINV);

    setup_pinv_system(8.0);
    HIDE_OUTPUT([&] {
        command("variable etot equal etotal");
        command(fmt::format("fix sym all symmetry {}", sym));
        // run 0 symmetrizes the random initial velocities
        command("run 0 post no");
    });
    const auto start = take_snapshot();

    HIDE_OUTPUT([&] { command("run 300 post no"); });
    const auto with_fix = take_snapshot();
    const double etot_with_fix = get_variable_value("etot");

    restore_snapshot(start);
    HIDE_OUTPUT([&] {
        command("unfix sym");
        command("run 300 post no");
    });
    const auto without_fix = take_snapshot();
    const double etot_without_fix = get_variable_value("etot");

    EXPECT_LT(max_position_difference(with_fix, without_fix), 1.0e-10);
    EXPECT_NEAR(etot_with_fix, etot_without_fix, 1.0e-10 * std::fabs(etot_without_fix));

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   Without the fix, round-off breaks the symmetry and it is lost after a
//   few thousand steps. This documents why the fix is needed at all.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, SymmetryLostWithoutFix)
{
    std::string sym = write_symfile("pinv_decay.json", JSON_PINV);

    setup_pinv_system(8.0);
    HIDE_OUTPUT([&] {
        command(fmt::format("fix sym all symmetry {}", sym));
        command("run 0 post no");
        command("unfix sym");
    });
    double start = 0.0;
    for (tagint a = 1; a <= 4; ++a)
        start = std::max(start, orbit_residual({a, a + 4}, {OP_IDENTITY, OP_INVERSION}));
    EXPECT_LT(start, 1.0e-12);

    HIDE_OUTPUT([&] { command("run 20000 post no"); });
    double end = 0.0;
    for (tagint a = 1; a <= 4; ++a)
        end = std::max(end, orbit_residual({a, a + 4}, {OP_IDENTITY, OP_INVERSION}));
    EXPECT_GT(end, 1.0e-6);

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   Three-fold axis in a hexagonal (prism) box: the operators have
//   off-diagonal integer matrices in fractional coordinates, so this checks
//   the fractional <-> Cartesian conversions of forces and velocities.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, HexagonalThreeFold)
{
    std::string sym =
        write_symfile("p3.json", R"({ "group_name": "P3", "ops": [ "x,y,z", "-y,x-y,z", "-x+y,-x,z" ] })");
    const std::vector<TestOp> ops = {OP_IDENTITY,
                                     {{{0, -1, 0}, {1, -1, 0}, {0, 0, 1}}, {0, 0, 0}},
                                     {{{-1, 1, 0}, {-1, 0, 0}, {0, 0, 1}}, {0, 0, 0}}};

    HIDE_OUTPUT([&] {
        command("units lj");
        command("atom_style atomic");
        command("boundary p p p");
        command("atom_modify map array");
        command("region box prism 0 4.0 0 3.4641016151377544 0 4.0 -2.0 0 0");
        command("create_box 1 box");
        // fractional (0.3,0.1,0.5) and (0.2,0.55,0.05) plus their C3 images
        command("create_atoms 1 single 1.0 0.346410161514 2.0 units box");
        command("create_atoms 1 single 3.2 0.692820323028 2.0 units box");
        command("create_atoms 1 single 1.8 2.424871130596 2.0 units box");
        command("create_atoms 1 single -0.5 1.905255888326 0.2 units box");
        command("create_atoms 1 single 0.6 2.078460969083 0.2 units box");
        command("create_atoms 1 single -0.1 2.944486372867 0.2 units box");
        command("mass 1 1.0");
        command("pair_style lj/cut 2.5");
        command("pair_coeff * * 1.0 1.0");
        command("pair_modify shift yes");
        command("velocity all create 0.3 4711 mom yes");
        command("timestep 0.002");
        command("variable etot equal etotal");
        command("fix 1 all nve");
        command(fmt::format("fix sym all symmetry {}", sym));
        command("run 0 post no");
    });
    // the input coordinates are rounded to 12 digits, the fix corrects that
    EXPECT_LT(orbit_residual({1, 2, 3}, ops), 1.0e-14);
    EXPECT_LT(orbit_residual({4, 5, 6}, ops), 1.0e-14);
    const double etot_start = get_variable_value("etot");

    HIDE_OUTPUT([&] { command("run 2000 post no"); });
    EXPECT_LT(orbit_residual({1, 2, 3}, ops), 1.0e-12);
    EXPECT_LT(orbit_residual({4, 5, 6}, ops), 1.0e-12);
    EXPECT_NEAR(get_variable_value("etot"), etot_start, 5.0e-4);

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   Two-fold screw axis along z: an operator with a nonzero translation.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, ScrewAxisNonzeroTranslation)
{
    std::string sym =
        write_symfile("p21.json", R"({ "group_name": "P2_1", "ops": [ "x,y,z", "-x,-y,z+1/2" ] })");
    const std::vector<TestOp> ops = {OP_IDENTITY,
                                     {{{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}}, {0, 0, 0.5}}};

    HIDE_OUTPUT([&] {
        command("units lj");
        command("atom_style atomic");
        command("boundary p p p");
        command("atom_modify map array");
        command("region box block 0 6 0 6 0 6");
        command("create_box 1 box");
        command("create_atoms 1 single 1.0 1.5 1.0");
        command("create_atoms 1 single 5.0 4.5 4.0");
        command("create_atoms 1 single 2.5 4.0 2.2");
        command("create_atoms 1 single 3.5 2.0 5.2");
        command("mass 1 1.0");
        command("pair_style lj/cut 2.5");
        command("pair_coeff * * 1.0 1.0");
        command("pair_modify shift yes");
        command("velocity all create 0.5 4711 mom yes");
        command("variable etot equal etotal");
        command("fix 1 all nve");
        command(fmt::format("fix sym all symmetry {}", sym));
        command("run 0 post no");
    });
    const double etot_start = get_variable_value("etot");

    HIDE_OUTPUT([&] { command("run 2000 post no"); });
    EXPECT_LT(orbit_residual({1, 2}, ops), 1.0e-12);
    EXPECT_LT(orbit_residual({3, 4}, ops), 1.0e-12);
    EXPECT_NEAR(get_variable_value("etot"), etot_start, 1.0e-3);

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   Screw axis in a box of two unit cells along z: the translation by one
//   unit cell must be part of the group, which the supercell entry adds.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, SupercellScrewAxis)
{
    std::string bad =
        write_symfile("p21_bad.json", R"({ "group_name": "P2_1", "ops": [ "x,y,z", "-x,-y,z+1/4" ] })");
    TEST_FAILURE("the operators do not form a group: the combination of operators 2 and 2 is "
                 "missing",
                 command("atom_modify map array");
                 command("region box block 0 6 0 6 0 12");
                 command("create_box 1 box");
                 command(fmt::format("fix sym all symmetry {}", bad)););
    LAMMPS_NS::platform::unlink(bad);

    std::string sym = write_symfile(
        "p21_super.json",
        R"({ "group_name": "P2_1", "ops": [ "x,y,z", "-x,-y,z+1/2" ], "supercell": [1, 1, 2] })");
    const std::vector<TestOp> ops = {OP_IDENTITY,
                                     {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, {0, 0, 0.5}},
                                     {{{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}}, {0, 0, 0.25}},
                                     {{{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}}, {0, 0, 0.75}}};

    HIDE_OUTPUT([&] {
        command("create_atoms 1 single 1.0 1.5 1.0");
        command("create_atoms 1 single 5.0 4.5 4.0");
        command("create_atoms 1 single 1.0 1.5 7.0");
        command("create_atoms 1 single 5.0 4.5 10.0");
        command("mass 1 1.0");
        command("pair_style lj/cut 2.5");
        command("pair_coeff * * 1.0 1.0");
        command("velocity all create 0.5 4711 mom yes");
        command("fix 1 all nve");
    });
    auto output = CAPTURE_OUTPUT([&] {
        command(fmt::format("fix sym all symmetry {}", sym));
        command("run 1000 post no");
    });
    EXPECT_THAT(output, ContainsRegex("with 4 operators, 4 atoms in 1 orbits"));
    EXPECT_LT(orbit_residual({1, 3, 2, 4}, ops), 1.0e-12);

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   Special position with an orbit of intermediate size: in Pmm2 an atom on
//   the y = 0 mirror plane has one image (from the two-fold axis and the
//   x mirror), so the orbit has 2 members for 4 operators.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, IntermediateWyckoffOrbit)
{
    std::string sym = write_symfile("pmm2.json",
        R"({ "group_name": "Pmm2", "ops": [ "x,y,z", "-x,-y,z", "-x,y,z", "x,-y,z" ] })");

    HIDE_OUTPUT([&] {
        command("units lj");
        command("atom_style atomic");
        command("boundary p p p");
        command("atom_modify map array");
        command("region box block 0 8 0 8 0 8");
        command("create_box 1 box");
        // orbit of size 2 on the y = 0 mirror plane
        command("create_atoms 1 single 1.2 0.0 2.4");
        command("create_atoms 1 single 6.8 0.0 2.4");
        // general position orbit of size 4
        command("create_atoms 1 single 1.0 1.5 3.0");
        command("create_atoms 1 single 7.0 6.5 3.0");
        command("create_atoms 1 single 7.0 1.5 3.0");
        command("create_atoms 1 single 1.0 6.5 3.0");
        command("mass 1 1.0");
        command("pair_style lj/cut 2.5");
        command("pair_coeff * * 1.0 1.0");
        command("pair_modify shift yes");
        command("velocity all create 0.7 12345");
        command("fix 1 all nve");
        command("variable etot equal etotal");
    });
    auto output = CAPTURE_OUTPUT([&] {
        command(fmt::format("fix sym all symmetry {}", sym));
        command("run 0 post no");
    });
    EXPECT_THAT(output, ContainsRegex("6 atoms in 2 orbits .1 on special positions. with 5 "
                                      "independent coordinates"));
    const double etot_start = get_variable_value("etot");
    HIDE_OUTPUT([&] { command("run 1000 post no"); });

    const std::vector<TestOp> ops = {OP_IDENTITY,
                                     {{{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}}, {0, 0, 0}},
                                     {{{-1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, {0, 0, 0}},
                                     {{{1, 0, 0}, {0, -1, 0}, {0, 0, 1}}, {0, 0, 0}}};
    EXPECT_LT(orbit_residual({3, 4, 5, 6}, ops), 1.0e-12);
    EXPECT_LT(orbit_residual({1, 2, 2, 1}, ops), 1.0e-12);
    // atom 1 is free to move in the mirror plane, but not across it
    double x1[3];
    get_pos(1, x1);
    EXPECT_NEAR(x1[1] / 8.0 - std::round(x1[1] / 8.0), 0.0, 1.0e-12);
    EXPECT_NEAR(get_variable_value("etot"), etot_start, 1.0e-3);

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   Special position: an atom on a mirror plane must stay on the mirror
//   plane throughout the run, and an atom that is placed within tol of the
//   mirror is moved onto it.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, WyckoffSiteOnMirrorIsPreserved)
{
    const double L = 8.0;
    std::string sym = write_symfile("pm_wyckoff.json", JSON_PM);

    HIDE_OUTPUT([&] {
        command("units lj");
        command("atom_style atomic");
        command("boundary p p p");
        command("atom_modify map array");
        command(fmt::format("region box block 0 {} 0 {} 0 {}", L, L, L));
        command("create_box 1 box");
        // atom 2 is the mirror image of atom 1 but off by less than tol
        command("create_atoms 1 single 1.2 2.5 3.0");
        command(fmt::format("create_atoms 1 single 1.2 {} 3.0", L - 2.5 + 2.0e-6));
        // atom 3 is on the mirror plane at y = L/2 but off by less than tol
        command(fmt::format("create_atoms 1 single 1.5 {} 2.0", 0.5 * L + 2.0e-6));
        command("mass 1 1.0");
        command("pair_style lj/cut 2.5");
        command("pair_coeff * * 1.0 1.0");
        command("velocity all create 0.7 12345 mom yes rot yes");
        command("fix 1 all nve");
        command(fmt::format("fix sym all symmetry {}", sym));
        command("run 0 post no");
    });

    double x1[3], x2[3], x3[3];
    get_pos(1, x1);
    get_pos(2, x2);
    get_pos(3, x3);
    EXPECT_NEAR(x3[1], 0.5 * L, 1.0e-14);
    EXPECT_NEAR(x1[1] + x2[1], L, 1.0e-14);

    HIDE_OUTPUT([&] { command("run 50 post no"); });

    get_pos(1, x1);
    get_pos(2, x2);
    get_pos(3, x3);
    EXPECT_NEAR(x3[1], 0.5 * L, 1.0e-12) << "atom on the y=L/2 mirror plane moved off it";
    EXPECT_NEAR(x1[1] + x2[1], L, 1.0e-12);
    EXPECT_NEAR(x1[0], x2[0], 1.0e-12);
    EXPECT_NEAR(x1[2], x2[2], 1.0e-12);

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   Orbits found from positions in an FCC lattice: with P-1 the 8 atoms
//   on inversion centers (corners, face and edge centers of the box at
//   multiples of half the box length) are on special positions
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, OrbitsOfFccLattice)
{
    std::string sym = write_symfile("pinv_fcc.json", JSON_PINV);

    auto output = CAPTURE_OUTPUT([&] {
        command("units lj");
        command("lattice fcc 0.8442");
        command("region box block 0 4 0 4 0 4");
        command("create_box 1 box");
        command("create_atoms 1 box");
        command("mass 1 1.0");
        command("pair_style lj/cut 2.5");
        command("pair_coeff * * 1.0 1.0");
        command(fmt::format("fix sym all symmetry {}", sym));
        command("run 0 post no");
    });
    EXPECT_THAT(output, ContainsRegex("256 atoms in 132 orbits .8 on special positions. with 372 "
                                      "independent coordinates"));

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   Only atoms in the fix group are constrained, but the images of all
//   group atoms must be in the group
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, GroupSelectsAtoms)
{
    std::string sym = write_symfile("pinv_group.json", JSON_PINV);

    setup_pinv_system(8.0);
    HIDE_OUTPUT([&] {
        command("create_atoms 1 single 4.6 2.5 6.0");
        command("group sym id 1:8");
        command("velocity all create 0.7 12345 mom yes rot yes");
    });
    auto output = CAPTURE_OUTPUT([&] {
        command(fmt::format("fix sym sym symmetry {}", sym));
        command("run 200 post no");
    });
    EXPECT_THAT(output, ContainsRegex("8 atoms in 4 orbits"));
    // atom 9 moves freely and breaks the symmetry of the forces, so the
    // orbit atoms stay symmetric only because of the fix
    for (tagint a = 1; a <= 4; ++a)
        EXPECT_LT(orbit_residual({a, a + 4}, {OP_IDENTITY, OP_INVERSION}), 1.0e-12);

    HIDE_OUTPUT([&] {
        command("unfix sym");
        command("group part id 1:7");
        command(fmt::format("fix sym part symmetry {}", sym));
    });
    TEST_FAILURE("no atom of group part at the image of atom ID 4 under symmetry operator 2",
                 command("run 0 post no"););

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   Errors for atoms that are not symmetric
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, AsymmetricAtomsTriggerErrors)
{
    std::string sym = write_symfile("pinv_err.json", JSON_PINV);

    setup_pinv_system(8.0);
    HIDE_OUTPUT([&] { command(fmt::format("fix sym all symmetry {}", sym)); });

    // atom 6 moved by more than tol
    HIDE_OUTPUT([&] { command("set atom 6 x 5.1"); });
    TEST_FAILURE("no atom of group all at the image of atom ID 2 under symmetry operator 2",
                 command("run 0 post no"););

    // missing image
    HIDE_OUTPUT([&] { command("set atom 6 x 5.0"); command("delete_atoms group all compress no"); });
    HIDE_OUTPUT([&] {
        command("create_atoms 1 single 1.2 1.5 1.8");
        command("create_atoms 1 single 6.8 6.5 6.2");
        command("create_atoms 1 single 3.0 1.2 0.9");
    });
    TEST_FAILURE("no atom of group all at the image of atom ID 3 under symmetry operator 2",
                 command("run 0 post no"););

    // overlapping atoms
    HIDE_OUTPUT([&] { command("create_atoms 1 single 5.0 6.8 7.1"); command("create_atoms 1 single 5.0 6.8 7.1"); });
    TEST_FAILURE("more than one atom at the image of atom ID 3 under symmetry operator 2",
                 command("run 0 post no"););

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   All atoms of an orbit must have the same type.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, MixedTypesInOrbitTriggersError)
{
    std::string sym = write_symfile("pinv_types.json", JSON_PINV);

    HIDE_OUTPUT([&] {
        command("units lj");
        command("atom_style atomic");
        command("region box block 0 8 0 8 0 8");
        command("create_box 2 box");
        command("create_atoms 1 single 1.2 1.5 1.8");
        command("create_atoms 1 single 3.0 1.2 0.9");
        command("create_atoms 1 single 1.5 3.0 1.1");
        command("create_atoms 1 single 0.9 1.1 3.2");
        command("create_atoms 1 single 6.8 6.5 6.2");
        command("create_atoms 1 single 5.0 6.8 7.1");
        command("create_atoms 2 single 6.5 5.0 6.9");
        command("create_atoms 1 single 7.1 6.9 4.8");
        command("mass * 1.0");
        command("pair_style lj/cut 2.5");
        command("pair_coeff * * 1.0 1.0");
    });

    TEST_FAILURE("atom ID 7 at the image of atom ID 3 under symmetry operator 2 has a different "
                 "atom type: 2 instead of 1",
                 command(fmt::format("fix sym all symmetry {}", sym));
                 command("run 0 post no"););

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   Metric check: a four-fold axis along z requires a box with equal edge
//   lengths along x and y.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, MetricMismatchTriggersError)
{
    std::string sym = write_symfile(
        "p4.json", R"({ "group_name": "P4", "ops": [ "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z" ] })");

    HIDE_OUTPUT([&] {
        command("units lj");
        command("atom_style atomic");
        command("boundary p p p");
        command("region box block 0 6 0 8 0 6");
        command("create_box 1 box");
        command("create_atoms 1 single 0 0 1");
        command("mass 1 1.0");
        command("pair_style lj/cut 2.5");
        command("pair_coeff * * 1.0 1.0");
    });

    TEST_FAILURE("symmetry operator 2 is not compatible with the simulation box: a = 6, b = 8",
                 command(fmt::format("fix sym all symmetry {}", sym));
                 command("run 0 post no"););

    // same operators in a square box are accepted
    HIDE_OUTPUT([&] {
        command("change_box all y final 0 6");
        command(fmt::format("fix sym all symmetry {}", sym));
        command("run 0 post no");
    });

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   Errors and warnings for the content of the symmetry file
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, SymmetryFileErrors)
{
    HIDE_OUTPUT([&] {
        command("region box block 0 8 0 8 0 8");
        command("create_box 1 box");
    });

    struct Case {
        const char *json;
        const char *error;
    };
    const Case cases[] = {
        {"{ not valid json", "Error parsing symmetry file"},
        {R"({ "group_name": "P1" })", "the \"ops\" entry is missing"},
        {R"({ "ops": [ ] })", "\"ops\" must be a non-empty list"},
        {R"({ "ops": [ "-x,-y,-z" ] })", "must include the identity x,y,z"},
        {R"({ "ops": [ "x,y,z", "-x,-y,-z+1", "-x,-y,-z" ] })", "operators 2 and 3 in the symmetry "
                                                                 "file are the same"},
        {R"({ "ops": [ "x,y,z", "-y,x,z" ] })", "the combination of operators 2 and 2 is missing"},
        {R"({ "ops": [ "x,y,z", "x+y,y,z" ] })", "the combination of operators 2 and 2 is missing"},
        {R"({ "ops": [ "x,y,z", "2*x,y,z" ] })", "operator 2 is not a symmetry operator: the "
                                                 "determinant of its matrix is 2"},
        {R"({ "ops": [ "x,y,z", "-x,-y" ] })", "operator '-x,-y' must have 3 comma separated"},
        {R"({ "ops": [ "x,y,z", { "R": [[-1,0,0],[0,-1,0],[0,0,-1]] } ] })",
         "operator 2: must be a string like"},
        {R"({ "ops": [ "x,y,z", { "R": [[-1,0,0],[0,-1.5,0],[0,0,-1]], "t": [0,0,0] } ] })",
         "operator 2: matrix elements must be integers, got -1.5"},
        {R"({ "ops": [ "x,y,z", "-y,x-y,z", "-x+y,-x,z" ], "supercell": [2, 1, 1] })",
         "operator 2 is not compatible with a supercell of 2 x 1 x 1 unit cells"},
        {R"({ "ops": [ "x,y,z" ], "supercell": [0, 1, 1] })", "must be a list of 3 positive"},
    };
    for (const auto &c : cases) {
        std::string sym = write_symfile("error.json", c.json);
        TEST_FAILURE(c.error, command(fmt::format("fix sym all symmetry {}", sym)););
        LAMMPS_NS::platform::unlink(sym);
    }

    // the old format with explicit orbits triggers a warning
    std::string sym = write_symfile("orbits.json", R"({ "ops": [ "x,y,z" ], "orbits": [ ] })");
    auto output = CAPTURE_OUTPUT([&] { command(fmt::format("fix sym all symmetry {}", sym)); });
    EXPECT_THAT(output, ContainsRegex("ignoring unknown entry \"orbits\" in symmetry file"));
    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   RESPA support: with run_style respa, per-level force symmetrization
//   keeps the orbit invariants holding to FP precision throughout the run.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, RespaPreservesSymmetry)
{
    std::string sym = write_symfile("pinv_respa.json", JSON_PINV);

    setup_pinv_system(8.0);
    HIDE_OUTPUT([&] {
        // 2-level RESPA: bond on inner, pair on outer (no bonds in this
        // system but the integrator pattern still exercises per-level
        // force evaluation and the post_force_respa hook).
        command("run_style respa 2 2 bond 1 pair 2");
        command(fmt::format("fix sym all symmetry {}", sym));
        command("run 25 post no");
    });

    for (tagint a = 1; a <= 4; ++a)
        EXPECT_LT(orbit_residual({a, a + 4}, {OP_IDENTITY, OP_INVERSION}), 1.0e-12);

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   The temperature must be computed from the independent degrees of
//   freedom: the free coordinates of the orbit representatives minus the
//   translations of the whole system the group allows.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, TemperatureUsesIndependentDof)
{
    // P-1: 4 general orbits with 3 free coordinates each, no translation
    std::string sym = write_symfile("pinv_dof.json", JSON_PINV);

    setup_pinv_system(8.0);
    HIDE_OUTPUT([&] {
        command("thermo_modify norm no");
        command("variable dof equal 2.0*ke/temp");
        command(fmt::format("fix sym all symmetry {}", sym));
        command("run 0 post no");
    });
    EXPECT_NEAR(get_variable_value("dof"), 12.0, 1.0e-10);
    LAMMPS_NS::platform::unlink(sym);

    // Pm: one general orbit (3), one atom on the mirror (2), minus the two
    // translations parallel to the mirror plane
    sym = write_symfile("pm_dof.json", JSON_PM);
    HIDE_OUTPUT([&] {
        command("clear");
        command("region box block 0 8 0 8 0 8");
        command("create_box 1 box");
        command("create_atoms 1 single 1.2 2.5 3.0");
        command("create_atoms 1 single 1.2 5.5 3.0");
        command("create_atoms 1 single 1.5 4.0 2.0");
        command("mass 1 1.0");
        command("pair_style lj/cut 2.5");
        command("pair_coeff * * 1.0 1.0");
        command("velocity all create 0.7 12345");
        command("thermo_modify norm no");
        command("variable dof equal 2.0*ke/temp");
        command(fmt::format("fix sym all symmetry {}", sym));
        command("run 0 post no");
    });
    EXPECT_NEAR(get_variable_value("dof"), 3.0, 1.0e-10);

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   Minimization with symmetrized forces keeps the symmetry
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, MinimizeKeepsSymmetry)
{
    std::string sym = write_symfile("pinv_min.json", JSON_PINV);

    setup_pinv_system(8.0);
    HIDE_OUTPUT([&] {
        command("pair_modify shift yes");
        command("variable pe equal pe");
        command(fmt::format("fix sym all symmetry {}", sym));
        command("run 0 post no");
    });
    const double pe_start = get_variable_value("pe");
    HIDE_OUTPUT([&] { command("minimize 0.0 1.0e-8 1000 10000"); });
    EXPECT_LT(get_variable_value("pe"), pe_start - 0.1);
    for (tagint a = 1; a <= 4; ++a)
        EXPECT_LT(orbit_residual({a, a + 4}, {OP_IDENTITY, OP_INVERSION}), 1.0e-12);

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   Four-fold axis with a changing box: under anisotropic box relaxation
//   and pressure control the symmetric stress keeps the box square.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, FourFoldAxisWithBoxChange)
{
    std::string sym = write_symfile(
        "p4_box.json", R"({ "group_name": "P4", "ops": [ "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z" ] })");
    const std::vector<TestOp> ops = {OP_IDENTITY,
                                     {{{0, -1, 0}, {1, 0, 0}, {0, 0, 1}}, {0, 0, 0}},
                                     {{{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}}, {0, 0, 0}},
                                     {{{0, 1, 0}, {-1, 0, 0}, {0, 0, 1}}, {0, 0, 0}}};

    HIDE_OUTPUT([&] {
        command("units lj");
        command("atom_style atomic");
        command("atom_modify map array");
        command("region box block 0 4 0 4 0 3");
        command("create_box 1 box");
        // fractional (0.2,0.1,0.5) and (0.5,0.35,0.1) and their C4 images
        command("create_atoms 1 single 0.8 0.4 1.5");
        command("create_atoms 1 single 3.6 0.8 1.5");
        command("create_atoms 1 single 3.2 3.6 1.5");
        command("create_atoms 1 single 0.4 3.2 1.5");
        command("create_atoms 1 single 2.0 1.4 0.3");
        command("create_atoms 1 single 2.6 2.0 0.3");
        command("create_atoms 1 single 2.0 2.6 0.3");
        command("create_atoms 1 single 1.4 2.0 0.3");
        command("mass 1 1.0");
        command("pair_style lj/cut 2.5");
        command("pair_coeff * * 1.0 1.0");
        command("pair_modify shift yes");
        command(fmt::format("fix sym all symmetry {}", sym));
        command("fix relax all box/relax aniso 0.0");
        command("minimize 0.0 1.0e-8 1000 10000");
    });
    EXPECT_LT(orbit_residual({1, 2, 3, 4}, ops), 1.0e-12);
    EXPECT_LT(orbit_residual({5, 6, 7, 8}, ops), 1.0e-12);
    EXPECT_NEAR(lmp->domain->xprd, lmp->domain->yprd, 1.0e-10);

    HIDE_OUTPUT([&] {
        command("unfix relax");
        command("velocity all create 0.3 4711");
        command("run 0 post no");
        command("fix 1 all npt temp 0.3 0.3 0.2 aniso 0.0 0.0 2.0");
        command(fmt::format("fix sym all symmetry {}", sym));
        command("run 500 post no");
    });
    EXPECT_LT(orbit_residual({1, 2, 3, 4}, ops), 1.0e-12);
    EXPECT_LT(orbit_residual({5, 6, 7, 8}, ops), 1.0e-12);
    EXPECT_NEAR(lmp->domain->xprd, lmp->domain->yprd, 1.0e-10);

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   2d systems: in-plane two-fold rotation; operators must keep z unchanged
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, TwoDimensions)
{
    std::string sym = write_symfile("p2_2d.json", R"({ "group_name": "p2", "ops": [ "x,y,z", "-x,-y,z" ] })");
    const std::vector<TestOp> ops = {OP_IDENTITY,
                                     {{{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}}, {0, 0, 0}}};

    HIDE_OUTPUT([&] {
        command("units lj");
        command("dimension 2");
        command("atom_style atomic");
        command("atom_modify map array");
        command("region box block 0 6 0 6 -0.5 0.5");
        command("create_box 1 box");
        command("create_atoms 1 single 1.0 1.5 0.0");
        command("create_atoms 1 single 2.5 4.0 0.0");
        command("create_atoms 1 single 5.0 4.5 0.0");
        command("create_atoms 1 single 3.5 2.0 0.0");
        command("mass 1 1.0");
        command("pair_style lj/cut 2.5");
        command("pair_coeff * * 1.0 1.0");
        command("pair_modify shift yes");
        command("velocity all create 0.5 4711");
        command("thermo_modify norm no");
        command("variable dof equal 2.0*ke/temp");
        command("fix 1 all nve");
        command("fix 2 all enforce2d");
        command(fmt::format("fix sym all symmetry {}", sym));
        command("run 0 post no");
    });
    // 2 general orbits with 2 free coordinates each, no translation
    EXPECT_NEAR(get_variable_value("dof"), 4.0, 1.0e-10);

    HIDE_OUTPUT([&] { command("run 1000 post no"); });
    EXPECT_LT(orbit_residual({1, 3}, ops), 1.0e-12);
    EXPECT_LT(orbit_residual({2, 4}, ops), 1.0e-12);
    LAMMPS_NS::platform::unlink(sym);

    sym = write_symfile("p2_z.json", JSON_PINV);
    TEST_FAILURE("symmetry operator 2 changes the z coordinate",
                 command(fmt::format("fix sym2 all symmetry {}", sym)););
    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   Fixes changing forces must be applied before fix symmetry
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, WarnAboutFixOrder)
{
    std::string sym = write_symfile("pinv_order.json", JSON_PINV);

    setup_pinv_system(8.0);
    auto output = CAPTURE_OUTPUT([&] {
        command(fmt::format("fix sym all symmetry {}", sym));
        command("fix kick all addforce 0.1 0.0 0.0");
        command("run 0 post no");
    });
    EXPECT_THAT(output, ContainsRegex("Fix kick .style addforce. may change forces after fix "
                                      "symmetry sym"));

    output = CAPTURE_OUTPUT([&] {
        command("unfix sym");
        command(fmt::format("fix sym all symmetry {}", sym));
        command("run 0 post no");
    });
    EXPECT_THAT(output, Not(ContainsRegex("may change forces after")));
    LAMMPS_NS::platform::unlink(sym);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
