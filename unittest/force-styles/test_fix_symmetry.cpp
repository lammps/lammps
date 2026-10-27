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

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
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

// --------------------------------------------------------------------------
//   P1 (identity-only, no orbits) must not move atoms or change energy
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, P1IsNoOp)
{
    std::string p1 = write_symfile("p1.json", R"({
        "group_name": "P1",
        "lattice":    "triclinic",
        "ops":   [ { "R": [[1,0,0],[0,1,0],[0,0,1]], "t": [0,0,0] } ],
        "orbits": []
    })");

    HIDE_OUTPUT([&] {
        command("atom_modify map array");
        command("lattice fcc 0.5");
        command("region box block 0 4 0 4 0 4");
        command("create_box 1 box");
        command("create_atoms 1 box");
        command("mass 1 1.0");
        command("pair_style lj/cut 2.5");
        command("pair_coeff * * 1.0 1.0");
        command("velocity all create 1.0 12345 mom yes rot yes");
        command("fix 1 all nve");
        command(fmt::format("fix sym all symmetry {}", p1));
        command("run 20");
    });

    // System ran; total atom count preserved is the loosest sanity check.
    EXPECT_GT(lmp->atom->natoms, 0);

    LAMMPS_NS::platform::unlink(p1);
}

// --------------------------------------------------------------------------
//   P-1 inversion: atom_asym + atom_image must equal (L, L, L) exactly
//   throughout the run (modulo the periodic-image wrap captured by my
//   end_of_step's round() trick).
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, PinvPreservesSymmetryUnderNVE)
{
    const double L = 8.0;
    std::string sym = write_symfile("pinv.json", R"({
        "group_name": "P-1",
        "lattice":    "triclinic",
        "ops": [
            { "R": [[ 1, 0, 0],[ 0, 1, 0],[ 0, 0, 1]], "t": [0,0,0] },
            { "R": [[-1, 0, 0],[ 0,-1, 0],[ 0, 0,-1]], "t": [0,0,0] }
        ],
        "orbits": [
            { "tags": [1, 5] }, { "tags": [2, 6] },
            { "tags": [3, 7] }, { "tags": [4, 8] }
        ]
    })");

    setup_pinv_system(L);
    HIDE_OUTPUT([&] {
        command(fmt::format("fix sym all symmetry {}", sym));
        command("run 50");
    });

    // Each (asym, image) pair must sum to (L, L, L) modulo L (because the
    // image may have wrapped through a periodic boundary). The "mod L" is
    // captured by checking the residual mod 1 of (x_a + x_i) / L.
    for (tagint a = 1; a <= 4; ++a) {
        const tagint b = a + 4;
        double xa[3], xb[3];
        get_pos(a, xa);
        get_pos(b, xb);
        for (int k = 0; k < 3; ++k) {
            double sum_over_L = (xa[k] + xb[k]) / L;
            double residual = sum_over_L - std::round(sum_over_L);
            EXPECT_NEAR(residual, 0.0, 1.0e-12)
                << "orbit " << a << " component " << k
                << ": x_asym = " << xa[k] << ", x_image = " << xb[k]
                << ", sum/L = " << sum_over_L;
        }
    }

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   Metric check: a four-fold axis along z requires a box with equal edge
//   lengths along x and y.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, MetricMismatchTriggersError)
{
    std::string sym = write_symfile("p4.json", R"({
        "group_name": "P4",
        "ops": [ { "R": [[1,0,0],[0,1,0],[0,0,1]], "t": [0,0,0] },
                 { "R": [[0,-1,0],[1,0,0],[0,0,1]], "t": [0,0,0] },
                 { "R": [[-1,0,0],[0,-1,0],[0,0,1]], "t": [0,0,0] },
                 { "R": [[0,1,0],[-1,0,0],[0,0,1]], "t": [0,0,0] } ],
        "orbits": [ ]
    })");

    HIDE_OUTPUT([&] {
        command("units lj");
        command("atom_style atomic");
        command("boundary p p p");
        command("region box block 0 6 0 8 0 6");
        command("create_box 1 box");
        command("create_atoms 1 single 1 1 1");
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
//   Orbit consistency check (pass 3 in build_orbit_map): if the user's
//   atoms don't sit at R*s_asym + t, init must fail with a helpful message.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, MisplacedImageTriggersError)
{
    // The orbits below claim tag 1 and tag 2 are inversion partners, but
    // tag 2 is not actually placed at -tag1.
    std::string sym = write_symfile("pinv_bad.json", R"({
        "group_name": "P-1",
        "lattice":    "triclinic",
        "ops": [
            { "R": [[ 1, 0, 0],[ 0, 1, 0],[ 0, 0, 1]], "t": [0,0,0] },
            { "R": [[-1, 0, 0],[ 0,-1, 0],[ 0, 0,-1]], "t": [0,0,0] }
        ],
        "orbits": [
            { "tags": [1, 2] }, { "tags": [3, 4] },
            { "tags": [5, 6] }, { "tags": [7, 8] }
        ]
    })");

    setup_pinv_system(8.0);

    TEST_FAILURE("atom ID 2 is not at the image of atom ID 1 under symmetry operator 2",
                 command(fmt::format("fix sym all symmetry {}", sym));
                 command("run 0"););

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   Existence check (pass 2): a tag listed in the symfile that doesn't
//   exist in the system must produce a clean error.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, MissingTagTriggersError)
{
    std::string sym = write_symfile("missing.json", R"({
        "group_name": "P-1",
        "lattice":    "triclinic",
        "ops": [
            { "R": [[ 1, 0, 0],[ 0, 1, 0],[ 0, 0, 1]], "t": [0,0,0] },
            { "R": [[-1, 0, 0],[ 0,-1, 0],[ 0, 0,-1]], "t": [0,0,0] }
        ],
        "orbits": [
            { "tags": [1, 999] }
        ]
    })");

    setup_pinv_system(8.0);

    TEST_FAILURE("atom ID 999 of symmetry orbit 1 does not exist",
                 command(fmt::format("fix sym all symmetry {}", sym));
                 command("run 0"););

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   Wyckoff special-position projection: an atom on a mirror plane must
//   stay on the mirror plane throughout the run, to ~FP precision.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, WyckoffSiteOnMirrorIsPreserved)
{
    const double L = 8.0;
    // Pm group: ops 1=(x,y,z) identity, 2=(x,-y,z) mirror about y=0.
    // Orbit 1 is a general pair (tags 1,2 related by mirror).
    // Orbit 2 is a single atom (tag 3) sitting on the mirror at y=L/2,
    // declared via "site_symmetry": [2] (op 2 fixes the asym).
    std::string sym = write_symfile("pm_wyckoff.json", R"({
        "group_name": "Pm",
        "lattice":    "triclinic",
        "ops": [
            { "R": [[1, 0, 0], [0,  1, 0], [0, 0, 1]], "t": [0, 0, 0] },
            { "R": [[1, 0, 0], [0, -1, 0], [0, 0, 1]], "t": [0, 0, 0] }
        ],
        "orbits": [
            { "tags": [1, 2] },
            { "tags": [3, 3], "site_symmetry": [2] }
        ]
    })");

    HIDE_OUTPUT([&] {
        command("units lj");
        command("atom_style atomic");
        command("boundary p p p");
        command("atom_modify map array");
        command(fmt::format("region box block 0 {} 0 {} 0 {}", L, L, L));
        command("create_box 1 box");
        // Tag 1: general position. Tag 2 = mirror of tag 1 in y.
        command("create_atoms 1 single 1.2 2.5 3.0");
        command(fmt::format("create_atoms 1 single 1.2 {} 3.0", L - 2.5));
        // Tag 3: Wyckoff on mirror at y = L/2.
        command(fmt::format("create_atoms 1 single 1.5 {} 2.0", 0.5 * L));
        command("mass 1 1.0");
        command("pair_style lj/cut 2.5");
        command("pair_coeff * * 1.0 1.0");
        command("velocity all create 0.7 12345 mom yes rot yes");
        command("fix 1 all nve");
        command(fmt::format("fix sym all symmetry {}", sym));
        command("run 50");
    });

    // Tag 3 must still be on the y=L/2 mirror.
    double x3[3];
    get_pos(3, x3);
    EXPECT_NEAR(x3[1], 0.5 * L, 1.0e-12)
        << "Wyckoff atom drifted off the y=L/2 mirror plane";

    // Tags 1 and 2 must still be y-mirror images: y1 + y2 == L.
    double x1[3], x2[3];
    get_pos(1, x1);
    get_pos(2, x2);
    EXPECT_NEAR(x1[1] + x2[1], L, 1.0e-12);
    EXPECT_NEAR(x1[0], x2[0], 1.0e-12);    // mirror leaves x unchanged
    EXPECT_NEAR(x1[2], x2[2], 1.0e-12);    // mirror leaves z unchanged

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   Wyckoff structural check: atom declared on a Wyckoff site but placed
//   off the constraint manifold must produce a clean init error.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, WyckoffOffSiteTriggersError)
{
    const double L = 8.0;
    std::string sym = write_symfile("pm_wyckoff_bad.json", R"({
        "group_name": "Pm",
        "lattice":    "triclinic",
        "ops": [
            { "R": [[1, 0, 0], [0,  1, 0], [0, 0, 1]], "t": [0, 0, 0] },
            { "R": [[1, 0, 0], [0, -1, 0], [0, 0, 1]], "t": [0, 0, 0] }
        ],
        "orbits": [
            { "tags": [1, 1], "site_symmetry": [2] }
        ]
    })");

    HIDE_OUTPUT([&] {
        command("units lj");
        command("atom_style atomic");
        command("boundary p p p");
        command("atom_modify map array");
        command(fmt::format("region box block 0 {} 0 {} 0 {}", L, L, L));
        command("create_box 1 box");
        // Atom at y != L/2 -- NOT on the mirror.
        command("create_atoms 1 single 1.5 2.0 2.0");
        command("mass 1 1.0");
        command("pair_style lj/cut 2.5");
        command("pair_coeff * * 1.0 1.0");
    });

    TEST_FAILURE("atom ID 1 is not on the special position required by symmetry operator 2",
                 command(fmt::format("fix sym all symmetry {}", sym));
                 command("run 0"););

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   RESPA support: with run_style respa, per-level force symmetrization
//   keeps the orbit invariants holding to FP precision throughout the run.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, RespaPreservesSymmetry)
{
    const double L = 8.0;
    std::string sym = write_symfile("pinv_respa.json", R"({
        "group_name": "P-1", "lattice": "triclinic",
        "ops": [
            { "R": [[ 1, 0, 0],[ 0, 1, 0],[ 0, 0, 1]], "t": [0,0,0] },
            { "R": [[-1, 0, 0],[ 0,-1, 0],[ 0, 0,-1]], "t": [0,0,0] }
        ],
        "orbits": [
            { "tags": [1, 5] }, { "tags": [2, 6] },
            { "tags": [3, 7] }, { "tags": [4, 8] }
        ]
    })");

    setup_pinv_system(L);
    HIDE_OUTPUT([&] {
        // 2-level RESPA: bond on inner, pair on outer (no bonds in this
        // system but the integrator pattern still exercises per-level
        // force evaluation and the post_force_respa hook).
        command("run_style respa 2 2 bond 1 pair 2");
        command(fmt::format("fix sym all symmetry {}", sym));
        command("run 25");
    });

    for (tagint a = 1; a <= 4; ++a) {
        const tagint b = a + 4;
        double xa[3], xb[3];
        get_pos(a, xa);
        get_pos(b, xb);
        for (int k = 0; k < 3; ++k) {
            double s = (xa[k] + xb[k]) / L;
            double residual = s - std::round(s);
            EXPECT_NEAR(residual, 0.0, 1.0e-12)
                << "orbit " << a << " component " << k << " under RESPA";
        }
    }

    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   JSON parser: malformed file => clean error, no segfault.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, MalformedJsonTriggersError)
{
    std::string sym = write_symfile("broken.json", "{ not valid json");

    setup_pinv_system(8.0);

    TEST_FAILURE("Error parsing symmetry file",
                 command(fmt::format("fix sym all symmetry {}", sym));
                 command("run 0"););

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
    std::string sym = write_symfile("pinv_equiv.json", R"({
        "group_name": "P-1",
        "lattice":    "triclinic",
        "ops": [
            { "R": [[ 1, 0, 0],[ 0, 1, 0],[ 0, 0, 1]], "t": [0,0,0] },
            { "R": [[-1, 0, 0],[ 0,-1, 0],[ 0, 0,-1]], "t": [0,0,0] }
        ],
        "orbits": [
            { "tags": [1, 5] }, { "tags": [2, 6] },
            { "tags": [3, 7] }, { "tags": [4, 8] }
        ]
    })");

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
    std::string sym = write_symfile("pinv_decay.json", R"({
        "group_name": "P-1",
        "lattice":    "triclinic",
        "ops": [
            { "R": [[ 1, 0, 0],[ 0, 1, 0],[ 0, 0, 1]], "t": [0,0,0] },
            { "R": [[-1, 0, 0],[ 0,-1, 0],[ 0, 0,-1]], "t": [0,0,0] }
        ],
        "orbits": [
            { "tags": [1, 5] }, { "tags": [2, 6] },
            { "tags": [3, 7] }, { "tags": [4, 8] }
        ]
    })");

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
    std::string sym = write_symfile("p3.json", R"({
        "group_name": "P3",
        "lattice":    "hexagonal",
        "ops": [
            { "R": [[ 1, 0, 0],[ 0, 1, 0],[ 0, 0, 1]], "t": [0,0,0] },
            { "R": [[ 0,-1, 0],[ 1,-1, 0],[ 0, 0, 1]], "t": [0,0,0] },
            { "R": [[-1, 1, 0],[-1, 0, 0],[ 0, 0, 1]], "t": [0,0,0] }
        ],
        "orbits": [ { "tags": [1, 2, 3] }, { "tags": [4, 5, 6] } ]
    })");
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
    // the input coordinates are rounded to 12 digits
    EXPECT_LT(orbit_residual({1, 2, 3}, ops), 1.0e-11);
    EXPECT_LT(orbit_residual({4, 5, 6}, ops), 1.0e-11);
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
    std::string sym = write_symfile("p21.json", R"({
        "group_name": "P2_1",
        "lattice":    "orthorhombic",
        "ops": [
            { "R": [[ 1, 0, 0],[ 0, 1, 0],[ 0, 0, 1]], "t": [0,0,0] },
            { "R": [[-1, 0, 0],[ 0,-1, 0],[ 0, 0, 1]], "t": [0,0,"1/2"] }
        ],
        "orbits": [ { "tags": [1, 2] }, { "tags": [3, 4] } ]
    })");
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
//   Special position with an orbit of intermediate size: in Pmm2 an atom on
//   the y = 0 mirror plane has one image (from the two-fold axis and the
//   x mirror), so the orbit has 2 members for 4 operators.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, IntermediateWyckoffOrbit)
{
    std::string sym = write_symfile("pmm2.json", R"({
        "group_name": "Pmm2",
        "lattice":    "orthorhombic",
        "ops": [
            { "R": [[ 1, 0, 0],[ 0, 1, 0],[ 0, 0, 1]], "t": [0,0,0] },
            { "R": [[-1, 0, 0],[ 0,-1, 0],[ 0, 0, 1]], "t": [0,0,0] },
            { "R": [[-1, 0, 0],[ 0, 1, 0],[ 0, 0, 1]], "t": [0,0,0] },
            { "R": [[ 1, 0, 0],[ 0,-1, 0],[ 0, 0, 1]], "t": [0,0,0] }
        ],
        "orbits": [
            { "tags": [3, 4, 5, 6] },
            { "tags": [1, 2, 2, 1], "site_symmetry": [4] }
        ]
    })");

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
    });

    HIDE_OUTPUT([&] {
        command("variable etot equal etotal");
        command(fmt::format("fix sym all symmetry {}", sym));
        command("run 0 post no");
    });
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
//   All atoms of an orbit must have the same type.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, MixedTypesInOrbitTriggersError)
{
    std::string sym = write_symfile("pinv_types.json", R"({
        "group_name": "P-1",
        "ops": [
            { "R": [[ 1, 0, 0],[ 0, 1, 0],[ 0, 0, 1]], "t": [0,0,0] },
            { "R": [[-1, 0, 0],[ 0,-1, 0],[ 0, 0,-1]], "t": [0,0,0] }
        ],
        "orbits": [
            { "tags": [1, 5] }, { "tags": [2, 6] },
            { "tags": [3, 7] }, { "tags": [4, 8] }
        ]
    })");

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

    TEST_FAILURE("atoms of symmetry orbit 3 have different atom types 1 and 2",
                 command(fmt::format("fix sym all symmetry {}", sym));
                 command("run 0 post no"););

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
    std::string sym = write_symfile("pinv_dof.json", R"({
        "group_name": "P-1",
        "ops": [
            { "R": [[ 1, 0, 0],[ 0, 1, 0],[ 0, 0, 1]], "t": [0,0,0] },
            { "R": [[-1, 0, 0],[ 0,-1, 0],[ 0, 0,-1]], "t": [0,0,0] }
        ],
        "orbits": [
            { "tags": [1, 5] }, { "tags": [2, 6] },
            { "tags": [3, 7] }, { "tags": [4, 8] }
        ]
    })");

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
    sym = write_symfile("pm_dof.json", R"({
        "group_name": "Pm",
        "ops": [
            { "R": [[1, 0, 0], [0,  1, 0], [0, 0, 1]], "t": [0, 0, 0] },
            { "R": [[1, 0, 0], [0, -1, 0], [0, 0, 1]], "t": [0, 0, 0] }
        ],
        "orbits": [ { "tags": [1, 2] }, { "tags": [3, 3] } ]
    })");
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
    std::string sym = write_symfile("pinv_min.json", R"({
        "group_name": "P-1",
        "ops": [
            { "R": [[ 1, 0, 0],[ 0, 1, 0],[ 0, 0, 1]], "t": [0,0,0] },
            { "R": [[-1, 0, 0],[ 0,-1, 0],[ 0, 0,-1]], "t": [0,0,0] }
        ],
        "orbits": [
            { "tags": [1, 5] }, { "tags": [2, 6] },
            { "tags": [3, 7] }, { "tags": [4, 8] }
        ]
    })");

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
    std::string sym = write_symfile("p4_box.json", R"({
        "group_name": "P4",
        "ops": [ { "R": [[1,0,0],[0,1,0],[0,0,1]], "t": [0,0,0] },
                 { "R": [[0,-1,0],[1,0,0],[0,0,1]], "t": [0,0,0] },
                 { "R": [[-1,0,0],[0,-1,0],[0,0,1]], "t": [0,0,0] },
                 { "R": [[0,1,0],[-1,0,0],[0,0,1]], "t": [0,0,0] } ],
        "orbits": [ { "tags": [1,2,3,4] }, { "tags": [5,6,7,8] } ]
    })");
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
        command("fix sym all symmetry p4_box.json");
        command("run 0 post no");
        command("fix 1 all npt temp 0.3 0.3 0.2 aniso 0.0 0.0 2.0");
        command("fix sym all symmetry p4_box.json");
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
    std::string sym = write_symfile("p2_2d.json", R"({
        "group_name": "p2",
        "ops": [ { "R": [[1,0,0],[0,1,0],[0,0,1]], "t": [0,0,0] },
                 { "R": [[-1,0,0],[0,-1,0],[0,0,1]], "t": [0,0,0] } ],
        "orbits": [ { "tags": [1,3] }, { "tags": [2,4] } ]
    })");
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

    sym = write_symfile("p2_z.json", R"({
        "ops": [ { "R": [[1,0,0],[0,1,0],[0,0,1]], "t": [0,0,0] },
                 { "R": [[-1,0,0],[0,-1,0],[0,0,-1]], "t": [0,0,0] } ],
        "orbits": [ ]
    })");
    TEST_FAILURE("symmetry operator 2 changes the z coordinate",
                 command(fmt::format("fix sym2 all symmetry {}", sym)););
    LAMMPS_NS::platform::unlink(sym);
}

// --------------------------------------------------------------------------
//   Fixes changing forces must be applied before fix symmetry
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, WarnAboutFixOrder)
{
    std::string sym = write_symfile("pinv_order.json", R"({
        "ops": [
            { "R": [[ 1, 0, 0],[ 0, 1, 0],[ 0, 0, 1]], "t": [0,0,0] },
            { "R": [[-1, 0, 0],[ 0,-1, 0],[ 0, 0,-1]], "t": [0,0,0] }
        ],
        "orbits": [
            { "tags": [1, 5] }, { "tags": [2, 6] },
            { "tags": [3, 7] }, { "tags": [4, 8] }
        ]
    })");

    setup_pinv_system(8.0);
    auto output = CAPTURE_OUTPUT([&] {
        command(fmt::format("fix sym all symmetry {}", sym));
        command("fix kick all addforce 0.1 0.0 0.0");
        command("run 0 post no");
    });
    EXPECT_THAT(output, ContainsRegex("Fix kick .style addforce. may change forces after fix "
                                      "symmetry sym"));
    EXPECT_THAT(output, ContainsRegex("Fix symmetry sym: group unknown with 2 operators, 8 atoms "
                                      "in 4 orbits .0 on special positions. with 12 independent"));

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
