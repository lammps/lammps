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
#include "lammps.h"
#include "platform.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <fstream>
#include <string>

bool verbose = false;

using LAMMPS_NS::Atom;
using LAMMPS_NS::tagint;

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
//   Box validator: declared "cubic" lattice with a non-cubic box must error.
// --------------------------------------------------------------------------
TEST_F(FixSymmetryTest, BoxMismatchTriggersError)
{
    std::string sym = write_symfile("cubic.json", R"({
        "group_name": "test",
        "lattice":    "cubic",
        "ops":   [ { "R": [[1,0,0],[0,1,0],[0,0,1]], "t": [0,0,0] } ],
        "orbits": []
    })");

    HIDE_OUTPUT([&] {
        command("units lj");
        command("atom_style atomic");
        command("boundary p p p");
        command("atom_modify map array");
        command("region box block 0 6 0 8 0 6");
        command("create_box 1 box");
        command("create_atoms 1 single 1 1 1");
        command("mass 1 1.0");
        command("pair_style lj/cut 2.5");
        command("pair_coeff * * 1.0 1.0");
    });

    TEST_FAILURE("cubic lattice requires xprd == yprd == zprd",
                 command(fmt::format("fix sym all symmetry {}", sym));
                 command("run 0"););

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

    TEST_FAILURE("differs from R\\*s_asym",
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

    TEST_FAILURE("atom tag 999 .* does not exist",
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

    TEST_FAILURE("declared on Wyckoff site .* but residual",
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

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
