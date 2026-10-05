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

// Unit tests for the PERI package pair styles (peri/pmb, peri/lps, peri/ves,
// peri/eps) and their accelerator variants.  Unlike the generic force-style
// tester these tests are self-contained: instead of comparing against a stored
// YAML reference they check closed-form physical invariants and use the base
// (non-accelerated) style as the oracle for the accelerator variants.  No YAML
// is processed.

#include "lammps.h"

#include "atom.h"
#include "compute.h"
#include "info.h"
#include "input.h"
#include "modify.h"
#include "utils.h"

#include "fmt/format.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <functional>
#include <mpi.h>
#include <string>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

using LAMMPS_NS::Info;
using LAMMPS_NS::LAMMPS;

namespace LAMMPS_NS {

// per-atom results of a peri run, gathered by atom tag (these tests run on a
// single MPI rank, so nlocal == natoms and tags are dense 1..natoms)
struct PeriData {
    int natoms = 0;
    std::vector<double> z;          // reference-frame z coordinate, by tag-1
    std::vector<int> type;          // atom type, by tag-1
    std::vector<double> damage;     // fraction of broken bonds, by tag-1
    std::vector<double> s0;         // per-atom critical stretch, by tag-1
    std::vector<double> smin;       // per-atom minimum stretch, by tag-1
    std::vector<double> fx, fy, fz; // per-atom force, by tag-1
    double max_damage = 0.0;
};

class PeriTest : public ::testing::Test {
protected:
    LAMMPS *lmp = nullptr;

    void TearDown() override { destroy(); }

    void destroy()
    {
        if (lmp) {
            if (!verbose) ::testing::internal::CaptureStdout();
            delete lmp;
            if (!verbose) ::testing::internal::GetCapturedStdout();
            lmp = nullptr;
        }
    }

    // create a fresh LAMMPS instance with optional extra command-line args
    // (used to request an accelerator suffix, e.g. {"-sf","omp","-pk","omp","1"})
    void create(const std::vector<std::string> &extra = {})
    {
        destroy();
        LAMMPS::argv args = {"PeriTest", "-log", "none", "-nocite"};
        for (const auto &e : extra)
            args.push_back(e);
        if (!verbose) {
            args.push_back("-screen");
            args.push_back("none");
            ::testing::internal::CaptureStdout();
        }
        lmp = new LAMMPS(args, MPI_COMM_WORLD);
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void command(const std::string &line) { lmp->input->one(line); }

    bool has_pair(const std::string &name)
    {
        Info info(lmp);
        return info.has_style("pair", name);
    }

    // Build a two-type peri bar (type 1 for z<4.5, type 2 above), build the
    // bond network, optionally apply a uniform z-stretch about the interface,
    // and run one step so the bond-breaking criterion is exercised.
    //
    // pair_setup must issue the pair_style and pair_coeff commands; it receives
    // the horizon to use.  stretch is the engineering strain applied in z.
    void build_bar(const std::function<void(double)> &pair_setup, double stretch)
    {
        const double horizon = 3.01;
        command("units lj");
        command("boundary s s s");
        command("atom_style peri");
        command("atom_modify map array");
        command("lattice sc 1.0");
        command("region box block 0 3 0 3 0 10 units box");
        command("create_box 2 box");
        command("create_atoms 1 box");
        command("region upper block INF INF INF INF 4.5 INF units box");
        command("set region upper type 2");

        pair_setup(horizon);

        command("set group all density 1.0");
        command("set group all volume 1.0");
        command("compute dmg all damage/atom");
        command("fix 1 all nve");
        command("timestep 1.0e-4");

        // build the bond network at the reference configuration
        command("run 0 post no");

        if (stretch != 0.0) {
            // displacement(z) = stretch*(z-4.5): a uniform z-strain about the
            // type interface so that bonds crossing it are stretched
            double dlo = stretch * (0.0 - 4.5);
            double dhi = stretch * (10.0 - 4.5);
            command(fmt::format("displace_atoms all ramp z {} {} z 0.0 10.0 units box", dlo, dhi));
        }
        command("run 1 post no");
    }

    PeriData extract()
    {
        PeriData d;
        auto *atom = lmp->atom;
        d.natoms   = (int)atom->natoms;
        d.z.assign(d.natoms, 0.0);
        d.type.assign(d.natoms, 0);
        d.damage.assign(d.natoms, 0.0);
        d.s0.assign(d.natoms, 0.0);
        d.smin.assign(d.natoms, 0.0);
        d.fx.assign(d.natoms, 0.0);
        d.fy.assign(d.natoms, 0.0);
        d.fz.assign(d.natoms, 0.0);

        auto *cdmg = lmp->modify->get_compute_by_id("dmg");
        cdmg->compute_peratom();
        double *dmg = cdmg->vector_atom;

        double **x0  = atom->x0;
        double **f   = atom->f;
        int *type    = atom->type;
        double *s0   = atom->s0;
        double *smin = atom->smin;
        tagint *tag  = atom->tag;
        int nlocal   = atom->nlocal;
        for (int i = 0; i < nlocal; ++i) {
            int t = (int)tag[i] - 1;
            if (t < 0 || t >= d.natoms) continue;
            d.z[t]       = x0[i][2];
            d.type[t]    = type[i];
            d.damage[t]  = dmg[i];
            d.s0[t]      = s0[i];
            d.smin[t]    = smin[i];
            d.fx[t]      = f[i][0];
            d.fy[t]      = f[i][1];
            d.fz[t]      = f[i][2];
            d.max_damage = std::max(d.max_damage, dmg[i]);
        }
        return d;
    }

    // pair_coeff helpers for each model: weak s00 (s00_12) on the 1-2 interface
    std::function<void(double)> pmb(double s00_bulk, double s00_12, double alpha)
    {
        return [this,s00_bulk,alpha,s00_12](double h) {
            command("pair_style peri/pmb");
            command(fmt::format("pair_coeff * * 1.0 {} {} {}", h, s00_bulk, alpha));
            command(fmt::format("pair_coeff 1 2 1.0 {} {} {}", h, s00_12, alpha));
        };
    }
    std::function<void(double)> lps(double s00_bulk, double s00_12, double alpha)
    {
        return [this,s00_bulk,alpha,s00_12](double h) {
            command("pair_style peri/lps");
            command(fmt::format("pair_coeff * * 1.0 1.0 {} {} {}", h, s00_bulk, alpha));
            command(fmt::format("pair_coeff 1 2 1.0 1.0 {} {} {}", h, s00_12, alpha));
        };
    }
};

// ------------------------------------------------------------------------
// Issue #984: with a weaker critical stretch (s00) on the 1-2 interface, the
// interface bonds must preferentially break.  The pre-fix code collapsed the
// critical stretch into a per-particle scalar (max over bonds) so the weak
// interface s00 was ignored and nothing broke -> this test fails on that code.
// ------------------------------------------------------------------------
TEST_F(PeriTest, pmb_interface_fracture_984)
{
    const double s00_bulk = 0.1, s00_weak = 0.001, alpha = 0.0, strain = 0.01;

    // weak interface: only the 1-2 bonds should break
    create();
    build_bar(pmb(s00_bulk, s00_weak, alpha), strain);
    PeriData weak = extract();

    // average damage close to the interface vs far away in the bulk
    double near_sum = 0.0, far_sum = 0.0;
    int near_n = 0, far_n = 0;
    for (int t = 0; t < weak.natoms; ++t) {
        double dz = std::fabs(weak.z[t] - 4.5);
        if (dz < 1.5) {
            near_sum += weak.damage[t];
            ++near_n;
        } else if (dz > 3.5) {
            far_sum += weak.damage[t];
            ++far_n;
        }
    }
    ASSERT_GT(near_n, 0);
    ASSERT_GT(far_n, 0);

    // the weak interface fractures (this is the regression assertion) ...
    EXPECT_GT(weak.max_damage, 0.0);
    EXPECT_GT(near_sum / near_n, 0.0);
    // ... while the bulk far from the interface does not
    EXPECT_DOUBLE_EQ(far_sum / far_n, 0.0);

    // control: a uniform s00 (no weak interface) does not break at this strain
    create();
    build_bar(pmb(s00_bulk, s00_bulk, alpha), strain);
    PeriData uniform = extract();
    EXPECT_DOUBLE_EQ(uniform.max_damage, 0.0);
}

TEST_F(PeriTest, lps_interface_fracture_984)
{
    const double s00_bulk = 0.1, s00_weak = 0.001, alpha = 0.0, strain = 0.01;
    create();
    build_bar(lps(s00_bulk, s00_weak, alpha), strain);
    PeriData weak = extract();

    double near_sum = 0.0;
    int near_n      = 0;
    for (int t = 0; t < weak.natoms; ++t)
        if (std::fabs(weak.z[t] - 4.5) < 1.5) {
            near_sum += weak.damage[t];
            ++near_n;
        }
    ASSERT_GT(near_n, 0);
    EXPECT_GT(weak.max_damage, 0.0);
    EXPECT_GT(near_sum / near_n, 0.0);
}

// ------------------------------------------------------------------------
// Below the critical stretch nothing breaks (validates the smin first-step
// initialization: no spurious bond breaking on the first step).
// ------------------------------------------------------------------------
TEST_F(PeriTest, no_break_below_threshold)
{
    create();
    build_bar(pmb(0.1, 0.1, 0.25), 0.01); // strain 0.01 < s00 0.1
    PeriData d = extract();
    EXPECT_DOUBLE_EQ(d.max_damage, 0.0);
}

// ------------------------------------------------------------------------
// The per-atom diagnostic s0 must equal s00 - alpha*smin for a single-material
// system (validates the Option-B display value and the new smin field).
// ------------------------------------------------------------------------
TEST_F(PeriTest, s0_equals_critical_stretch)
{
    const double s00 = 0.1, alpha = 0.25;
    create();
    build_bar(pmb(s00, s00, alpha), 0.005);
    PeriData d = extract();

    int checked = 0;
    for (int t = 0; t < d.natoms; ++t) {
        // skip atoms whose smin is still the (no-bond) sentinel
        if (d.smin[t] > 0.5e308 || d.smin[t] < -0.5e308) continue;
        EXPECT_NEAR(d.s0[t], s00 - alpha * d.smin[t], 1.0e-12);
        ++checked;
    }
    EXPECT_GT(checked, 0);
}

// ------------------------------------------------------------------------
// Total linear momentum is conserved under fix nve with no external forces.
// ------------------------------------------------------------------------
TEST_F(PeriTest, momentum_conservation)
{
    create();
    command("units lj");
    command("boundary s s s");
    command("atom_style peri");
    command("atom_modify map array");
    command("lattice sc 1.0");
    command("region box block 0 4 0 4 0 4 units box");
    command("create_box 1 box");
    command("create_atoms 1 box");
    command("pair_style peri/pmb");
    command("pair_coeff * * 1.0 3.01 0.1 0.25");
    command("set group all density 1.0");
    command("set group all volume 1.0");
    command("velocity all create 0.05 12345 loop geom");
    command("fix 1 all nve");
    command("timestep 1.0e-3");
    command("run 0 post no");

    auto px = [&]() {
        double **v    = lmp->atom->v;
        double *rmass = lmp->atom->rmass;
        int nlocal    = lmp->atom->nlocal;
        double p[3]   = {0.0, 0.0, 0.0};
        for (int i = 0; i < nlocal; ++i) {
            p[0] += rmass[i] * v[i][0];
            p[1] += rmass[i] * v[i][1];
            p[2] += rmass[i] * v[i][2];
        }
        return std::array<double, 3>{p[0], p[1], p[2]};
    };
    auto p0 = px();
    command("run 50 post no");
    auto p1 = px();
    for (int k = 0; k < 3; ++k)
        EXPECT_NEAR(p0[k], p1[k], 1.0e-10);
}

// ------------------------------------------------------------------------
// Accelerator consistency: the /omp variant must reproduce the base style's
// per-atom forces and damage.  Skipped when the OPENMP package is not built.
// ------------------------------------------------------------------------
TEST_F(PeriTest, omp_consistency)
{
    {
        create();
        bool have = has_pair("peri/pmb/omp");
        destroy();
        if (!have) GTEST_SKIP() << "peri/pmb/omp not available";
    }

    const double s00_bulk = 0.1, s00_weak = 0.001, alpha = 0.0, strain = 0.01;

    create();
    build_bar(pmb(s00_bulk, s00_weak, alpha), strain);
    PeriData base = extract();

    create({"-sf", "omp", "-pk", "omp", "1"});
    build_bar(pmb(s00_bulk, s00_weak, alpha), strain);
    PeriData omp = extract();

    ASSERT_EQ(base.natoms, omp.natoms);
    for (int t = 0; t < base.natoms; ++t) {
        EXPECT_NEAR(base.fx[t], omp.fx[t], 1.0e-10);
        EXPECT_NEAR(base.fy[t], omp.fy[t], 1.0e-10);
        EXPECT_NEAR(base.fz[t], omp.fz[t], 1.0e-10);
        EXPECT_DOUBLE_EQ(base.damage[t], omp.damage[t]);
    }
}

// ------------------------------------------------------------------------
// Accelerator consistency for KOKKOS.  No peri/*/kk styles exist yet, so this
// is skipped; it is in place so that a future KOKKOS port is validated against
// the base style automatically.
// ------------------------------------------------------------------------
TEST_F(PeriTest, kokkos_consistency)
{
    create();
    bool have = has_pair("peri/pmb/kk");
    destroy();
    if (!have) GTEST_SKIP() << "peri/pmb/kk not available (no KOKKOS port yet)";

    const double s00_bulk = 0.1, s00_weak = 0.001, alpha = 0.0, strain = 0.01;

    create();
    build_bar(pmb(s00_bulk, s00_weak, alpha), strain);
    PeriData base = extract();

    create({"-k", "on", "-sf", "kk"});
    build_bar(pmb(s00_bulk, s00_weak, alpha), strain);
    PeriData kk = extract();

    ASSERT_EQ(base.natoms, kk.natoms);
    for (int t = 0; t < base.natoms; ++t) {
        EXPECT_NEAR(base.fx[t], kk.fx[t], 1.0e-10);
        EXPECT_NEAR(base.fy[t], kk.fy[t], 1.0e-10);
        EXPECT_NEAR(base.fz[t], kk.fz[t], 1.0e-10);
        EXPECT_DOUBLE_EQ(base.damage[t], kk.damage[t]);
    }
}

} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = LAMMPS_NS::utils::split_words(var);
        for (const auto &arg : env)
            if (arg == "-v") verbose = true;
    }
    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
