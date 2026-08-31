/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under the
   GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// Verify that global energy and global virial from a KSpace style remain
// unchanged when the ENERGY_ONLY path skips force-producing work.  Also
// verify that requesting per-atom tallies with ENERGY_ONLY retains the full
// path.  Monte Carlo fixes use this interface for trial-move energies.

#include "../testing/core.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "kspace.h"
#include "utils.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <mpi.h>
#include <string>
#include <vector>

bool verbose = false;

namespace LAMMPS_NS {

class KSpaceEnergyOnlyTest : public LAMMPSTest {
protected:
    void build_system(const std::string &kspace_cmd, const std::string &modify_cmd, bool triclinic,
                      bool slab = false)
    {
        HIDE_OUTPUT([&] {
            command("clear");
            command("units real");
            command("atom_style charge");
            if (slab)
                command("boundary p p f");
            else
                command("boundary p p p");
            if (triclinic)
                command("region box prism 0 40 0 40 0 40 5 3 4");
            else
                command("region box block 0 40 0 40 0 40");
            command("create_box 1 box");

            uint64_t state    = 246813579ULL;
            auto next_uniform = [&state]() {
                state = state * 6364136223846793005ULL + 1442695040888963407ULL;
                return static_cast<double>((state >> 11) & ((1ULL << 53) - 1)) /
                       static_cast<double>(1ULL << 53);
            };

            for (int atom_id = 1; atom_id <= 100; ++atom_id) {
                const double x = 1.0 + 38.0 * next_uniform();
                const double y = 1.0 + 38.0 * next_uniform();
                const double z = 1.0 + 38.0 * next_uniform();
                command(fmt::format("create_atoms 1 single {:.15g} {:.15g} {:.15g} units box", x, y,
                                    z));
                const double charge = (atom_id % 2 == 1) ? 1.0 : -1.0;
                command(fmt::format("set atom {} charge {:.1f}", atom_id, charge));
            }

            command("mass 1 1.0");
            command("pair_style coul/long 8.0");
            command("pair_coeff * *");
            command("pair_modify table 0");
            command(kspace_cmd);
            if (!modify_cmd.empty()) command(modify_cmd);
            command("neighbor 2.0 bin");
            command("run 0 post no");
        });
    }

    void build_system_disp(const std::string &kspace_cmd, const std::string &mix_cmd)
    {
        HIDE_OUTPUT([&] {
            command("clear");
            command("units real");
            command("atom_style charge");
            command("boundary p p p");
            command("region box block 0 40 0 40 0 40");
            command("create_box 2 box");

            uint64_t state    = 246813579ULL;
            auto next_uniform = [&state]() {
                state = state * 6364136223846793005ULL + 1442695040888963407ULL;
                return static_cast<double>((state >> 11) & ((1ULL << 53) - 1)) /
                       static_cast<double>(1ULL << 53);
            };

            for (int atom_id = 1; atom_id <= 100; ++atom_id) {
                const double x = 1.0 + 38.0 * next_uniform();
                const double y = 1.0 + 38.0 * next_uniform();
                const double z = 1.0 + 38.0 * next_uniform();
                const int type = (atom_id % 2 == 1) ? 1 : 2;
                command(fmt::format("create_atoms {} single {:.15g} {:.15g} {:.15g} units box",
                                    type, x, y, z));
                const double charge = (atom_id % 2 == 1) ? 1.0 : -1.0;
                command(fmt::format("set atom {} charge {:.1f}", atom_id, charge));
            }

            command("mass * 1.0");
            command("pair_style lj/long/coul/long long long 8.0");
            command("pair_coeff 1 1 0.2 2.5");
            command("pair_coeff 2 2 0.1 3.0");
            if (!mix_cmd.empty()) command(mix_cmd);
            command(kspace_cmd);
            command("kspace_modify disp/auto yes");
            command("neighbor 2.0 bin");
            command("run 0 post no");
        });
    }

    void compute_and_snapshot(int eflag, int vflag, double &e, double v[6])
    {
        KSpace *ks = lmp->force->kspace;
        ks->compute(eflag, vflag);
        e = ks->energy;
        for (int i = 0; i < 6; ++i)
            v[i] = ks->virial[i];
    }

    std::vector<std::array<double, 3>> snapshot_forces()
    {
        std::vector<std::array<double, 3>> result(lmp->atom->nlocal);
        for (int i = 0; i < lmp->atom->nlocal; ++i)
            result[i] = {lmp->atom->f[i][0], lmp->atom->f[i][1], lmp->atom->f[i][2]};
        return result;
    }

    static double relerr(double a, double b)
    {
        const double den = std::max({std::fabs(a), std::fabs(b), 1.0});
        return std::fabs(a - b) / den;
    }

    void check_energy_only(const std::string &label)
    {
        double e_full, v_full[6];
        double e_only, v_only[6];

        compute_and_snapshot(ENERGY_GLOBAL, VIRIAL_PAIR, e_full, v_full);
        const auto forces_after_full = snapshot_forces();
        compute_and_snapshot(ENERGY_GLOBAL | ENERGY_ONLY, VIRIAL_PAIR, e_only, v_only);

        EXPECT_LT(relerr(e_full, e_only), 1.0e-12) << label << " energy";
        for (int i = 0; i < 6; ++i)
            EXPECT_LT(relerr(v_full[i], v_only[i]), 1.0e-12) << label << " virial[" << i << "]";

        ASSERT_EQ(forces_after_full.size(), static_cast<size_t>(lmp->atom->nlocal));
        for (int i = 0; i < lmp->atom->nlocal; ++i) {
            EXPECT_DOUBLE_EQ(forces_after_full[i][0], lmp->atom->f[i][0])
                << label << " fx[" << i << "]";
            EXPECT_DOUBLE_EQ(forces_after_full[i][1], lmp->atom->f[i][1])
                << label << " fy[" << i << "]";
            EXPECT_DOUBLE_EQ(forces_after_full[i][2], lmp->atom->f[i][2])
                << label << " fz[" << i << "]";
        }

        // Per-atom requests must retain the full path even when ENERGY_ONLY
        // is also present, because they require the force-producing data.
        double e_pa, v_pa[6];
        compute_and_snapshot(ENERGY_GLOBAL | ENERGY_ATOM | ENERGY_ONLY, VIRIAL_PAIR | VIRIAL_ATOM,
                             e_pa, v_pa);
        EXPECT_LT(relerr(e_full, e_pa), 1.0e-12) << label << " per-atom-fallback energy";
        for (int i = 0; i < 6; ++i)
            EXPECT_LT(relerr(v_full[i], v_pa[i]), 1.0e-12)
                << label << " per-atom-fallback virial[" << i << "]";

        if (verbose && (lmp->comm->me == 0))
            utils::print("{:<16} e_full={:.12g} e_only={:.12g} rel.err E={:.2e}\n", label, e_full,
                         e_only, relerr(e_full, e_only));
    }
};

TEST_F(KSpaceEnergyOnlyTest, Ewald)
{
    if (!info->has_style("kspace", "ewald")) GTEST_SKIP();
    if (!info->has_style("pair", "coul/long")) GTEST_SKIP();

    build_system("kspace_style ewald 1.0e-4", "", false);
    check_energy_only("ewald");

    build_system("kspace_style ewald 1.0e-4", "", true);
    check_energy_only("ewald/triclinic");
}

TEST_F(KSpaceEnergyOnlyTest, PPPMik)
{
    if (!info->has_style("kspace", "pppm")) GTEST_SKIP();
    if (!info->has_style("pair", "coul/long")) GTEST_SKIP();

    build_system("kspace_style pppm 1.0e-4", "", false);
    check_energy_only("pppm/ik");

    build_system("kspace_style pppm 1.0e-4", "", true);
    check_energy_only("pppm/ik/triclinic");
}

TEST_F(KSpaceEnergyOnlyTest, EwaldSlab)
{
    if (!info->has_style("kspace", "ewald")) GTEST_SKIP();
    if (!info->has_style("pair", "coul/long")) GTEST_SKIP();

    build_system("kspace_style ewald 1.0e-4", "kspace_modify slab 3.0", false, true);
    check_energy_only("ewald/slab");
}

TEST_F(KSpaceEnergyOnlyTest, PPPMSlab)
{
    if (!info->has_style("kspace", "pppm")) GTEST_SKIP();
    if (!info->has_style("pair", "coul/long")) GTEST_SKIP();

    build_system("kspace_style pppm 1.0e-4", "kspace_modify slab 3.0", false, true);
    check_energy_only("pppm/slab");
}

TEST_F(KSpaceEnergyOnlyTest, PPPMad)
{
    if (!info->has_style("kspace", "pppm")) GTEST_SKIP();
    if (!info->has_style("pair", "coul/long")) GTEST_SKIP();

    build_system("kspace_style pppm 1.0e-4", "kspace_modify diff ad", false);
    check_energy_only("pppm/ad");
}

TEST_F(KSpaceEnergyOnlyTest, PPPMcg)
{
    if (!info->has_style("kspace", "pppm/cg")) GTEST_SKIP();
    if (!info->has_style("pair", "coul/long")) GTEST_SKIP();

    build_system("kspace_style pppm/cg 1.0e-4", "", false);
    check_energy_only("pppm/cg");

    build_system("kspace_style pppm/cg 1.0e-4", "", true);
    check_energy_only("pppm/cg/triclinic");
}

TEST_F(KSpaceEnergyOnlyTest, PPPMstagger)
{
    if (!info->has_style("kspace", "pppm/stagger")) GTEST_SKIP();
    if (!info->has_style("pair", "coul/long")) GTEST_SKIP();

    build_system("kspace_style pppm/stagger 1.0e-4", "", false);
    check_energy_only("pppm/stagger");
}

TEST_F(KSpaceEnergyOnlyTest, PPPMdisp)
{
    if (!info->has_style("kspace", "pppm/disp")) GTEST_SKIP();
    if (!info->has_style("pair", "lj/long/coul/long")) GTEST_SKIP();

    build_system_disp("kspace_style pppm/disp 1.0e-4", "pair_modify mix geometric");
    check_energy_only("pppm/disp/geom");

    build_system_disp("kspace_style pppm/disp 1.0e-4", "pair_modify mix arithmetic");
    check_energy_only("pppm/disp/arith");

    build_system_disp("kspace_style pppm/disp 1.0e-4", "pair_coeff 1 2 0.15 2.7");
    check_energy_only("pppm/disp/none");
}

TEST_F(KSpaceEnergyOnlyTest, EwaldDisp)
{
    if (!info->has_style("kspace", "ewald/disp")) GTEST_SKIP();
    if (!info->has_style("pair", "lj/long/coul/long")) GTEST_SKIP();

    build_system_disp("kspace_style ewald/disp 1.0e-4", "pair_modify mix geometric");
    check_energy_only("ewald/disp");
}

} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = LAMMPS_NS::utils::split_words(var);
        for (const auto &arg : env) {
            if (arg == "-v") verbose = true;
        }
    }

    if ((argc > 1) && (std::string(argv[1]) == "-v")) verbose = true;

    const int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
