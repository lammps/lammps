/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// unit tests for fix atom/swap with changing boxes: dynamic triclinic NPT
// (bin vs nsq), expanding shrink-wrap with cross-rank migration, triclinic
// shrink-wrap idempotence with ncycles 0, and fix atom/swap defined before a
// flipping fix npt on a triclinic shrink-wrapped box.  can be run with any
// number of MPI ranks.

#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "exceptions.h"
#include "fix.h"
#include "input.h"
#include "lammps.h"
#include "modify.h"
#include "update.h"

#include "gtest/gtest.h"

#include "../testing/test_mpi_main.h"

#include <cmath>
#include <string>
#include <utility>
#include <vector>

namespace LAMMPS_NS {

namespace {
using Lines = std::vector<std::string>;

// binary fcc lattice with lj/cut and unequal per-type cutoffs
Lines lj_fcc(int n)
{
    const std::string box = fmt::format("region box block 0 {0} 0 {0} 0 {0}", n);
    return {"units lj",
            "atom_style atomic",
            "atom_modify map array",
            "lattice fcc 0.8442",
            box,
            "create_box 2 box",
            "create_atoms 1 box",
            "mass * 1.0",
            "set type 1 type/fraction 2 0.5 4321",
            "velocity all create 1.5 87287 loop geom",
            "pair_style lj/cut 2.5",
            "pair_coeff 1 1 1.0 1.0 2.5",
            "pair_coeff 2 2 0.8 1.1 2.0",
            "pair_coeff 1 2 0.9 1.05 2.8",
            "neighbor 0.3 bin"};
}
} // namespace

class FixAtomSwapTest : public ::testing::Test {
protected:
    LAMMPS *lmp = nullptr;

    void SetUp() override { lmp = create(MPI_COMM_WORLD); }

    void TearDown() override { destroy(lmp); }

    static LAMMPS *create(MPI_Comm comm)
    {
        LAMMPS::argv args = {"FixAtomSwap", "-log", "none", "-echo", "screen", "-nocite"};
        if (!verbose) ::testing::internal::CaptureStdout();
        auto *instance = new LAMMPS(args, comm);
        if (!verbose) ::testing::internal::GetCapturedStdout();
        return instance;
    }

    static void destroy(LAMMPS *&instance)
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        delete instance;
        instance = nullptr;
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void reset()
    {
        destroy(lmp);
        lmp = create(MPI_COMM_WORLD);
    }

    static std::string run(LAMMPS *instance, const Lines &lines)
    {
        ::testing::internal::CaptureStdout();
        try {
            for (const auto &line : lines)
                instance->input->one(line);
        } catch (LAMMPSException &) {
            auto out = ::testing::internal::GetCapturedStdout();
            utils::print(stdout, "{}", out);
            throw;
        }
        auto out = ::testing::internal::GetCapturedStdout();
        if (verbose) utils::print(stdout, "{}", out);
        return out;
    }

    std::string run(const Lines &lines) { return run(lmp, lines); }

    static double fix_vector(LAMMPS *instance, const std::string &id, int n)
    {
        auto *fix = instance->modify->get_fix_by_id(id);
        if (!fix) return -1.0;
        return fix->compute_vector(n);
    }

    double fix_vector(const std::string &id, int n) { return fix_vector(lmp, id, n); }

    // Cartesian coordinates indexed by global atom tag (3 * (natoms+1))
    static std::vector<double> coords_by_tag(LAMMPS *instance)
    {
        const auto natoms = (std::size_t)instance->atom->natoms;
        std::vector<double> local(3 * (natoms + 1), 0.0), global(3 * (natoms + 1), 0.0);
        for (int i = 0; i < instance->atom->nlocal; ++i) {
            const auto t     = (std::size_t)instance->atom->tag[i];
            local[3 * t]     = instance->atom->x[i][0];
            local[3 * t + 1] = instance->atom->x[i][1];
            local[3 * t + 2] = instance->atom->x[i][2];
        }
        MPI_Allreduce(local.data(), global.data(), (int)(3 * (natoms + 1)), MPI_DOUBLE, MPI_SUM,
                      instance->world);
        return global;
    }
};

// dynamic triclinic NPT: neighbor bin and nsq must agree on accepts and PE

TEST_F(FixAtomSwapTest, DynamicTriclinicNPTBinMatchesNSQ)
{
    auto run_style = [this](const std::string &neigh) {
        reset();
        Lines input;
        for (const auto &line : lj_fcc(4)) {
            if (utils::strmatch(line, "^region box"))
                input.emplace_back("region box prism 0 4 0 4 0 4 0.6 0.4 0.3");
            else if (utils::strmatch(line, "^neighbor"))
                input.emplace_back("neighbor 0.1 " + neigh);
            else
                input.push_back(line);
        }
        input.emplace_back("fix 1 all npt temp 1.5 1.5 0.5 tri 2.0 2.0 2.0");
        input.emplace_back("fix mc all atom/swap 5 10 29494 1.0 ke yes types 1 2");
        input.emplace_back("run 20 post no");
        run(input);
        const double accepts = fix_vector("mc", 1);
        const double pe      = lmp->modify->get_compute_by_id("thermo_pe")->compute_scalar();
        return std::make_pair(accepts, pe);
    };

    const auto bin = run_style("bin");
    const auto nsq = run_style("nsq");
    EXPECT_GT(bin.first, 0.0);
    EXPECT_EQ(bin.first, nsq.first);
    EXPECT_NEAR(bin.second, nsq.second, 1.0e-10 * (1.0 + std::fabs(nsq.second)));
}

// expanding orthogonal shrink-wrap with rotating velocities: all 27 atoms stay
// present while migrating across MPI subdomains, with no lost atoms

TEST_F(FixAtomSwapTest, ExpandingShrinkWrapNoLostAtoms)
{
    const Lines input = {"units lj",
                         "atom_style atomic",
                         "atom_modify map array",
                         "boundary s s s",
                         "lattice sc 1.0 origin 0.5 0.5 0.5",
                         "region box block 0 3 0 3 0 3",
                         "create_box 2 box",
                         "create_atoms 1 box",
                         "set type 1 type/fraction 2 0.5 4321",
                         "mass * 1.0",
                         "variable vx atom 2.0*(x-1.5)+4.0*(y-1.5)",
                         "variable vy atom 2.0*(y-1.5)-4.0*(x-1.5)",
                         "variable vz atom 2.0*(z-1.5)",
                         "velocity all set v_vx v_vy v_vz units box",
                         "pair_style lj/cut 2.5",
                         "pair_coeff * * 0.0 1.0",
                         "neighbor 0.3 bin",
                         "neigh_modify every 1 delay 0 check no",
                         "timestep 0.01",
                         "fix integrate all nve",
                         "fix mc all atom/swap 1 5 12345 1.0 types 1 2",
                         "run 0 post no"};
    run(input);
    ASSERT_EQ(lmp->atom->natoms, 27);
    const double volume_start = lmp->domain->xprd * lmp->domain->yprd * lmp->domain->zprd;
    run({"run 10 post no"});
    EXPECT_EQ(lmp->atom->natoms, 27);
    EXPECT_EQ(lmp->update->ntimestep, 10);
    int nlocal = lmp->atom->nlocal, count = 0;
    MPI_Allreduce(&nlocal, &count, 1, MPI_INT, MPI_SUM, lmp->world);
    EXPECT_EQ(count, 27);
    EXPECT_GT(lmp->domain->xprd * lmp->domain->yprd * lmp->domain->zprd, volume_start);
}

// triclinic shrink-wrap with ncycles 0 must match plain Verlet (expanding)

TEST_F(FixAtomSwapTest, TriclinicShrinkWrapIdempotentZeroCycles)
{
    const Lines setup = {"units lj",
                         "atom_style atomic",
                         "boundary s s p",
                         "atom_modify map array sort 0 0.0",
                         "region box prism -5 5 -5 5 -5 5 2.0 0.0 0.0 units box",
                         "create_box 2 box",
                         "create_atoms 1 single -3 -3 0 units box",
                         "create_atoms 2 single 3 3 0 units box",
                         "create_atoms 1 single -1 -1 1 units box",
                         "create_atoms 2 single 1 1 -1 units box",
                         "mass * 1",
                         "pair_style lj/cut 2.5",
                         "pair_coeff * * 0.0 1.0",
                         "neighbor 0.3 bin",
                         "neigh_modify every 1 delay 0 check no",
                         "group left id 1 3",
                         "group right id 2 4",
                         "group top id 2 3",
                         "group bottom id 1 4",
                         "velocity left set -2.0 0.0 0.0",
                         "velocity right set 2.0 0.0 0.0",
                         "velocity top set NULL 1.0 NULL",
                         "velocity bottom set NULL -1.0 NULL",
                         "timestep 0.0001",
                         "fix integrate all nve"};
    run(setup);
    run({"run 2 post no"});
    const double volume_start = lmp->domain->xprd * lmp->domain->yprd * lmp->domain->zprd;
    run({"run 8 post no"});
    ASSERT_EQ(lmp->atom->natoms, 4);
    ASSERT_EQ(lmp->update->ntimestep, 10);
    const double boxlo[3] = {lmp->domain->boxlo[0], lmp->domain->boxlo[1], lmp->domain->boxlo[2]};
    const double boxhi[3] = {lmp->domain->boxhi[0], lmp->domain->boxhi[1], lmp->domain->boxhi[2]};
    const double xy       = lmp->domain->xy;
    const double volume   = lmp->domain->xprd * lmp->domain->yprd * lmp->domain->zprd;
    const auto coords     = coords_by_tag(lmp);
    EXPECT_GT(volume, volume_start);

    reset();
    run(setup);
    run({"run 2 post no", "fix mc all atom/swap 1 0 12345 1.5 types 1 2", "run 8 post no"});
    ASSERT_EQ(lmp->atom->natoms, 4);
    ASSERT_EQ(lmp->update->ntimestep, 10);
    for (int d = 0; d < 3; ++d) {
        EXPECT_NEAR(lmp->domain->boxlo[d], boxlo[d], 1.0e-10 * (1.0 + std::fabs(boxlo[d])));
        EXPECT_NEAR(lmp->domain->boxhi[d], boxhi[d], 1.0e-10 * (1.0 + std::fabs(boxhi[d])));
    }
    EXPECT_NEAR(lmp->domain->xy, xy, 1.0e-10 * (1.0 + std::fabs(xy)));
    EXPECT_NEAR(lmp->domain->xprd * lmp->domain->yprd * lmp->domain->zprd, volume,
                1.0e-10 * volume);
    const auto actual = coords_by_tag(lmp);
    ASSERT_EQ(actual.size(), coords.size());
    for (std::size_t i = 0; i < coords.size(); ++i)
        EXPECT_NEAR(actual[i], coords[i], 1.0e-10 * (1.0 + std::fabs(coords[i])));
}

// fix atom/swap defined before a flipping fix npt on a triclinic shrink-wrapped
// box must preserve the tilt flip and match plain Verlet on 1 and N MPI ranks

TEST_F(FixAtomSwapTest, TriclinicShrinkWrapFlipBeforeNPT)
{
    const Lines setup = {"units lj",
                         "atom_style atomic",
                         "atom_modify map array sort 0 0.0",
                         "boundary p s s",
                         "lattice fcc 0.8442",
                         "region box prism 0 4 0 4 0 4 1.98 1.98 0.0",
                         "create_box 2 box",
                         "create_atoms 1 box",
                         "set type 1 type/fraction 2 0.5 4321",
                         "mass * 1.0",
                         "pair_style lj/cut 2.5",
                         "pair_coeff * * 1.0 1.0 2.5",
                         "velocity all create 1.5 87287 loop geom",
                         "neighbor 0.3 bin",
                         "neigh_modify every 1 delay 0 check no",
                         "timestep 0.005"};
    const std::string npt_cmd = "fix 1 all npt temp 1.5 1.5 0.5 x 20.0 20.0 0.2 flip yes";

    run(setup);
    ASSERT_EQ(lmp->atom->natoms, 324);
    ASSERT_GT(lmp->domain->xy, 0.0);
    ASSERT_GT(lmp->domain->xz, 0.0);
    run({npt_cmd, "run 20 post no"});
    ASSERT_EQ(lmp->atom->natoms, 324);
    int nlocal = lmp->atom->nlocal, count = 0;
    MPI_Allreduce(&nlocal, &count, 1, MPI_INT, MPI_SUM, lmp->world);
    ASSERT_EQ(count, 324);
    EXPECT_LT(lmp->domain->xy, 0.0);
    EXPECT_LT(lmp->domain->xz, 0.0);
    const double boxlo[3] = {lmp->domain->boxlo[0], lmp->domain->boxlo[1], lmp->domain->boxlo[2]};
    const double boxhi[3] = {lmp->domain->boxhi[0], lmp->domain->boxhi[1], lmp->domain->boxhi[2]};
    const double tilt[3]  = {lmp->domain->xy, lmp->domain->xz, lmp->domain->yz};
    const double pe       = lmp->modify->get_compute_by_id("thermo_pe")->compute_scalar();
    const auto coords     = coords_by_tag(lmp);

    reset();
    run(setup);
    ASSERT_GT(lmp->domain->xy, 0.0);
    ASSERT_GT(lmp->domain->xz, 0.0);
    run({"fix mc all atom/swap 1 0 12345 1.0 types 1 2", npt_cmd, "run 20 post no"});
    ASSERT_EQ(lmp->atom->natoms, 324);
    nlocal = lmp->atom->nlocal;
    count  = 0;
    MPI_Allreduce(&nlocal, &count, 1, MPI_INT, MPI_SUM, lmp->world);
    ASSERT_EQ(count, 324);
    EXPECT_LT(lmp->domain->xy, 0.0);
    EXPECT_LT(lmp->domain->xz, 0.0);
    for (int d = 0; d < 3; ++d) {
        EXPECT_NEAR(lmp->domain->boxlo[d], boxlo[d], 1.0e-10 * (1.0 + std::fabs(boxlo[d])));
        EXPECT_NEAR(lmp->domain->boxhi[d], boxhi[d], 1.0e-10 * (1.0 + std::fabs(boxhi[d])));
    }
    EXPECT_NEAR(lmp->domain->xy, tilt[0], 1.0e-10 * (1.0 + std::fabs(tilt[0])));
    EXPECT_NEAR(lmp->domain->xz, tilt[1], 1.0e-10 * (1.0 + std::fabs(tilt[1])));
    EXPECT_NEAR(lmp->domain->yz, tilt[2], 1.0e-10 * (1.0 + std::fabs(tilt[2])));
    const double actual_pe = lmp->modify->get_compute_by_id("thermo_pe")->compute_scalar();
    EXPECT_NEAR(actual_pe, pe, 1.0e-10 * (1.0 + std::fabs(pe)));
    const auto actual = coords_by_tag(lmp);
    ASSERT_EQ(actual.size(), coords.size());
    for (std::size_t i = 0; i < coords.size(); ++i)
        EXPECT_NEAR(actual[i], coords[i], 1.0e-10 * (1.0 + std::fabs(coords[i])));
}

} // namespace LAMMPS_NS
