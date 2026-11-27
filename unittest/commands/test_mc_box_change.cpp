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

// unit tests for MC pre_exchange fixes (atom/swap, mol/swap, neighbor/swap,
// gcmc, widom, charge/regulation) with changing boxes: expanding shrink-wrap
// boundaries without lost atoms, NPT neighbor bin vs nsq consistency, and
// fix atom/swap defined before a flipping fix npt on a triclinic shrink-wrapped
// box.  can be run with any number of MPI ranks.

#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "exceptions.h"
#include "fix.h"
#include "info.h"
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

Lines lj_fcc(int n, const std::string &atom_style = "atomic", double rho = 0.8442)
{
    const std::string box = fmt::format("region box block 0 {0} 0 {0} 0 {0}", n);
    Lines lines           = {"units lj",
                             "atom_style " + atom_style,
                             "atom_modify map array",
                             fmt::format("lattice fcc {}", rho),
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
    return lines;
}

Lines expanding_shrink_wrap_base(const std::string &atom_style = "atomic")
{
    Lines lines = {"units lj",
                   "atom_style " + atom_style,
                   "atom_modify map array",
                   "boundary s s s",
                   "region box block 0 10 0 10 0 10",
                   "create_box 2 box",
                   "mass * 1.0"};
    int tag     = 1;
    for (int ix = 0; ix < 3; ++ix) {
        for (int iy = 0; iy < 3; ++iy) {
            for (int iz = 0; iz < 3; ++iz, ++tag) {
                const int type = ((ix + iy + iz) % 2 == 0) ? 1 : 2;
                const double x = 1.0 + 4.0 * ix;
                const double y = 1.0 + 4.0 * iy;
                const double z = 1.0 + 4.0 * iz;
                lines.push_back(fmt::format("create_atoms {} single {} {} {}", type, x, y, z));
                const double dx = x - 5.0;
                const double dy = y - 5.0;
                const double dz = z - 5.0;
                const double vx = 0.05 * dx - 0.15 * dy + 0.04 * dz;
                const double vy = 0.15 * dx + 0.05 * dy - 0.03 * dz;
                const double vz = -0.04 * dx + 0.03 * dy + 0.05 * dz;
                lines.push_back(fmt::format("group g{} id {}", tag, tag));
                lines.push_back(fmt::format("velocity g{} set {} {} {}", tag, vx, vy, vz));
                lines.push_back(fmt::format("group g{} delete", tag));
            }
        }
    }
    lines.push_back("pair_style lj/cut 2.5");
    lines.push_back("pair_coeff * * 0.1 1.0 2.5");
    lines.push_back("neighbor 0.3 bin");
    lines.push_back("neigh_modify every 1 delay 0 check no");
    lines.push_back("timestep 0.1");
    lines.push_back("fix 1 all nve");
    return lines;
}
} // namespace

class MCBoxChangeTest : public ::testing::Test {
protected:
    LAMMPS *lmp = nullptr;

    void SetUp() override { lmp = create(MPI_COMM_WORLD); }

    void TearDown() override { destroy(lmp); }

    static LAMMPS *create(MPI_Comm comm)
    {
        LAMMPS::argv args = {"MCBoxChange", "-log", "none", "-echo", "screen", "-nocite"};
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

    int global_nlocal() const
    {
        int nlocal = lmp->atom->nlocal, count = 0;
        MPI_Allreduce(&nlocal, &count, 1, MPI_INT, MPI_SUM, lmp->world);
        return count;
    }

    double box_volume() const
    {
        return lmp->domain->xprd * lmp->domain->yprd * lmp->domain->zprd;
    }

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
        MPI_Allreduce(local.data(), global.data(), (int)local.size(), MPI_DOUBLE, MPI_SUM,
                      instance->world);
        return global;
    }

    std::vector<double> coords_by_tag() const { return coords_by_tag(lmp); }

    void check_expanding_shrink_wrap(const std::string &atom_style, const Lines &extra_setup,
                                     const Lines &mc_cmds)
    {
        reset();
        run(expanding_shrink_wrap_base(atom_style));
        run(extra_setup);
        ASSERT_EQ(lmp->atom->natoms, 27);
        const double vol0 = box_volume();
        run(mc_cmds);
        run({"run 20 post no"});
        ASSERT_EQ(lmp->atom->natoms, 27);
        ASSERT_EQ(global_nlocal(), 27);
        ASSERT_EQ(lmp->update->ntimestep, 20);
        EXPECT_GT(box_volume(), vol0);
    }
};

// expanding shrink-wrapped box must not lose atoms

TEST_F(MCBoxChangeTest, ExpandingShrinkWrapNoLostAtoms)
{
    check_expanding_shrink_wrap("atomic", {}, {"fix mc all atom/swap 1 0 12345 1.5 types 1 2"});
    check_expanding_shrink_wrap("atomic", {},
                                {"fix mc all gcmc 1 0 0 1 12345 1.5 -10.0 0.1"});
    check_expanding_shrink_wrap("atomic", {}, {"fix mc all widom 1 5 1 12345 1.5"});
    EXPECT_NEAR(fix_vector("mc", 2), box_volume(), 1.0e-10 * box_volume());
    check_expanding_shrink_wrap(
        "charge", {"set type 1 charge 1.0", "set type 2 charge -1.0"},
        {"fix mc all charge/regulation 1 2 pIp 3.0 pIm 3.0 onlysalt yes 1 -1 nevery 1 nmc 0 "
         "seed 12345 temp 1.5"});
}

// dynamic triclinic NPT: neighbor bin must match neighbor nsq for atom/swap,
// gcmc (full_energy), and widom (full_energy)

TEST_F(MCBoxChangeTest, DynamicTriclinicNPTBinMatchesNsq)
{
    auto run_triclinic = [&](const std::string &neigh_style, const Lines &mc_cmds, int metric_mode,
                             double rho = 0.8442) {
        reset();
        run(lj_fcc(3, "atomic", rho));
        run({"change_box all triclinic",
             "neighbor 0.3 " + neigh_style,
             "neigh_modify every 1 delay 0 check no",
             "timestep 0.005",
             "fix 1 all npt temp 1.5 1.5 0.5 tri 5.0 5.0 1.0"});
        run(mc_cmds);
        run({"run 25 post no"});
        double metric = 0.0;
        if (metric_mode == 1) {
            metric = fix_vector("mc", 1) + fix_vector("mc", 3) + fix_vector("mc", 5);
        } else {
            metric = fix_vector("mc", 1);
        }
        const double pe = lmp->modify->get_compute_by_id("thermo_pe")->compute_scalar();
        return std::make_pair(metric, pe);
    };

    // atom/swap
    {
        const Lines mc = {"fix mc all atom/swap 5 40 12345 1.5 types 1 2"};
        const auto bin = run_triclinic("bin", mc, 0);
        const auto nsq = run_triclinic("nsq", mc, 0);
        EXPECT_GT(bin.first, 0.0);
        EXPECT_DOUBLE_EQ(bin.first, nsq.first);
        EXPECT_NEAR(bin.second, nsq.second, 1.0e-10 * (1.0 + std::fabs(nsq.second)));
    }

    // gcmc (atomic + full_energy)
    {
        const Lines mc = {"fix mc all gcmc 5 10 20 1 12345 1.5 -2.0 0.2 full_energy"};
        const auto bin = run_triclinic("bin", mc, 1);
        const auto nsq = run_triclinic("nsq", mc, 1);
        EXPECT_GT(bin.first, 0.0);
        EXPECT_DOUBLE_EQ(bin.first, nsq.first);
        EXPECT_NEAR(bin.second, nsq.second, 1.0e-10 * (1.0 + std::fabs(nsq.second)));
    }

    // widom (atomic + full_energy): compare Boltzmann factor vector(1)
    {
        const Lines mc = {"fix mc all widom 5 40 1 12345 1.5 full_energy"};
        const auto bin = run_triclinic("bin", mc, 2, 0.35);
        const auto nsq = run_triclinic("nsq", mc, 2, 0.35);
        EXPECT_GT(bin.first, 0.0);
        EXPECT_NEAR(bin.first, nsq.first, 1.0e-10 * (1.0 + std::fabs(nsq.first)));
        EXPECT_NEAR(bin.second, nsq.second, 1.0e-10 * (1.0 + std::fabs(nsq.second)));
    }
}

// orthogonal NPT: neighbor bin must match neighbor nsq for charge/regulation
// (charge/regulation rejects triclinic boxes in init())

TEST_F(MCBoxChangeTest, ChargeRegulationOrthogonalNPTBinMatchesNsq)
{
    auto run_case = [&](const std::string &neigh_style) {
        reset();
        run(lj_fcc(3, "charge", 0.35));
        run({"set type 1 charge 1.0",
             "set type 2 charge -1.0",
             "neighbor 0.3 " + neigh_style,
             "neigh_modify every 1 delay 0 check no",
             "timestep 0.005",
             "fix 1 all npt temp 1.5 1.5 0.5 iso 2.0 2.0 1.0",
             "fix mc all charge/regulation 1 2 pIp 2.0 pIm 2.0 onlysalt yes 1 -1 nevery 5 "
             "nmc 30 seed 12345 temp 1.5",
             "run 25 post no"});
        const double accepted = fix_vector("mc", 1);
        const double pe       = lmp->modify->get_compute_by_id("thermo_pe")->compute_scalar();
        return std::make_pair(accepted, pe);
    };

    const auto bin = run_case("bin");
    const auto nsq = run_case("nsq");
    EXPECT_GT(bin.first, 0.0);
    EXPECT_DOUBLE_EQ(bin.first, nsq.first);
    EXPECT_NEAR(bin.second, nsq.second, 1.0e-10 * (1.0 + std::fabs(nsq.second)));
}

// fix mol/swap requires PKG_MOLECULE (atom_style molecular)

TEST_F(MCBoxChangeTest, MolSwapChangingBoxes)
{
    Info info(lmp);
    if (!info.has_style("atom", "molecular")) GTEST_SKIP();

    check_expanding_shrink_wrap("molecular", {"set atom * mol 1"},
                                {"fix mc all mol/swap 1 0 1 2 12345 1.5"});

    auto run_case = [&](const std::string &neigh_style) {
        reset();
        run(lj_fcc(3, "molecular"));
        run({"variable molid atom id",
             "set atom * mol v_molid",
             "change_box all triclinic",
             "neighbor 0.3 " + neigh_style,
             "neigh_modify every 1 delay 0 check no",
             "timestep 0.005",
             "fix 1 all npt temp 1.5 1.5 0.5 tri 5.0 5.0 1.0",
             "fix mc all mol/swap 5 40 1 2 12345 1.5",
             "run 25 post no"});
        const double accepted = fix_vector("mc", 1);
        const double pe       = lmp->modify->get_compute_by_id("thermo_pe")->compute_scalar();
        return std::make_pair(accepted, pe);
    };

    const auto bin = run_case("bin");
    const auto nsq = run_case("nsq");
    EXPECT_GT(bin.first, 0.0);
    EXPECT_DOUBLE_EQ(bin.first, nsq.first);
    EXPECT_NEAR(bin.second, nsq.second, 1.0e-10 * (1.0 + std::fabs(nsq.second)));
}

// fix neighbor/swap requires PKG_VORONOI (compute voronoi/atom)

TEST_F(MCBoxChangeTest, NeighborSwapChangingBoxes)
{
    Info info(lmp);
    if (!info.has_style("compute", "voronoi/atom")) GTEST_SKIP();

    check_expanding_shrink_wrap("atomic", {"compute voro all voronoi/atom neighbors yes"},
                                {"fix mc all neighbor/swap 1 0 12345 1.5 2.0 voro types 1 2"});

    auto run_case = [&](const std::string &neigh_style) {
        reset();
        run(lj_fcc(3));
        run({"change_box all triclinic",
             "neighbor 0.3 " + neigh_style,
             "neigh_modify every 1 delay 0 check no",
             "timestep 0.005",
             "compute voro all voronoi/atom neighbors yes",
             "fix 1 all npt temp 1.5 1.5 0.5 tri 5.0 5.0 1.0",
             "fix mc all neighbor/swap 5 40 12345 1.5 2.0 voro types 2 1",
             "run 25 post no"});
        const double accepted = fix_vector("mc", 1);
        const double pe       = lmp->modify->get_compute_by_id("thermo_pe")->compute_scalar();
        return std::make_pair(accepted, pe);
    };

    const auto bin = run_case("bin");
    const auto nsq = run_case("nsq");
    EXPECT_GT(bin.first, 0.0);
    EXPECT_DOUBLE_EQ(bin.first, nsq.first);
    EXPECT_NEAR(bin.second, nsq.second, 1.0e-10 * (1.0 + std::fabs(nsq.second)));
}

// fix atom/swap defined before a flipping fix npt on a triclinic shrink-wrapped
// box must preserve the tilt flip and match plain Verlet on 1 and N MPI ranks

TEST_F(MCBoxChangeTest, TriclinicShrinkWrapFlipBeforeNPT)
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
    ASSERT_EQ(global_nlocal(), 324);
    EXPECT_LT(lmp->domain->xy, 0.0);
    EXPECT_LT(lmp->domain->xz, 0.0);
    const double boxlo[3] = {lmp->domain->boxlo[0], lmp->domain->boxlo[1], lmp->domain->boxlo[2]};
    const double boxhi[3] = {lmp->domain->boxhi[0], lmp->domain->boxhi[1], lmp->domain->boxhi[2]};
    const double tilt[3]  = {lmp->domain->xy, lmp->domain->xz, lmp->domain->yz};
    const double pe       = lmp->modify->get_compute_by_id("thermo_pe")->compute_scalar();
    const auto coords     = coords_by_tag();

    reset();
    run(setup);
    ASSERT_GT(lmp->domain->xy, 0.0);
    ASSERT_GT(lmp->domain->xz, 0.0);
    run({"fix mc all atom/swap 1 0 12345 1.0 types 1 2", npt_cmd, "run 20 post no"});
    ASSERT_EQ(lmp->atom->natoms, 324);
    ASSERT_EQ(global_nlocal(), 324);
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
    const auto actual = coords_by_tag();
    ASSERT_EQ(actual.size(), coords.size());
    for (std::size_t i = 0; i < coords.size(); ++i)
        EXPECT_NEAR(actual[i], coords[i], 1.0e-10 * (1.0 + std::fabs(coords[i])));
}

} // namespace LAMMPS_NS
