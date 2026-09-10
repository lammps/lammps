/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   This file is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

// Verify that PPPM can project a source-group potential to a zero-charge
// sensor.  A symmetric finite difference of the reciprocal-space energy with
// respect to the sensor charge provides an independent reference; its only
// extra term is the analytically known neutralizing-background derivative.

#include "../testing/core.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "group.h"
#include "kspace.h"
#include "math_const.h"
#include "pppm_tip4p_xtb.h"
#include "pppm_xtb.h"
#include "gtest/gtest.h"

#include <cmath>
#include <cstdlib>
#include <mpi.h>
#include <string>
#include <vector>

using namespace LAMMPS_NS;
using namespace MathConst;

// Whether the shared test helpers should expose captured LAMMPS output.
bool verbose = false;

namespace {

class KSpaceGroupPotentialTest : public LAMMPSTest, public ::testing::WithParamInterface<bool> {};

#if defined(FFT_SINGLE)
// A larger charge displacement avoids cancellation when reciprocal-space
// energies are accumulated with single-precision FFTs.
constexpr double charge_epsilon      = 1.0e-3;
constexpr double potential_tolerance = 1.0e-5;
#else
constexpr double charge_epsilon      = 1.0e-5;
constexpr double potential_tolerance = 2.0e-7;
#endif

TEST_P(KSpaceGroupPotentialTest, PPPMXTBZeroChargeSensorMatchesEnergyDerivative)
{
    if (!info->has_style("kspace", "pppm/xtb")) GTEST_SKIP();
    if (!info->has_style("pair", "coul/long")) GTEST_SKIP();

    HIDE_OUTPUT([&] {
        command("units lj");
        command("atom_style charge");
        command("atom_modify map array");
        command("boundary p p p");
        if (GetParam())
            command("region box prism 0 12 0 11 0 10 2.0 -1.0 1.5 units box");
        else
            command("region box block 0 12 0 11 0 10 units box");
        command("create_box 1 box");
        command("create_atoms 1 single 2.1 3.2 4.3 units box");
        command("create_atoms 1 single 8.4 7.1 1.8 units box");
        command("create_atoms 1 single 5.7 2.6 8.2 units box");
        command("mass 1 1.0");
        command("set atom 1 charge 1.0");
        command("set atom 2 charge -0.25");
        command("set atom 3 charge 0.0");
        command("group source id 1 2");
        command("group sensor id 3");
        command("pair_style coul/long 4.5");
        command("pair_coeff * *");
        command("kspace_style pppm/xtb 1.0e-9");
        command("run 0 post no");
    });

    const int source_group = lmp->group->find("source");
    const int sensor_group = lmp->group->find("sensor");
    ASSERT_GE(source_group, 0);
    ASSERT_GE(sensor_group, 0);

    const int sensor = lmp->atom->map(3);
    ASSERT_GE(sensor, 0);
    ASSERT_LT(sensor, lmp->atom->nlocal);

    std::vector<double> potential(lmp->atom->nmax, 0.0);
    auto *xtb_kspace = dynamic_cast<PPPMXTB *>(lmp->force->kspace);
    ASSERT_NE(xtb_kspace, nullptr);
    xtb_kspace->compute_group_potential(potential.data(), lmp->group->bitmask[sensor_group],
                                        lmp->group->bitmask[source_group], false);
    EXPECT_NE(potential[sensor], 0.0);

    auto reciprocal_energy = [&](double sensor_charge) {
        lmp->atom->q[sensor] = sensor_charge;
        lmp->comm->forward_comm();
        lmp->force->kspace->qsum_qsq(0);
        lmp->force->kspace->compute(1, 0);
        return lmp->force->kspace->energy;
    };

    const double energy_plus  = reciprocal_energy(charge_epsilon);
    const double energy_minus = reciprocal_energy(-charge_epsilon);
    reciprocal_energy(0.0);

    const double derivative    = (energy_plus - energy_minus) / (2.0 * charge_epsilon);
    const double source_charge = 0.75;
    const double volume        = lmp->domain->xprd * lmp->domain->yprd * lmp->domain->zprd;
    const double alpha         = lmp->force->kspace->g_ewald;
    const double background    = -MY_PI * source_charge / (alpha * alpha * volume);
    const double projected     = lmp->force->qqrd2e * (potential[sensor] + background);

    EXPECT_NEAR(projected, derivative, potential_tolerance);
}

TEST_P(KSpaceGroupPotentialTest, PPPMTIP4PXTBUsesVirtualMChargeSite)
{
    if (!info->has_style("kspace", "pppm/tip4p/xtb")) GTEST_SKIP();
    if (!info->has_style("pair", "lj/cut/tip4p/long")) GTEST_SKIP();
    if (!info->has_style("bond", "harmonic")) GTEST_SKIP();
    if (!info->has_style("angle", "harmonic")) GTEST_SKIP();

    HIDE_OUTPUT([&] {
        command("units lj");
        command("atom_style full");
        command("atom_modify map array");
        command("boundary p p p");
        if (GetParam())
            command("region box prism 0 12 0 11 0 10 2.0 -1.0 1.5 units box");
        else
            command("region box block 0 12 0 11 0 10 units box");
        command("create_box 3 box bond/types 1 angle/types 1 "
                "extra/bond/per/atom 2 extra/angle/per/atom 1 extra/special/per/atom 2");
        command("create_atoms 1 single 3.0 3.0 3.0 units box");
        command("create_atoms 2 single 3.9572 3.0 3.0 units box");
        command("create_atoms 2 single 2.760013 3.926627 3.0 units box");
        command("create_atoms 3 single 8.4 7.1 1.8 units box");
        command("mass * 1.0");
        command("set atom 1 charge -1.1128");
        command("set atom 2*3 charge 0.5564");
        command("set atom 4 charge 0.0");
        command("create_bonds single/bond 1 1 2");
        command("create_bonds single/bond 1 1 3");
        command("create_bonds single/angle 1 2 1 3");
        command("bond_style harmonic");
        command("bond_coeff 1 1000.0 0.9572");
        command("angle_style harmonic");
        command("angle_coeff 1 100.0 104.52");
        command("group source id 1:3");
        command("group sensor id 4");
        command("pair_style lj/cut/tip4p/long 1 2 1 1 0.1546 4.5");
        command("pair_coeff * * 0.0 1.0");
        command("kspace_style pppm/tip4p/xtb 1.0e-9");
        command("run 0 post no");
    });

    const int source_group = lmp->group->find("source");
    const int sensor_group = lmp->group->find("sensor");
    const int sensor       = lmp->atom->map(4);
    ASSERT_GE(source_group, 0);
    ASSERT_GE(sensor_group, 0);
    ASSERT_GE(sensor, 0);
    ASSERT_LT(sensor, lmp->atom->nlocal);

    std::vector<double> potential(lmp->atom->nmax, 0.0);
    auto *xtb_kspace = dynamic_cast<PPPMTIP4PXTB *>(lmp->force->kspace);
    ASSERT_NE(xtb_kspace, nullptr);
    xtb_kspace->compute_group_potential(potential.data(), lmp->group->bitmask[sensor_group],
                                        lmp->group->bitmask[source_group], false);

    auto reciprocal_energy = [&](double sensor_charge) {
        lmp->atom->q[sensor] = sensor_charge;
        lmp->comm->forward_comm();
        lmp->force->kspace->qsum_qsq(0);
        lmp->force->kspace->compute(1, 0);
        return lmp->force->kspace->energy;
    };

    const double derivative =
        (reciprocal_energy(charge_epsilon) - reciprocal_energy(-charge_epsilon)) /
        (2.0 * charge_epsilon);
    reciprocal_energy(0.0);

    EXPECT_NEAR(lmp->force->qqrd2e * potential[sensor], derivative, potential_tolerance);
}

INSTANTIATE_TEST_SUITE_P(OrthogonalAndTriclinic, KSpaceGroupPotentialTest,
                         ::testing::Values(false, true),
                         [](const ::testing::TestParamInfo<bool> &info) {
                             return info.param ? "Triclinic" : "Orthogonal";
                         });

} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);
    const int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
