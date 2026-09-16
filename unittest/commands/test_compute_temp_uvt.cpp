/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute.h"
#include "fix.h"
#include "modify.h"
#include "../testing/core.h"
#include "gtest/gtest.h"

#define protected public
#include "../../src/EXTRA-FIX/fix_uvt.h"
#undef protected

#include <mpi.h>

bool verbose = false;

namespace LAMMPS_NS {

class ComputeTempUVTTest : public LAMMPSTest {
 protected:
  void SetUp() override
  {
    LAMMPSTest::SetUp();
    if (!Info::has_package("EXTRA-FIX")) GTEST_SKIP();
    command("units lj");
    command("atom_style atomic");
    command("atom_modify map yes");
    command("region box block 0 4 0 4 0 4");
    command("create_box 1 box");
    command("create_atoms 1 single 1 1 1");
    command("create_atoms 1 single 3 3 3");
    command("mass 1 1.0");
    command("pair_style zero 0.5");
    command("pair_coeff * *");
    command("velocity all set 1 2 3");
    command("variable deriv equal 0.0");
    command("fix cp all uvt temp 1 1 0.5 mu 0 0 0.5 ne 1 ne_velocity 2 dedn v_deriv");
    command("compute nuclear all temp");
    command("compute combined all temp/uvt cp");
    command("run 0 post no");
  }
};

TEST_F(ComputeTempUVTTest, CombinedScalarAndNuclearTensor)
{
  auto *temp = lmp->modify->get_compute_by_id("combined");
  auto *nuclear = lmp->modify->get_compute_by_id("nuclear");
  auto *fix = lmp->modify->get_fix_by_id("cp");
  int dim;
  auto *ne_dot = static_cast<double *>(fix->extract("ne_dot", dim));
  auto *ne_mass = static_cast<double *>(fix->extract("ne_mass", dim));
  ASSERT_NE(ne_dot, nullptr);
  ASSERT_NE(ne_mass, nullptr);
  EXPECT_DOUBLE_EQ(*ne_mass, 0.75);
  EXPECT_DOUBLE_EQ(temp->dof, 4.0);
  // Two unit-mass atoms each have |v|^2 = 14; Ne contributes 0.75 * 2^2.
  EXPECT_DOUBLE_EQ(temp->compute_scalar(), 31.0 / 4.0);
  EXPECT_DOUBLE_EQ(temp->compute_scalar(), 31.0 / 4.0);
  *ne_dot = 4.0;
  EXPECT_DOUBLE_EQ(temp->compute_scalar(), 40.0 / 4.0);

  temp->compute_vector();
  nuclear->compute_vector();
  for (int i = 0; i < 6; ++i) EXPECT_DOUBLE_EQ(temp->vector[i], nuclear->vector[i]);
}

TEST_F(ComputeTempUVTTest, FixVectorReportsCombinedTemperature)
{
  auto *fix = lmp->modify->get_fix_by_id("cp");
  ASSERT_EQ(fix->size_vector, 19);
  EXPECT_DOUBLE_EQ(fix->compute_vector(0), 1.0);
  EXPECT_DOUBLE_EQ(fix->compute_vector(5), 1.5);
  EXPECT_DOUBLE_EQ(fix->compute_vector(6), 0.0);
  EXPECT_DOUBLE_EQ(fix->compute_vector(1), 31.0 / 4.0);
  EXPECT_EQ(fix->get_thermo_colname(1), "f_cp:T_ins");
  command("velocity all set 2 0 0");
  EXPECT_DOUBLE_EQ(fix->compute_vector(1), 11.0 / 4.0);
  command("thermo_style custom step f_cp[1] temp f_cp[2]");
  command("thermo_modify norm yes");
  EXPECT_NO_THROW(command("run 0 post no"));
  command("run 5 post no");
  EXPECT_NEAR(fix->compute_vector(1),
              lmp->modify->get_compute_by_id("cp_temp")->compute_scalar(), 1.0e-12);
}

TEST_F(ComputeTempUVTTest, PhysicalOutputIndicesAreIndependentOfChainLength)
{
  for (const int chain : {1, 3, 5}) {
    command("unfix cp");
    command("fix cp all uvt temp 1 1 0.5 mu 0 0 0.5 ne 1 ne_velocity 2 "
            "dedn v_deriv tchain " + std::to_string(chain));
    command("run 0 post no");
    auto *fix = dynamic_cast<FixUVT *>(lmp->modify->get_fix_by_id("cp"));
    ASSERT_NE(fix, nullptr);
    EXPECT_EQ(fix->size_vector, 7 + 4 * chain);
    EXPECT_DOUBLE_EQ(fix->compute_vector(0), 1.0);
    EXPECT_DOUBLE_EQ(fix->compute_vector(1), 31.0 / 4.0);
    EXPECT_DOUBLE_EQ(fix->compute_vector(2), 2.0);
    EXPECT_DOUBLE_EQ(fix->compute_vector(3), 0.0);
    EXPECT_DOUBLE_EQ(fix->compute_vector(4), 0.0);
    EXPECT_DOUBLE_EQ(fix->compute_vector(5), 1.5);
    EXPECT_DOUBLE_EQ(fix->compute_vector(6), 0.0);
    for (int n = 0; n < 4 * chain; ++n) {
      EXPECT_DOUBLE_EQ(fix->compute_vector(n + 7), fix->FixNH::compute_vector(n));
      EXPECT_EQ(fix->get_thermo_colname(n + 7), fix->FixNH::get_thermo_colname(n));
    }
  }
}

TEST_F(ComputeTempUVTTest, FixUsesCombinedTemperatureAndEnergy)
{
  auto *fix = dynamic_cast<FixUVT *>(lmp->modify->get_fix_by_id("cp"));
  ASSERT_NE(fix, nullptr);
  EXPECT_STREQ(fix->temperature->style, "temp/uvt");
  EXPECT_DOUBLE_EQ(fix->tdof, 4.0);
  EXPECT_DOUBLE_EQ(fix->t_current, 31.0 / 4.0);
  EXPECT_DOUBLE_EQ(fix->ke_target, 4.0);
  EXPECT_DOUBLE_EQ(fix->eta_mass[0], 1.0);

  // The first thermostat potential already includes the electronic DOF.
  fix->eta[0] = 0.25;
  EXPECT_DOUBLE_EQ(fix->compute_scalar(), 2.5);

  // The virtual velocity-scaling hook must act on both kinds of velocity.
  fix->factor_eta = 0.5;
  fix->nh_v_temp();
  EXPECT_DOUBLE_EQ(fix->temperature->compute_scalar(), 31.0 / 16.0);
  int dim;
  EXPECT_DOUBLE_EQ(*static_cast<double *>(fix->extract("ne_dot", dim)), 1.0);
}

TEST_F(ComputeTempUVTTest, RestartPreservesCombinedTemperatureAndChainEnergy)
{
  command("run 5 post no");
  auto *fix = dynamic_cast<FixUVT *>(lmp->modify->get_fix_by_id("cp"));
  ASSERT_NE(fix, nullptr);
  const double energy = fix->compute_scalar();
  const double temp = fix->temperature->compute_scalar();
  const double chain_mass = fix->eta_mass[0];
  command("write_restart temp_uvt.restart");
  command("clear");
  command("read_restart temp_uvt.restart");
  command("variable deriv equal 0.0");
  command("fix cp all uvt temp 1 1 0.5 mu 0 0 0.5 ne 1 dedn v_deriv");
  command("run 0 post no");
  platform::unlink("temp_uvt.restart");
  fix = dynamic_cast<FixUVT *>(lmp->modify->get_fix_by_id("cp"));
  ASSERT_NE(fix, nullptr);
  EXPECT_DOUBLE_EQ(fix->eta_mass[0], chain_mass);
  EXPECT_NEAR(fix->compute_scalar(), energy, 1.0e-12);
  EXPECT_NEAR(fix->temperature->compute_scalar(), temp, 1.0e-12);
  EXPECT_NEAR(fix->compute_vector(1), temp, 1.0e-12);
}

TEST_F(ComputeTempUVTTest, FixedMassesArePreservedDuringIntegration)
{
  auto *fix = dynamic_cast<FixUVT *>(lmp->modify->get_fix_by_id("cp"));
  ASSERT_NE(fix, nullptr);
  int dim;
  auto *mass = static_cast<double *>(fix->extract("ne_mass", dim));
  *mass = 1.5;
  fix->eta_mass[0] = 2.0;
  fix->eta_mass_flag = 0;
  fix->t_current = fix->temperature->compute_scalar();
  fix->initial_integrate(0);
  fix->final_integrate();
  EXPECT_DOUBLE_EQ(*mass, 1.5);
  EXPECT_DOUBLE_EQ(fix->eta_mass[0], 2.0);
  EXPECT_NEAR(fix->t_current, fix->temperature->compute_scalar(), 1.0e-12);
}

TEST_F(ComputeTempUVTTest, ThermostatRejectsIncompatibleTemperatureComputes)
{
  command("fix_modify cp temp nuclear");
  EXPECT_ANY_THROW(command("run 0 post no"));
  command("fix_modify cp temp combined");
  EXPECT_NO_THROW(command("run 0 post no"));

  command("fix other all uvt temp 1 1 0.5 mu 0 0 0.5 ne 1 dedn v_deriv");
  command("fix_modify cp temp other_temp");
  EXPECT_ANY_THROW(command("run 0 post no"));
}

TEST_F(ComputeTempUVTTest, DynamicDOFAndExtraDOF)
{
  auto *temp = lmp->modify->get_compute_by_id("combined");
  command("compute_modify combined dynamic/dof yes");
  temp->setup();
  command("create_atoms 1 single 2 2 2");
  command("velocity all set 1 2 3");
  // Ne_mass stays at its prior value until FixUVT updates it.
  EXPECT_DOUBLE_EQ(temp->compute_scalar(), 45.0 / 7.0);
  EXPECT_DOUBLE_EQ(temp->dof, 7.0);
  EXPECT_DOUBLE_EQ(temp->compute_scalar(), 45.0 / 7.0);

  command("compute_modify combined extra/dof 9");
  EXPECT_DOUBLE_EQ(temp->compute_scalar(), 45.0);
  EXPECT_DOUBLE_EQ(temp->dof, 1.0);
  command("compute_modify combined extra/dof 10");
  EXPECT_ANY_THROW(temp->compute_scalar());
}

TEST_F(ComputeTempUVTTest, ValidatesFixAndGroup)
{
  EXPECT_ANY_THROW(command("compute bad all temp/uvt"));
  EXPECT_ANY_THROW(command("compute bad all temp/uvt cp extra"));
  command("compute missing all temp/uvt absent");
  EXPECT_ANY_THROW(lmp->modify->get_compute_by_id("missing")->init());
  command("fix other all nve");
  command("compute wrong all temp/uvt other");
  EXPECT_ANY_THROW(lmp->modify->get_compute_by_id("wrong")->init());
  command("group subset id 1");
  command("compute mismatch subset temp/uvt cp");
  EXPECT_ANY_THROW(lmp->modify->get_compute_by_id("mismatch")->init());

  command("unfix cp");
  EXPECT_ANY_THROW(lmp->modify->get_compute_by_id("combined")->init());
  command("fix cp all uvt temp 1 1 0.5 mu 0 0 0.5 ne 1 ne_velocity 3 dedn v_deriv");
  EXPECT_NO_THROW(lmp->modify->get_compute_by_id("combined")->init());
  command("uncompute missing");
  command("uncompute wrong");
  command("uncompute mismatch");
  command("unfix other");
  command("run 0 post no");
  EXPECT_DOUBLE_EQ(lmp->modify->get_compute_by_id("combined")->compute_scalar(), 34.75 / 4.0);
}

}    // namespace LAMMPS_NS

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleMock(&argc, argv);
  int result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
