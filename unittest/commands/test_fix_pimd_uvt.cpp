/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with the U.S. Government retains
   certain rights in this software.
------------------------------------------------------------------------- */

#define LAMMPS_LIB_MPI 1

#include "lammps.h"
#include "platform.h"

#include "../testing/core.h"
#include "gtest/gtest.h"
#include "pimd_test_utils.h"

#include <cmath>
#include <mpi.h>
#include <vector>

bool verbose = false;

namespace LAMMPS_NS {

class FixPIMDUVTSerialTest : public LAMMPSTest {
 protected:
  void setup_quadratic_system()
  {
    command("units lj");
    command("atom_style atomic");
    command("atom_modify map yes");
    command("boundary p p p");
    command("lattice fcc 0.8442");
    command("region box block 0 3 0 3 0 3");
    command("create_box 1 box");
    command("create_atoms 1 box");
    command("mass 1 1.0");
    command("pair_style zero 2.5");
    command("pair_coeff * *");
    command("neighbor 0.3 bin");
    command("neigh_modify every 1 delay 0 check yes");
    command("velocity all create 1.0 492845 mom yes rot no dist gaussian");
    command("timestep 0.005");
  }

  void setup_uvt_fix(const char *style, double ne, double ne_velocity, double mu_start = 2.0,
                     double mu_stop = 2.0, double mu_damp = 0.5)
  {
    command("variable k_quad equal 5.0");
    command("variable N0_quad equal 1.0");
    command(pimd_test::quadratic_dedn_variable());
    command(std::string("fix cp all ") + style +
            " method tp-rpmd ensemble uvt thermostat NHC temp 1.0 Tdamp 0.5 "
            "tchain 3 tloop 1 mu " +
            std::to_string(mu_start) + " " + std::to_string(mu_stop) + " " +
            std::to_string(mu_damp) + " ne " + std::to_string(ne) + " ne_velocity " +
            std::to_string(ne_velocity) + " dedn v_dEdN");
  }

  double fix_value(const char *id, int index)
  {
    return pimd_test::fix_value(lmp, id, index);
  }
};

TEST_F(FixPIMDUVTSerialTest, RestartRestoresElectronState)
{
  const auto uvt = pimd_test::uvt_vector_indices();
  double ne_before = 0.0;
  double nedot_before = 0.0;
  double dedn_before = 0.0;
  double mu_before = 0.0;

  setup_quadratic_system();
  setup_uvt_fix("pimd/uvt", 1.8, 0.0);
  command("run 20 post no");

  ne_before = fix_value("cp", uvt.ne);
  nedot_before = fix_value("cp", uvt.ne_dot);
  dedn_before = fix_value("cp", uvt.dedn);
  mu_before = fix_value("cp", uvt.mu);

  command("write_restart pimd_uvt.restart");
  command("clear");
  command("read_restart pimd_uvt.restart");
  setup_uvt_fix("pimd/uvt", 1.8, 0.0);
  command("run 0 post no");
  platform::unlink("pimd_uvt.restart");

  EXPECT_NEAR(fix_value("cp", uvt.ne), ne_before, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", uvt.ne_dot), nedot_before, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", uvt.dedn), dedn_before, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", uvt.mu), mu_before, 1.0e-10);
}

TEST_F(FixPIMDUVTSerialTest, QuadraticToyFixedPointHasZeroElectronicForce)
{
  const auto uvt = pimd_test::uvt_vector_indices();

  setup_quadratic_system();
  setup_uvt_fix("pimd/uvt", 1.4, 0.0);
  command("run 0 post no");

  EXPECT_NEAR(fix_value("cp", uvt.ne), 1.4, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", uvt.ne_dot), 0.0, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", uvt.dedn), 2.0, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", uvt.mu), 2.0, 1.0e-12);
}

TEST_F(FixPIMDUVTSerialTest, QuadraticToyPhysicsAveragesConverge)
{
  const auto uvt = pimd_test::uvt_vector_indices();

  command("units lj");
  command("atom_style atomic");
  command("atom_modify map yes");
  command("boundary p p p");
  command("lattice sc 0.7");
  command("region box block 0 2 0 2 0 2");
  command("create_box 1 box");
  command("create_atoms 1 box");
  command("mass 1 1.0");
  command("pair_style zero 2.5");
  command("pair_coeff * *");
  command("neighbor 0.3 bin");
  command("neigh_modify every 1 delay 0 check yes");
  command("velocity all create 1.0 492845 mom yes rot no dist gaussian");
  command("timestep 0.005");
  setup_uvt_fix("pimd/uvt", 1.8, 0.0);
  command(("fix avg all ave/time 1 20000 20000 f_cp[" +
           std::to_string(pimd_test::lammps_fix_index(uvt.ne)) + "] f_cp[" +
           std::to_string(pimd_test::lammps_fix_index(uvt.ne_dot)) + "] f_cp[" +
           std::to_string(pimd_test::lammps_fix_index(uvt.dedn)) + "]")
              .c_str());
  command("run 20000 post no");

  EXPECT_NEAR(fix_value("avg", 0), 1.4, 1.0e-1);
  EXPECT_NEAR(fix_value("avg", 1), 0.0, 1.0e-1);
  EXPECT_NEAR(fix_value("avg", 2), 2.0, 1.0e-1);
}

TEST_F(FixPIMDUVTSerialTest, P1RetainsElectronicThermostat)
{
  setup_quadratic_system();
  setup_uvt_fix("pimd/uvt", 1.8, 0.4);
  command("run 20 post no");

  EXPECT_TRUE(std::isfinite(fix_value("cp", 10)));
  EXPECT_TRUE(std::isfinite(fix_value("cp", 13)));
  EXPECT_GT(std::fabs(fix_value("cp", 13)), 1.0e-12);
}

TEST(FixPIMDUVTMPI, PartitionedRunExercisesBeadExpansion)
{
  int nprocs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (nprocs != 2) GTEST_SKIP() << "This test requires exactly 2 MPI ranks";

  const auto uvt = pimd_test::uvt_vector_indices();

  const char *args[] = {"LAMMPS_test", "-log", "none", "-partition", "2x1", "-echo",
                        "screen",      "-nocite",       "-in",        "none", nullptr};
  char **argv = (char **) args;
  int argc = (sizeof(args) / sizeof(char *)) - 1;

  void *lmp = nullptr;
  ASSERT_NO_THROW(lmp = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr));
  ASSERT_NE(lmp, nullptr);

  auto command = [lmp](const char *line) { lammps_command(lmp, line); };
  auto fix_value = [lmp](const char *id, int index) { return pimd_test::fix_value(lmp, id, index); };

  command("units lj");
  command("atom_style atomic");
  command("atom_modify map yes");
  command("boundary p p p");
  command("lattice sc 0.7");
  command("region box block 0 2 0 2 0 2");
  command("create_box 1 box");
  command("create_atoms 1 box");
  command("mass 1 1.0");
  command("pair_style zero 2.5");
  command("pair_coeff * *");
  command("neighbor 0.3 bin");
  command("neigh_modify every 1 delay 0 check yes");
  command("timestep 0.002");
  command("variable beadshift universe 0.0 0.15");
  command("displace_atoms all move ${beadshift} 0.0 0.0 units box");
  command("velocity all create 0.8 97531 mom yes rot no dist gaussian");
  command("variable k_quad equal 3.0");
  command("variable N0_quad equal 1.2");
  command(pimd_test::quadratic_dedn_variable().c_str());
  command("fix cp all pimd/uvt method tp-rpmd ensemble uvt thermostat NHC temp 0.8 Tdamp 0.2 "
          "tchain 3 tloop 1 mu 1.5 1.5 0.2 ne 1.8 ne_velocity 0.0 dedn v_dEdN");
  command("run 1 post no");

  EXPECT_TRUE(std::isfinite(fix_value("cp", uvt.ne)));
  EXPECT_TRUE(std::isfinite(fix_value("cp", uvt.ne_dot)));
  EXPECT_TRUE(std::isfinite(fix_value("cp", uvt.dedn)));
  EXPECT_NEAR(fix_value("cp", uvt.mu), 1.5, 1.0e-12);
  EXPECT_NE(fix_value("cp", uvt.dedn), 0.0);

  lammps_close(lmp);
}

TEST(FixPIMDUVTMPI, BeadAveragedQuadraticDerivativeAtFixedPoint)
{
  int nprocs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (nprocs != 2) GTEST_SKIP() << "This test requires exactly 2 MPI ranks for bead averaging";

  const auto uvt = pimd_test::uvt_vector_indices();

  const char *args[] = {"LAMMPS_test", "-log", "none", "-partition", "2x1", "-echo",
                        "screen",      "-nocite",       "-in",        "none", nullptr};
  char **argv = (char **) args;
  int argc = (sizeof(args) / sizeof(char *)) - 1;

  void *lmp = nullptr;
  ASSERT_NO_THROW(lmp = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr));
  ASSERT_NE(lmp, nullptr);

  auto command = [lmp](const char *line) { lammps_command(lmp, line); };
  auto fix_value = [lmp](const char *id, int index) { return pimd_test::fix_value(lmp, id, index); };

  command("units lj");
  command("atom_style atomic");
  command("atom_modify map yes");
  command("boundary p p p");
  command("lattice sc 0.7");
  command("region box block 0 2 0 2 0 2");
  command("create_box 1 box");
  command("create_atoms 1 box");
  command("mass 1 1.0");
  command("pair_style zero 2.5");
  command("pair_coeff * *");
  command("neighbor 0.3 bin");
  command("neigh_modify every 1 delay 0 check yes");
  command("timestep 0.002");
  command("variable beadshift universe 0.0 0.15");
  command("displace_atoms all move ${beadshift} 0.0 0.0 units box");
  command("velocity all create 0.8 97531 mom yes rot no dist gaussian");
  command("variable N0_quad universe 1.0 1.4");
  command("variable k_quad equal 3.0");
  command(pimd_test::quadratic_dedn_variable().c_str());
  command("fix cp all pimd/uvt method tp-rpmd ensemble uvt thermostat NHC temp 0.8 Tdamp 0.2 "
          "tchain 3 tloop 1 mu 1.5 1.5 0.2 ne 1.7 ne_velocity 0.0 dedn v_dEdN");
  command("run 0 post no");

  EXPECT_NEAR(fix_value("cp", uvt.ne), 1.7, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", uvt.ne_dot), 0.0, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", uvt.dedn), 1.5, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", uvt.mu), 1.5, 1.0e-12);

  lammps_close(lmp);
}

TEST(FixPIMDUVTMPI, TPRPMDCentroidThermostatIsSuppressed)
{
  int nprocs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (nprocs != 2) GTEST_SKIP() << "This test requires exactly 2 MPI ranks for 2-bead TP-RPMD";

  const auto uvt = pimd_test::uvt_vector_indices();
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const char *args[] = {"LAMMPS_test", "-log", "none", "-partition", "2x1", "-echo",
                        "screen",      "-nocite",       "-in",        "none", nullptr};
  char **argv = (char **) args;
  int argc = (sizeof(args) / sizeof(char *)) - 1;

  void *lmp = nullptr;
  ASSERT_NO_THROW(lmp = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr));
  ASSERT_NE(lmp, nullptr);

  auto command = [lmp](const char *line) { lammps_command(lmp, line); };
  auto fix_value = [lmp](const char *id, int index) { return pimd_test::fix_value(lmp, id, index); };

  command("units lj");
  command("atom_style atomic");
  command("atom_modify map yes");
  command("boundary p p p");
  command("lattice sc 0.7");
  command("region box block 0 2 0 2 0 2");
  command("create_box 1 box");
  command("create_atoms 1 box");
  command("mass 1 1.0");
  command("pair_style zero 2.5");
  command("pair_coeff * *");
  command("neighbor 0.3 bin");
  command("neigh_modify every 1 delay 0 check yes");
  command("timestep 0.002");
  command("variable beadshift universe 0.0 0.15");
  command("displace_atoms all move ${beadshift} 0.0 0.0 units box");
  command("velocity all create 0.8 97531 mom yes rot no dist gaussian");
  command("variable k_quad equal 3.0");
  command("variable N0_quad equal 1.2");
  command(pimd_test::quadratic_dedn_variable().c_str());
  command("fix cp all pimd/uvt method tp-rpmd ensemble uvt thermostat NHC temp 0.8 Tdamp 0.2 "
          "tchain 3 tloop 1 mu 1.5 1.5 0.2 ne 1.8 ne_velocity 0.0 dedn v_dEdN");
  command("run 10");

  if (rank == 0) {
    for (int i = 10; i <= 15; ++i) EXPECT_NEAR(fix_value("cp", i), 0.0, 1.0e-14) << "index=" << i;
  } else {
    for (int i = 10; i <= 15; ++i) EXPECT_TRUE(std::isfinite(fix_value("cp", i))) << "index=" << i;
  }

  EXPECT_TRUE(std::isfinite(fix_value("cp", uvt.ne_dot)));

  lammps_close(lmp);
}

TEST(FixPIMDUVTMPI, P4LongTimeConvergence)
{
  int nprocs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (nprocs != 4) GTEST_SKIP() << "This test requires exactly 4 MPI ranks for 4-bead TP-RPMD";

  const auto uvt = pimd_test::uvt_vector_indices();

  const char *args[] = {"LAMMPS_test", "-log", "none", "-partition", "4x1", "-echo",
                        "screen",      "-nocite",       "-in",        "none", nullptr};
  char **argv = (char **) args;
  int argc = (sizeof(args) / sizeof(char *)) - 1;

  void *lmp = nullptr;
  ASSERT_NO_THROW(lmp = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr));
  ASSERT_NE(lmp, nullptr);

  auto command = [lmp](const char *line) { lammps_command(lmp, line); };
  auto fix_value = [lmp](const char *id, int index) { return pimd_test::fix_value(lmp, id, index); };

  command("units lj");
  command("atom_style atomic");
  command("atom_modify map yes");
  command("boundary p p p");
  command("lattice sc 0.7");
  command("region box block 0 2 0 2 0 2");
  command("create_box 1 box");
  command("create_atoms 1 box");
  command("mass 1 1.0");
  command("pair_style zero 2.5");
  command("pair_coeff * *");
  command("neighbor 0.3 bin");
  command("neigh_modify every 1 delay 0 check yes");
  command("timestep 0.005");
  command("variable beadshift universe 0.0 0.1 0.2 0.3");
  command("displace_atoms all move ${beadshift} 0.0 0.0 units box");
  command("velocity all create 1.0 492845 mom yes rot no dist gaussian");
  command("variable k_quad equal 5.0");
  command("variable N0_quad equal 1.0");
  command(pimd_test::quadratic_dedn_variable().c_str());
  command("fix cp all pimd/uvt method tp-rpmd ensemble uvt thermostat NHC temp 1.0 Tdamp 0.5 "
          "tchain 3 tloop 1 mu 2.0 2.0 0.5 ne 1.8 ne_velocity 0.0 dedn v_dEdN");
  command(("fix avg all ave/time 1 20000 20000 f_cp[" +
           std::to_string(pimd_test::lammps_fix_index(uvt.ne)) + "] f_cp[" +
           std::to_string(pimd_test::lammps_fix_index(uvt.ne_dot)) + "] f_cp[" +
           std::to_string(pimd_test::lammps_fix_index(uvt.dedn)) + "]")
              .c_str());
  command("run 20000 post no");

  EXPECT_NEAR(fix_value("avg", 0), 1.4, 1.0e-1);
  EXPECT_NEAR(fix_value("avg", 1), 0.0, 1.0e-1);
  EXPECT_NEAR(fix_value("avg", 2), 2.0, 1.0e-1);

  lammps_close(lmp);
}

}    // namespace LAMMPS_NS

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleMock(&argc, argv);

  if (const char *var = getenv("TEST_ARGS")) {
    std::vector<std::string> env = LAMMPS_NS::utils::split_words(var);
    for (auto arg : env) {
      if (arg == "-v") verbose = true;
    }
  }

  if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

  int rv = RUN_ALL_TESTS();
  MPI_Finalize();
  return rv;
}
