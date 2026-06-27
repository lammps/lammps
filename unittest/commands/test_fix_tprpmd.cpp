/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with the U.S. Government retains
   certain rights in this software.
------------------------------------------------------------------------- */

#define LAMMPS_LIB_MPI 1

#include "lammps.h"
#include "library.h"

#include "platform.h"

#include "../testing/core.h"
#include "gtest/gtest.h"

#include <cmath>
#include <mpi.h>

bool verbose = false;

namespace LAMMPS_NS {

class FixTPRPMDSerialTest : public LAMMPSTest {
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

  void setup_quadratic_fix()
  {
    command("variable k_quad equal 5.0");
    command("variable N0_quad equal 1.0");
    command("variable dEdN equal v_k_quad*(f_cp[23]-v_N0_quad)");
    command("fix cp all tprpmd method tp-rpmd ensemble uvt thermostat NHC temp 1.0 Tdamp 0.5 "
            "tchain 3 tloop 1 mu 2.0 2.0 0.5 ne 1.8 ne_velocity 0.0 dedn v_dEdN");
  }

  double fix_value(const char *id, int index)
  {
    void *ptr = lammps_extract_fix(lmp, id, LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, index, 0);
    double value = *(double *) ptr;
    lammps_free(ptr);
    return value;
  }
};

TEST_F(FixTPRPMDSerialTest, RestartRestoresElectronState)
{
  double ne_before = 0.0;
  double nedot_before = 0.0;
  double dedn_before = 0.0;
  double mu_before = 0.0;

  setup_quadratic_system();
  setup_quadratic_fix();
  command("run 20 post no");

  ne_before = fix_value("cp", 22);
  nedot_before = fix_value("cp", 23);
  dedn_before = fix_value("cp", 24);
  mu_before = fix_value("cp", 25);

  command("write_restart tprpmd.restart");
  command("clear");
  command("read_restart tprpmd.restart");
  command("variable k_quad equal 5.0");
  command("variable N0_quad equal 1.0");
  command("variable dEdN equal v_k_quad*(f_cp[23]-v_N0_quad)");
  command("fix cp all tprpmd method tp-rpmd ensemble uvt thermostat NHC temp 1.0 Tdamp 0.5 "
          "tchain 3 tloop 1 mu 2.0 2.0 0.5 ne 1.8 ne_velocity 0.0 dedn v_dEdN");
  command("run 0 post no");
  platform::unlink("tprpmd.restart");

  EXPECT_NEAR(fix_value("cp", 22), ne_before, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", 23), nedot_before, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", 24), dedn_before, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", 25), mu_before, 1.0e-10);
}

TEST_F(FixTPRPMDSerialTest, QuadraticToyFixedPointHasZeroElectronicForce)
{
  setup_quadratic_system();
  command("variable k_quad equal 5.0");
  command("variable N0_quad equal 1.0");
  command("variable dEdN equal v_k_quad*(f_cp[23]-v_N0_quad)");
  command("fix cp all tprpmd method tp-rpmd ensemble uvt thermostat NHC temp 1.0 Tdamp 0.5 "
          "tchain 3 tloop 1 mu 2.0 2.0 0.5 ne 1.4 ne_velocity 0.0 dedn v_dEdN");
  command("run 0 post no");

  EXPECT_NEAR(fix_value("cp", 22), 1.4, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", 23), 0.0, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", 24), 2.0, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", 25), 2.0, 1.0e-12);
}

TEST_F(FixTPRPMDSerialTest, QuadraticToyPhysicsAveragesConverge)
{
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
  command("variable k_quad equal 5.0");
  command("variable N0_quad equal 1.0");
  command("variable dEdN equal v_k_quad*(f_cp[23]-v_N0_quad)");
  command("fix cp all tprpmd method tp-rpmd ensemble uvt thermostat NHC temp 1.0 Tdamp 0.5 "
          "tchain 3 tloop 1 mu 2.0 2.0 0.5 ne 1.8 ne_velocity 0.0 dedn v_dEdN");
  command("fix avg all ave/time 1 20000 20000 f_cp[23] f_cp[24] f_cp[25]");
  command("run 20000 post no");

  EXPECT_NEAR(fix_value("avg", 0), 1.4, 1.0e-1);
  EXPECT_NEAR(fix_value("avg", 1), 0.0, 1.0e-1);
  EXPECT_NEAR(fix_value("avg", 2), 2.0, 1.0e-1);
}

TEST_F(FixTPRPMDSerialTest, TPRPMDP1RetainsElectronicThermostat)
{
  setup_quadratic_system();
  command("variable k_quad equal 5.0");
  command("variable N0_quad equal 1.0");
  command("variable dEdN equal v_k_quad*(f_cp[23]-v_N0_quad)");
  command("fix cp all tprpmd method tp-rpmd ensemble uvt thermostat NHC temp 1.0 Tdamp 0.5 "
          "tchain 3 tloop 1 mu 2.0 2.0 0.5 ne 1.8 ne_velocity 0.4 dedn v_dEdN");
  command("run 20 post no");

  EXPECT_TRUE(std::isfinite(fix_value("cp", 10)));
  EXPECT_TRUE(std::isfinite(fix_value("cp", 13)));
  EXPECT_GT(std::fabs(fix_value("cp", 13)), 1.0e-12);
}

TEST(FixTPRPMDMPI, PartitionedRunExercisesBeadExpansion)
{
  int nprocs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (nprocs != 2) GTEST_SKIP() << "This test requires exactly 2 MPI ranks";

  const char *args[] = {"LAMMPS_test", "-log", "none", "-partition", "2x1", "-echo",
                        "screen",      "-nocite",       "-in",        "none", nullptr};
  char **argv = (char **) args;
  int argc = (sizeof(args) / sizeof(char *)) - 1;

  void *lmp = nullptr;
  ASSERT_NO_THROW(lmp = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr));
  ASSERT_NE(lmp, nullptr);

  auto command = [lmp](const char *line) { lammps_command(lmp, line); };
  auto fix_value = [lmp](const char *id, int index) {
    void *ptr = lammps_extract_fix(lmp, id, LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, index, 0);
    double value = *(double *) ptr;
    lammps_free(ptr);
    return value;
  };

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
  command("variable dEdN equal v_k_quad*(f_cp[23]-v_N0_quad)");
  command("fix cp all tprpmd method tp-rpmd ensemble uvt thermostat NHC temp 0.8 Tdamp 0.2 "
          "tchain 3 tloop 1 mu 1.5 1.5 0.2 ne 1.8 ne_velocity 0.0 dedn v_dEdN");
  command("run 1 post no");

  EXPECT_TRUE(std::isfinite(fix_value("cp", 22)));
  EXPECT_TRUE(std::isfinite(fix_value("cp", 23)));
  EXPECT_TRUE(std::isfinite(fix_value("cp", 24)));
  EXPECT_NEAR(fix_value("cp", 25), 1.5, 1.0e-12);
  EXPECT_NE(fix_value("cp", 24), 0.0);

  lammps_close(lmp);
}

TEST(FixTPRPMDMPI, BeadAveragedQuadraticDerivativeAtFixedPoint)
{
  int nprocs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (nprocs != 2) GTEST_SKIP() << "This test requires exactly 2 MPI ranks for bead averaging";

  const char *args[] = {"LAMMPS_test", "-log", "none", "-partition", "2x1", "-echo",
                        "screen",      "-nocite",       "-in",        "none", nullptr};
  char **argv = (char **) args;
  int argc = (sizeof(args) / sizeof(char *)) - 1;

  void *lmp = nullptr;
  ASSERT_NO_THROW(lmp = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr));
  ASSERT_NE(lmp, nullptr);

  auto command = [lmp](const char *line) { lammps_command(lmp, line); };
  auto fix_value = [lmp](const char *id, int index) {
    void *ptr = lammps_extract_fix(lmp, id, LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, index, 0);
    double value = *(double *) ptr;
    lammps_free(ptr);
    return value;
  };

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
  command("variable dEdN equal v_k_quad*(f_cp[23]-v_N0_quad)");
  command("fix cp all tprpmd method tp-rpmd ensemble uvt thermostat NHC temp 0.8 Tdamp 0.2 "
          "tchain 3 tloop 1 mu 1.5 1.5 0.2 ne 1.7 ne_velocity 0.0 dedn v_dEdN");
  command("run 0 post no");

  EXPECT_NEAR(fix_value("cp", 22), 1.7, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", 23), 0.0, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", 24), 1.5, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", 25), 1.5, 1.0e-12);

  lammps_close(lmp);
}

TEST(FixTPRPMDMPI, TPRPMDCentroidThermostatIsSuppressed)
{
  int nprocs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (nprocs != 2) GTEST_SKIP() << "This test requires exactly 2 MPI ranks for 2-bead TP-RPMD";

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const char *args[] = {"LAMMPS_test", "-log", "none", "-partition", "2x1", "-echo",
                        "screen",      "-nocite",       "-in",        "none", nullptr};
  char **argv = (char **) args;
  int argc = (sizeof(args) / sizeof(char *)) - 1;

  void *lmp = nullptr;
  ASSERT_NO_THROW(lmp = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr));
  ASSERT_NE(lmp, nullptr);

  auto command = [lmp](const char *line) { lammps_command(lmp, line); };
  auto fix_value = [lmp](const char *id, int index) {
    void *ptr = lammps_extract_fix(lmp, id, LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, index, 0);
    double value = *(double *) ptr;
    lammps_free(ptr);
    return value;
  };

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
  command("variable dEdN equal v_k_quad*(f_cp[23]-v_N0_quad)");
  command("fix cp all tprpmd method tp-rpmd ensemble uvt thermostat NHC temp 0.8 Tdamp 0.2 "
          "tchain 3 tloop 1 mu 1.5 1.5 0.2 ne 1.8 ne_velocity 0.0 dedn v_dEdN");
  command("run 10");

  if (rank == 0) {
    EXPECT_NEAR(fix_value("cp", 10), 0.0, 1.0e-14);
    EXPECT_NEAR(fix_value("cp", 11), 0.0, 1.0e-14);
    EXPECT_NEAR(fix_value("cp", 12), 0.0, 1.0e-14);
    EXPECT_NEAR(fix_value("cp", 13), 0.0, 1.0e-14);
    EXPECT_NEAR(fix_value("cp", 14), 0.0, 1.0e-14);
    EXPECT_NEAR(fix_value("cp", 15), 0.0, 1.0e-14);
  } else {
    EXPECT_TRUE(std::isfinite(fix_value("cp", 10)));
    EXPECT_TRUE(std::isfinite(fix_value("cp", 11)));
    EXPECT_TRUE(std::isfinite(fix_value("cp", 12)));
    EXPECT_TRUE(std::isfinite(fix_value("cp", 13)));
    EXPECT_TRUE(std::isfinite(fix_value("cp", 14)));
    EXPECT_TRUE(std::isfinite(fix_value("cp", 15)));
  }

  lammps_close(lmp);
}

TEST(FixTPRPMDMPI, P4LongTimeConvergence)
{
  int nprocs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (nprocs != 4) GTEST_SKIP() << "This test requires exactly 4 MPI ranks for 4-bead TP-RPMD";

  const char *args[] = {"LAMMPS_test", "-log", "none", "-partition", "4x1", "-echo",
                        "screen",      "-nocite",       "-in",        "none", nullptr};
  char **argv = (char **) args;
  int argc = (sizeof(args) / sizeof(char *)) - 1;

  void *lmp = nullptr;
  ASSERT_NO_THROW(lmp = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr));
  ASSERT_NE(lmp, nullptr);

  auto command = [lmp](const char *line) { lammps_command(lmp, line); };
  auto fix_value = [lmp](const char *id, int index) {
    void *ptr = lammps_extract_fix(lmp, id, LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, index, 0);
    double value = *(double *) ptr;
    lammps_free(ptr);
    return value;
  };

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
  command("variable dEdN equal v_k_quad*(f_cp[23]-v_N0_quad)");
  command("fix cp all tprpmd method tp-rpmd ensemble uvt thermostat NHC temp 1.0 Tdamp 0.5 "
          "tchain 3 tloop 1 mu 2.0 2.0 0.5 ne 1.8 ne_velocity 0.0 dedn v_dEdN");
  command("fix avg all ave/time 1 20000 20000 f_cp[23] f_cp[24] f_cp[25]");
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
