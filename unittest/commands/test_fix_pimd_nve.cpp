/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with the U.S. Government retains
   certain rights in this software.
------------------------------------------------------------------------- */

#define LAMMPS_LIB_MPI 1

#include "atom.h"
#include "lammps.h"
#include "platform.h"

#include "../testing/core.h"
#include "gtest/gtest.h"
#include "pimd_test_utils.h"

#include <cmath>
#include <mpi.h>
#include <string>
#include <vector>

bool verbose = false;

namespace LAMMPS_NS {

class FixPIMDNVESerialTest : public LAMMPSTest {
 protected:
  void setup_zero_pair_system()
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
  }

  void setup_nve_fix(const char *method = "nmpimd")
  {
    command(std::string("fix cp all pimd/nve method ") + method + " temp 1.0");
  }

  double fix_value(const char *id, int index)
  {
    return pimd_test::fix_value(lmp, id, index);
  }
};

TEST_F(FixPIMDNVESerialTest, DoesNotMoveAtomsOutsideFixGroup)
{
  command("units lj");
  command("atom_style atomic");
  command("atom_modify map yes");
  command("boundary p p p");
  command("region box block 0 4 0 4 0 4");
  command("create_box 2 box");
  command("create_atoms 1 single 1.0 1.0 1.0");
  command("create_atoms 2 single 2.5 2.5 2.5");
  command("mass * 1.0");
  command("group mobile type 1");
  command("pair_style zero 2.5");
  command("pair_coeff * *");
  command("neighbor 0.3 bin");
  command("neigh_modify every 1 delay 0 check yes");
  command("velocity all set 0.2 0.1 -0.05");
  command("timestep 0.002");

  auto *atom = lmp->atom;
  double initial[3] = {0.0, 0.0, 0.0};
  bool found = false;
  for (int i = 0; i < atom->nlocal; ++i) {
    if (atom->type[i] == 2) {
      initial[0] = atom->x[i][0];
      initial[1] = atom->x[i][1];
      initial[2] = atom->x[i][2];
      found = true;
      break;
    }
  }
  ASSERT_TRUE(found);

  command("fix cp mobile pimd/nve temp 1.0");
  command("run 10 post no");

  found = false;
  for (int i = 0; i < atom->nlocal; ++i) {
    if (atom->type[i] == 2) {
      EXPECT_NEAR(atom->x[i][0], initial[0], 1.0e-12);
      EXPECT_NEAR(atom->x[i][1], initial[1], 1.0e-12);
      EXPECT_NEAR(atom->x[i][2], initial[2], 1.0e-12);
      found = true;
      break;
    }
  }
  ASSERT_TRUE(found);
}

TEST_F(FixPIMDNVESerialTest, P1StandaloneRunProducesFiniteVector)
{
  setup_zero_pair_system();
  setup_nve_fix();
  command("run 1 post no");

  for (int i = 0; i < 10; ++i) {
    EXPECT_TRUE(std::isfinite(fix_value("cp", i))) << "index=" << i;
  }
  EXPECT_GT(fix_value("cp", 0), 0.0);
}

TEST_F(FixPIMDNVESerialTest, PIMDStyleParsesAndRuns)
{
  setup_zero_pair_system();
  setup_nve_fix("pimd");
  command("run 10 post no");

  for (int i = 0; i < 10; ++i) {
    EXPECT_TRUE(std::isfinite(fix_value("cp", i))) << "index=" << i;
  }
}

TEST_F(FixPIMDNVESerialTest, RejectsCMDMethod)
{
  setup_zero_pair_system();
  EXPECT_ANY_THROW(command("fix cp all pimd/nve method cmd temp 1.0"));
}

TEST_F(FixPIMDNVESerialTest, RestartRestoresNuclearPrefix)
{
  std::vector<double> before(10);

  setup_zero_pair_system();
  setup_nve_fix();
  command("run 10 post no");
  for (int i = 0; i < 10; ++i) before[i] = fix_value("cp", i);

  command("write_restart pimd_nve.restart");
  command("clear");
  command("read_restart pimd_nve.restart");
  setup_nve_fix();
  command("run 0 post no");
  platform::unlink("pimd_nve.restart");

  for (int i = 0; i < 10; ++i) {
    EXPECT_NEAR(fix_value("cp", i), before[i], 1.0e-10) << "index=" << i;
  }
}

TEST(FixPIMDNVEMPI, PartitionedRunExercisesBeadExpansion)
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
  command("fix cp all pimd/nve temp 0.8");
  command("run 10 post no");

  for (int i = 0; i < 10; ++i) {
    EXPECT_TRUE(std::isfinite(fix_value("cp", i))) << "index=" << i;
  }

  lammps_close(lmp);
}

TEST(FixPIMDNVEMPI, PIMDPartitionedRunSupportsOneRankPerBead)
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
  command("velocity all create 0.8 13579 mom yes rot no dist gaussian");
  command("fix cp all pimd/nve method pimd temp 0.8");
  command("run 10 post no");

  for (int i = 0; i < 10; ++i) {
    EXPECT_TRUE(std::isfinite(fix_value("cp", i))) << "index=" << i;
  }

  lammps_close(lmp);
}

TEST(FixPIMDNVEMPI, ConservesTotalEnergyOverShortRun)
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
  command("fix cp all pimd/nve temp 0.8");
  command("run 0 post no");

  double energy_before = fix_value("cp", 3);
  command("run 100 post no");
  double energy_after = fix_value("cp", 3);

  EXPECT_NEAR(energy_after, energy_before, 1.0e-8);

  lammps_close(lmp);
}

TEST(FixPIMDNVEMPI, MultiRankPerBeadRunProducesFiniteVector)
{
  int nprocs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (nprocs != 4) GTEST_SKIP() << "This test requires exactly 4 MPI ranks";

  const char *args[] = {"LAMMPS_test", "-log", "none", "-partition", "2x2", "-echo",
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
  command("fix cp all pimd/nve temp 0.8");
  command("run 10 post no");

  for (int i = 0; i < 10; ++i) {
    EXPECT_TRUE(std::isfinite(fix_value("cp", i))) << "index=" << i;
  }

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
