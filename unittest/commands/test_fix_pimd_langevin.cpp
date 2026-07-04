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
#include "pimd_langevin_test_utils.h"

#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <mpi.h>
#include <string>
#include <vector>

bool verbose = false;

namespace LAMMPS_NS {
namespace {

using pimd_langevin_test::fix_value;
using pimd_langevin_test::fix_vector;
using pimd_langevin_test::fix_vector_size;

bool has_exact_procs(int expected_procs)
{
  int nprocs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  return nprocs == expected_procs;
}

void setup_single_atom_zero_pair(const std::function<void(const std::string &)> &command)
{
  command("units lj");
  command("atom_style atomic");
  command("atom_modify map yes");
  command("boundary p p p");
  command("region box block 0 4 0 4 0 4");
  command("create_box 1 box");
  command("create_atoms 1 single 1.0 1.0 1.0");
  command("mass 1 1.0");
  command("pair_style zero 2.5");
  command("pair_coeff * *");
  command("neighbor 0.3 bin");
  command("neigh_modify every 1 delay 0 check yes");
  command("velocity all set 0.2 0.1 -0.05");
  command("timestep 0.002");
}

void setup_two_atom_group_system(const std::function<void(const std::string &)> &command)
{
  command("units lj");
  command("atom_style atomic");
  command("atom_modify map yes");
  command("boundary p p p");
  command("region box block 0 4 0 4 0 4");
  command("create_box 2 box");
  command("create_atoms 1 single 1.0 1.0 1.0");
  command("create_atoms 2 single 3.0 3.0 3.0");
  command("mass 1 1.0");
  command("mass 2 1.0");
  command("group real type 1");
  command("pair_style zero 2.5");
  command("pair_coeff * *");
  command("neighbor 0.3 bin");
  command("neigh_modify every 1 delay 0 check yes");
  command("velocity all set 0.0 0.0 0.0");
  command("set atom 1 vx 0.2 vy 0.0 vz 0.0");
  command("set atom 2 vx 0.4 vy 0.0 vz 0.0");
  command("timestep 0.002");
}

void setup_partition_zero_pair_system(void *lmp, bool two_types = false)
{
  auto command = [lmp](const char *line) { lammps_command(lmp, line); };

  command("units lj");
  command("atom_style atomic");
  command("atom_modify map yes");
  command("boundary p p p");
  command("region box block 0 4 0 4 0 4");
  if (two_types) {
    command("create_box 2 box");
  } else {
    command("create_box 1 box");
  }
  command("create_atoms 1 single 1.0 1.0 1.0");
  if (two_types) command("create_atoms 2 single 2.5 2.5 2.5");
  command("mass 1 1.0");
  if (two_types) command("mass 2 1.0");
  command("pair_style zero 2.5");
  command("pair_coeff * *");
  command("neighbor 0.3 bin");
  command("neigh_modify every 1 delay 0 check yes");
  command("velocity all set 0.2 0.1 -0.05");
  command("timestep 0.002");
  command("variable beadshift universe 0.0 0.15");
  command("displace_atoms all move ${beadshift} 0.0 0.0 units box");
}

void setup_partition_lj_system(void *lmp)
{
  auto command = [lmp](const char *line) { lammps_command(lmp, line); };

  command("units lj");
  command("atom_style atomic");
  command("atom_modify map yes");
  command("boundary p p p");
  command("lattice sc 0.7");
  command("region box block 0 2 0 2 0 2");
  command("create_box 1 box");
  command("create_atoms 1 box");
  command("mass 1 1.0");
  command("pair_style lj/cut 2.5");
  command("pair_coeff * * 1.0 1.0 2.5");
  command("neighbor 0.3 bin");
  command("neigh_modify every 1 delay 0 check yes");
  command("timestep 0.002");
  command("variable beadshift universe 0.0 0.15");
  command("displace_atoms all move ${beadshift} 0.0 0.0 units box");
  command("velocity all create 0.8 97531 mom yes rot no dist gaussian");
}

void setup_barostat_system(void *lmp)
{
  auto command = [lmp](const char *line) { lammps_command(lmp, line); };

  command("units lj");
  command("atom_style atomic");
  command("atom_modify map yes");
  command("boundary p p p");
  command("lattice sc 0.72");
  command("region box block 0 3 0 3 0 3");
  command("create_box 1 box");
  command("create_atoms 1 box");
  command("mass 1 1.0");
  command("pair_style lj/cut 2.5");
  command("pair_coeff * * 1.0 1.0 2.5");
  command("neighbor 0.3 bin");
  command("neigh_modify every 1 delay 0 check yes");
  command("timestep 0.001");
  command("variable beadshift universe 0.0 0.08 0.16 0.24");
  command("displace_atoms all move ${beadshift} 0.0 0.0 units box");
  command("velocity all create 0.8 24680 mom yes rot no dist gaussian");
}

void cleanup_partition_restart(const std::string &prefix, int count)
{
  char suffix[8];
  for (int i = 1; i <= count; ++i) {
    snprintf(suffix, sizeof(suffix), ".%02d", i);
    platform::unlink((prefix + suffix).c_str());
  }
}

void *open_partitioned_lammps(int /*expected_procs*/, const char *partition)
{
  const char *args[] = {"LAMMPS_test", "-log", "none", "-partition", partition,
                        "-echo",       "screen", "-nocite",    "-in",
                        "none",        nullptr};
  char **argv = (char **) args;
  const int argc = (sizeof(args) / sizeof(char *)) - 1;

  void *lmp = nullptr;
  EXPECT_NO_THROW(lmp = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr));
  EXPECT_NE(lmp, nullptr);
  return lmp;
}

void expect_all_finite(void *lmp, const char *fix_id)
{
  const auto values = fix_vector(lmp, fix_id);
  for (size_t i = 0; i < values.size(); ++i) {
    EXPECT_TRUE(std::isfinite(values[i])) << "index=" << i;
  }
}

void expect_reproducible_nvt_partitioned(const std::function<void(void *)> &setup_commands,
                                         const std::function<void(void *)> &fix_commands,
                                         int expected_procs, const char *partition,
                                         double tolerance = 1.0e-12)
{
  std::vector<double> first;

  {
    void *lmp = open_partitioned_lammps(expected_procs, partition);
    ASSERT_NE(lmp, nullptr);
    setup_commands(lmp);
    fix_commands(lmp);
    lammps_command(lmp, "run 2 post no");
    first = fix_vector(lmp, "cp");
    lammps_close(lmp);
  }

  {
    void *lmp = open_partitioned_lammps(expected_procs, partition);
    ASSERT_NE(lmp, nullptr);
    setup_commands(lmp);
    fix_commands(lmp);
    lammps_command(lmp, "run 2 post no");
    const auto second = fix_vector(lmp, "cp");
    ASSERT_EQ(first.size(), second.size());
    for (size_t i = 0; i < first.size(); ++i) {
      EXPECT_NEAR(first[i], second[i], tolerance) << "index=" << i;
    }
    lammps_close(lmp);
  }
}

}    // namespace

class FixPIMDLangevinSerialTest : public LAMMPSTest {
 protected:
  int fix_vector_size(const char *id) { return pimd_langevin_test::fix_vector_size(lmp, id); }
  double fix_value(const char *id, int index) { return pimd_langevin_test::fix_value(lmp, id, index); }
};

TEST_F(FixPIMDLangevinSerialTest, NMPIMDNVEBAOABP1)
{
  setup_single_atom_zero_pair([this](const std::string &line) { command(line); });
  command("fix cp all pimd/langevin method nmpimd ensemble nve integrator baoab temp 1.0 "
          "thermostat PILE_L 1234");
  command("run 0 post no");

  ASSERT_EQ(fix_vector_size("cp"), pimd_langevin_test::kNuclearPrefixScalars);
  EXPECT_NEAR(fix_value("cp", pimd_langevin_test::KE_BEAD), 0.02625, 1.0e-12);
  EXPECT_NEAR(fix_value("cp", pimd_langevin_test::SE_BEAD), 0.0, 1.0e-12);
  EXPECT_NEAR(fix_value("cp", pimd_langevin_test::PE_BEAD), 0.0, 1.0e-12);
  EXPECT_NEAR(fix_value("cp", pimd_langevin_test::TOTAL_ENERGY), 0.02625, 1.0e-12);
  EXPECT_NEAR(fix_value("cp", pimd_langevin_test::T_PRIM), 1.5, 1.0e-12);
  EXPECT_NEAR(fix_value("cp", pimd_langevin_test::T_VIR), 0.0, 1.0e-12);
  EXPECT_NEAR(fix_value("cp", pimd_langevin_test::T_CV), 1.5, 1.0e-12);
  EXPECT_NEAR(fix_value("cp", pimd_langevin_test::P_PRIM), 1.0 / 64.0, 1.0e-12);
  EXPECT_NEAR(fix_value("cp", pimd_langevin_test::P_MD), 0.02625 / 192.0, 1.0e-12);
  EXPECT_NEAR(fix_value("cp", pimd_langevin_test::P_CV), 0.02625 / 96.0, 1.0e-12);
}

TEST_F(FixPIMDLangevinSerialTest, GroupRestrictedEstimators)
{
  setup_two_atom_group_system([this](const std::string &line) { command(line); });
  command("fix cp real pimd/langevin method nmpimd ensemble nve integrator baoab temp 2.0 "
          "thermostat PILE_L 1234");
  command("run 0 post no");

  EXPECT_NEAR(fix_value("cp", pimd_langevin_test::KE_BEAD), 0.02, 1.0e-12);
  EXPECT_NEAR(fix_value("cp", pimd_langevin_test::T_PRIM), 3.0, 1.0e-12);
  EXPECT_NEAR(fix_value("cp", pimd_langevin_test::P_PRIM), 0.03125, 1.0e-12);

  command("clear");
  setup_two_atom_group_system([this](const std::string &line) { command(line); });
  command("fix cp all pimd/langevin method nmpimd ensemble nve integrator baoab temp 2.0 "
          "thermostat PILE_L 1234");
  command("run 0 post no");

  EXPECT_NEAR(fix_value("cp", pimd_langevin_test::KE_BEAD), 0.10, 1.0e-12);
  EXPECT_NEAR(fix_value("cp", pimd_langevin_test::T_PRIM), 6.0, 1.0e-12);
  EXPECT_NEAR(fix_value("cp", pimd_langevin_test::P_PRIM), 0.0625, 1.0e-12);
}

TEST_F(FixPIMDLangevinSerialTest, NuclearPrefixLayout)
{
  setup_single_atom_zero_pair([this](const std::string &line) { command(line); });
  command("fix cp all pimd/langevin method nmpimd ensemble nve integrator baoab temp 1.0 "
          "thermostat PILE_L 1234");
  command("run 0 post no");
  EXPECT_EQ(fix_vector_size("cp"), pimd_langevin_test::kNuclearPrefixScalars);

  command("clear");
  setup_barostat_system(lmp);
  command("fix cp all pimd/langevin method nmpimd ensemble npt integrator obabo temp 0.8 "
          "thermostat PILE_L 1234 tau 1.0 iso 1.0 barostat BZP taup 1.0 fixcom no");
  command("run 0 post no");
  EXPECT_EQ(fix_vector_size("cp"), pimd_langevin_test::kIsoBarostatVectorSize);

  command("clear");
  setup_barostat_system(lmp);
  command("fix cp all pimd/langevin method nmpimd ensemble npt integrator obabo temp 0.8 "
          "thermostat PILE_L 1234 tau 1.0 x 1.0 y 1.0 z 1.0 barostat BZP taup 1.0 fixcom no");
  command("run 0 post no");
  EXPECT_EQ(fix_vector_size("cp"), pimd_langevin_test::kAnisoBarostatVectorSize);
}

TEST(FixPIMDLangevinMPI, NMPIMDNVEOBABOP2)
{
  if (!has_exact_procs(2)) GTEST_SKIP() << "This test requires exactly 2 MPI ranks";
  void *lmp = open_partitioned_lammps(2, "2x1");
  ASSERT_NE(lmp, nullptr);

  setup_partition_zero_pair_system(lmp);
  lammps_command(lmp, "fix cp all pimd/langevin method nmpimd ensemble nve integrator obabo "
                     "temp 1.0 thermostat PILE_L 1234 fixcom no");
  lammps_command(lmp, "run 0 post no");
  double before = fix_value(lmp, "cp", pimd_langevin_test::TOTAL_ENERGY);
  lammps_command(lmp, "run 10 post no");
  double after = fix_value(lmp, "cp", pimd_langevin_test::TOTAL_ENERGY);

  EXPECT_EQ(fix_vector_size(lmp, "cp"), pimd_langevin_test::kNuclearPrefixScalars);
  expect_all_finite(lmp, "cp");
  EXPECT_NEAR(after, before, 1.0e-8);

  lammps_close(lmp);
}

TEST(FixPIMDLangevinMPI, NMPIMDNVTBAOABP2)
{
  if (!has_exact_procs(2)) GTEST_SKIP() << "This test requires exactly 2 MPI ranks";
  expect_reproducible_nvt_partitioned(
      [](void *lmp) { setup_partition_zero_pair_system(lmp); },
      [](void *lmp) {
        lammps_command(lmp, "fix cp all pimd/langevin method nmpimd ensemble nvt integrator "
                           "baoab temp 1.0 thermostat PILE_L 1234 tau 1.0 fixcom no");
      },
      2, "2x1");
}

TEST(FixPIMDLangevinMPI, NMPIMDNVTOBABOMultiRankPerBead)
{
  if (!has_exact_procs(4)) GTEST_SKIP() << "This test requires exactly 4 MPI ranks";
  void *lmp = open_partitioned_lammps(4, "2x2");
  ASSERT_NE(lmp, nullptr);

  setup_partition_lj_system(lmp);
  lammps_command(lmp, "fix cp all pimd/langevin method nmpimd ensemble nvt integrator obabo "
                     "temp 0.8 sp 0.2 thermostat PILE_L 1234 tau 1.0 fixcom no");
  lammps_command(lmp, "run 2 post no");

  EXPECT_EQ(fix_vector_size(lmp, "cp"), pimd_langevin_test::kNuclearPrefixScalars);
  expect_all_finite(lmp, "cp");
  EXPECT_GT(lammps_get_natoms(lmp), 0.0);

  lammps_close(lmp);
}

TEST(FixPIMDLangevinMPI, PIMDNVEOBABOP2)
{
  if (!has_exact_procs(2)) GTEST_SKIP() << "This test requires exactly 2 MPI ranks";
  void *lmp = open_partitioned_lammps(2, "2x1");
  ASSERT_NE(lmp, nullptr);

  setup_partition_zero_pair_system(lmp);
  lammps_command(lmp, "fix cp all pimd/langevin method pimd ensemble nve integrator obabo "
                     "temp 1.0 thermostat PILE_L 1234 fixcom no");
  lammps_command(lmp, "run 2 post no");

  EXPECT_EQ(fix_vector_size(lmp, "cp"), pimd_langevin_test::kNuclearPrefixScalars);
  expect_all_finite(lmp, "cp");
  EXPECT_GT(std::abs(fix_value(lmp, "cp", pimd_langevin_test::SE_BEAD)), 0.0);

  lammps_close(lmp);
}

TEST(FixPIMDLangevinMPI, PIMDNVTBAOABP2)
{
  if (!has_exact_procs(2)) GTEST_SKIP() << "This test requires exactly 2 MPI ranks";
  expect_reproducible_nvt_partitioned(
      [](void *lmp) { setup_partition_zero_pair_system(lmp); },
      [](void *lmp) {
        lammps_command(lmp, "fix cp all pimd/langevin method pimd ensemble nvt integrator "
                           "baoab temp 1.0 thermostat PILE_L 1234 tau 1.0 fixcom no");
      },
      2, "2x1");
}

TEST(FixPIMDLangevinMPI, NMPIMDNPTBZP)
{
  if (!has_exact_procs(4)) GTEST_SKIP() << "This test requires exactly 4 MPI ranks";
  void *lmp = open_partitioned_lammps(4, "4x1");
  ASSERT_NE(lmp, nullptr);

  setup_barostat_system(lmp);
  lammps_command(lmp, "fix cp all pimd/langevin method nmpimd ensemble npt integrator obabo "
                     "temp 0.8 thermostat PILE_L 1234 tau 1.0 iso 1.0 barostat BZP taup 1.0 "
                     "fixcom no");
  lammps_command(lmp, "run 2 post no");

  EXPECT_EQ(fix_vector_size(lmp, "cp"), pimd_langevin_test::kIsoBarostatVectorSize);
  expect_all_finite(lmp, "cp");

  lammps_close(lmp);
}

TEST(FixPIMDLangevinMPI, NMPIMDNPTMTTK)
{
  if (!has_exact_procs(4)) GTEST_SKIP() << "This test requires exactly 4 MPI ranks";
  void *lmp = open_partitioned_lammps(4, "4x1");
  ASSERT_NE(lmp, nullptr);

  setup_barostat_system(lmp);
  lammps_command(lmp, "fix cp all pimd/langevin method nmpimd ensemble npt integrator obabo "
                     "temp 0.8 thermostat PILE_L 1234 tau 1.0 iso 1.0 barostat MTTK taup 1.0 "
                     "fixcom no");
  lammps_command(lmp, "run 2 post no");

  EXPECT_EQ(fix_vector_size(lmp, "cp"), pimd_langevin_test::kIsoBarostatVectorSize);
  expect_all_finite(lmp, "cp");

  lammps_close(lmp);
}

TEST(FixPIMDLangevinRestartTest, NMPIMDNVTRestartContinuation)
{
  if (!has_exact_procs(2)) GTEST_SKIP() << "This test requires exactly 2 MPI ranks";
  void *lmp = open_partitioned_lammps(2, "2x1");
  ASSERT_NE(lmp, nullptr);

  setup_partition_zero_pair_system(lmp);
  lammps_command(lmp, "variable beadid world 01 02");
  lammps_command(lmp, "fix cp all pimd/langevin method nmpimd ensemble nvt integrator baoab "
                     "temp 1.0 thermostat PILE_L 1234 tau 1.0 fixcom no");
  lammps_command(lmp, "run 2 post no");
  lammps_command(lmp, "write_restart pimd_langevin_nmpimd_nvt.restart.${beadid}");
  lammps_command(lmp, "clear");
  lammps_command(lmp, "read_restart pimd_langevin_nmpimd_nvt.restart.${beadid}");
  lammps_command(lmp, "fix cp all pimd/langevin method nmpimd ensemble nvt integrator baoab "
                     "temp 1.0 thermostat PILE_L 1234 tau 1.0 fixcom no");
  lammps_command(lmp, "run 1 post no");

  EXPECT_EQ(fix_vector_size(lmp, "cp"), pimd_langevin_test::kNuclearPrefixScalars);
  expect_all_finite(lmp, "cp");

  lammps_close(lmp);
  cleanup_partition_restart("pimd_langevin_nmpimd_nvt.restart", 2);
}

TEST(FixPIMDLangevinRestartTest, PIMDNVTRestartContinuation)
{
  if (!has_exact_procs(2)) GTEST_SKIP() << "This test requires exactly 2 MPI ranks";
  void *lmp = open_partitioned_lammps(2, "2x1");
  ASSERT_NE(lmp, nullptr);

  setup_partition_zero_pair_system(lmp);
  lammps_command(lmp, "variable beadid world 01 02");
  lammps_command(lmp, "fix cp all pimd/langevin method pimd ensemble nvt integrator baoab "
                     "temp 1.0 thermostat PILE_L 1234 tau 1.0 fixcom no");
  lammps_command(lmp, "run 2 post no");
  lammps_command(lmp, "write_restart pimd_langevin_pimd_nvt.restart.${beadid}");
  lammps_command(lmp, "clear");
  lammps_command(lmp, "read_restart pimd_langevin_pimd_nvt.restart.${beadid}");
  lammps_command(lmp, "fix cp all pimd/langevin method pimd ensemble nvt integrator baoab "
                     "temp 1.0 thermostat PILE_L 1234 tau 1.0 fixcom no");
  lammps_command(lmp, "run 1 post no");

  EXPECT_EQ(fix_vector_size(lmp, "cp"), pimd_langevin_test::kNuclearPrefixScalars);
  expect_all_finite(lmp, "cp");

  lammps_close(lmp);
  cleanup_partition_restart("pimd_langevin_pimd_nvt.restart", 2);
}

TEST(FixPIMDLangevinRestartTest, NMPIMDNPTBZPRestartContinuation)
{
  const auto iso = pimd_langevin_test::iso_barostat_indices();
  if (!has_exact_procs(4)) GTEST_SKIP() << "This test requires exactly 4 MPI ranks";

  void *lmp = open_partitioned_lammps(4, "4x1");
  ASSERT_NE(lmp, nullptr);

  setup_barostat_system(lmp);
  lammps_command(lmp, "variable beadid world 01 02 03 04");
  lammps_command(lmp, "fix cp all pimd/langevin method nmpimd ensemble npt integrator obabo "
                     "temp 0.8 thermostat PILE_L 1234 tau 1.0 iso 1.0 barostat BZP taup 1.0 "
                     "fixcom no");
  lammps_command(lmp, "run 2 post no");
  const double vw_before = fix_value(lmp, "cp", iso.vw0);
  lammps_command(lmp, "write_restart pimd_langevin_npt_bzp.restart.${beadid}");
  lammps_command(lmp, "clear");
  lammps_command(lmp, "read_restart pimd_langevin_npt_bzp.restart.${beadid}");
  lammps_command(lmp, "fix cp all pimd/langevin method nmpimd ensemble npt integrator obabo "
                     "temp 0.8 thermostat PILE_L 1234 tau 1.0 iso 1.0 barostat BZP taup 1.0 "
                     "fixcom no");
  lammps_command(lmp, "run 0 post no");

  EXPECT_EQ(fix_vector_size(lmp, "cp"), pimd_langevin_test::kIsoBarostatVectorSize);
  EXPECT_NEAR(fix_value(lmp, "cp", iso.vw0), vw_before, 1.0e-12);
  expect_all_finite(lmp, "cp");

  lammps_close(lmp);
  cleanup_partition_restart("pimd_langevin_npt_bzp.restart", 4);
}

TEST(FixPIMDLangevinExpectedErrorTest, PIMDRejectsMultipleRanksPerBead)
{
  if (!has_exact_procs(4)) GTEST_SKIP() << "This test requires exactly 4 MPI ranks";
  void *lmp = open_partitioned_lammps(4, "2x2");
  ASSERT_NE(lmp, nullptr);

  setup_partition_lj_system(lmp);
  lammps_command(lmp, "fix cp all pimd/langevin method pimd ensemble nvt integrator obabo "
                     "temp 0.8 thermostat PILE_L 1234 tau 1.0 fixcom no");
  lammps_command(lmp, "run 0 post no");

  EXPECT_EQ(lammps_has_error(lmp), 1);
  EXPECT_THAT(pimd_langevin_test::last_error(lmp),
              ContainsRegex(".*Method pimd only supports a single processor per bead.*"));
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
