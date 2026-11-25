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
#include "modify.h"

#define protected public
#include "../../src/REPLICA/fix_pimd_langevin.h"
#undef protected

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

LAMMPS *as_lammps(void *lmp)
{
  return static_cast<LAMMPS *>(lmp);
}

bool check_langevin_state_collective(void *lmp, const char *stage)
{
  auto *fix = dynamic_cast<FixPIMDLangevin *>(as_lammps(lmp)->modify->get_fix_by_id("cp"));
  int local_ok = fix != nullptr;
  if (fix) {
    for (double value : fix_vector(lmp, "cp"))
      if (!std::isfinite(value)) local_ok = 0;
    // Retain the regression check for initialized normal-mode coefficients.
    for (double *coefficients : {fix->tau_k, fix->c1_k, fix->c2_k}) {
      if (!coefficients)
        local_ok = 0;
      else
        for (int mode = 0; mode < fix->np; ++mode)
          if (!std::isfinite(coefficients[mode])) local_ok = 0;
    }
  }

  int global_ok;
  MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (!global_ok) ADD_FAILURE() << "Invalid Langevin state at " << stage;
  return global_ok != 0;
}

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
  char suffix[16];
  for (int i = 1; i <= count; ++i) {
    snprintf(suffix, sizeof(suffix), ".%02d", i);
    platform::unlink((prefix + suffix).c_str());
  }
}

struct PartitionedLammpsHandle {
  void *lmp = nullptr;

  ~PartitionedLammpsHandle()
  {
    if (lmp) {
      lammps_close(lmp);
    }
  }
};

struct RestartFileGuard {
  std::string prefix;
  int count = 0;

  RestartFileGuard(std::string prefix_in, int count_in) :
      prefix(std::move(prefix_in)), count(count_in)
  {
    cleanup_partition_restart(prefix, count);
  }

  ~RestartFileGuard() { cleanup_partition_restart(prefix, count); }
};

void *open_partitioned_lammps(const char *partition)
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

}    // namespace

class FixPIMDLangevinSerialTest : public LAMMPSTest {
 protected:
  int fix_vector_size(const char *id) { return pimd_langevin_test::fix_vector_size(lmp, id); }
  double fix_value(const char *id, int index) { return pimd_langevin_test::fix_value(lmp, id, index); }
};

TEST_F(FixPIMDLangevinSerialTest, InputValidation)
{

  setup_single_atom_zero_pair([this](const std::string &line) { command(line); });
  command("fix cp all pimd/langevin temp 2.0 thermostat PILE_L 1234 fmass 0.5 "
          "ensemble nve sp 1.5 tau 0.7 fmmode physical fixcom no integrator obabo "
          "method nmpimd");
  auto *fix = dynamic_cast<FixPIMDLangevin *>(lmp->modify->get_fix_by_id("cp"));
  ASSERT_NE(fix, nullptr);
  EXPECT_DOUBLE_EQ(fix->temp, 2.0);
  EXPECT_DOUBLE_EQ(fix->fmass, 0.5);
  EXPECT_DOUBLE_EQ(fix->sp, 1.5);
  EXPECT_DOUBLE_EQ(fix->tau, 0.7);
  EXPECT_EQ(fix->removecomflag, 0);
  command("run 2 post no");
  EXPECT_TRUE(std::isfinite(fix_value("cp", pimd_langevin_test::TOTAL_ENERGY)));

  command("clear");

  for (const char *options : {"method cmd", "removecom yes", "temp",
                              "thermostat PILE_L", "unknown 1"}) {
    setup_single_atom_zero_pair([this](const std::string &line) { command(line); });
    // universe_all() tears down output state; isolate these errors in a child process.
    EXPECT_EXIT({
      try {
        command(std::string("fix cp all pimd/langevin ensemble nve ") + options);
      } catch (const LAMMPSException &) {
        std::_Exit(0);
      }
      std::_Exit(1);
    }, ::testing::ExitedWithCode(0), "") << options;
    command("clear");
  }

  for (const char *options : {"esynch yes temp 1.0", "temp 1.0 esynch no"}) {
    setup_single_atom_zero_pair([this](const std::string &line) { command(line); });
    command(std::string("fix cp all pimd/langevin/bosonic method pimd ensemble nve ") + options);
    command("run 2 post no");
    EXPECT_EQ(fix_vector_size("cp"), 6);
    for (int i = 0; i < 6; ++i) EXPECT_TRUE(std::isfinite(fix_value("cp", i)));
    command("clear");
  }

  setup_single_atom_zero_pair([this](const std::string &line) { command(line); });
  EXPECT_EXIT({
    try {
      command("fix cp all pimd/langevin/bosonic method pimd ensemble nve esynch invalid");
    } catch (const LAMMPSException &e) {
      std::_Exit(std::string(e.what()).find("esynch parameter can only receive yes or no") !=
                         std::string::npos ? 0 : 1);
    }
    std::_Exit(1);
  }, ::testing::ExitedWithCode(0), "");
}

TEST_F(FixPIMDLangevinSerialTest, SingleAtomEstimators)
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

TEST(FixPIMDLangevinMPI, NormalModeIntegration)
{
  if (!has_exact_procs(2)) GTEST_SKIP() << "This test requires exactly 2 MPI ranks";
  for (bool thermostat : {false, true}) {
    PartitionedLammpsHandle handle;
    handle.lmp = open_partitioned_lammps("2x1");
    ASSERT_NE(handle.lmp, nullptr);
    setup_partition_zero_pair_system(handle.lmp);
    const std::string settings = thermostat ? "nvt integrator baoab" : "nve integrator obabo";
    lammps_command(handle.lmp,
                   ("fix cp all pimd/langevin method nmpimd ensemble " + settings +
                    " temp 1.0 thermostat PILE_L 1234 tau 1.0 fixcom no").c_str());
    lammps_command(handle.lmp, "run 0 post no");
    const double before = fix_value(handle.lmp, "cp", pimd_langevin_test::TOTAL_ENERGY);
    lammps_command(handle.lmp, "run 10 post no");
    EXPECT_EQ(fix_vector_size(handle.lmp, "cp"), pimd_langevin_test::kNuclearPrefixScalars);
    expect_all_finite(handle.lmp, "cp");
    if (thermostat)
      ASSERT_TRUE(check_langevin_state_collective(handle.lmp, "after_run"));
    else
      EXPECT_NEAR(fix_value(handle.lmp, "cp", pimd_langevin_test::TOTAL_ENERGY), before, 1.0e-8);
  }
}

TEST(FixPIMDLangevinMPI, NMPIMDNVTOBABOMultiRankPerBead)
{
  if (!has_exact_procs(4)) GTEST_SKIP() << "This test requires exactly 4 MPI ranks";
  PartitionedLammpsHandle handle;
  handle.lmp = open_partitioned_lammps("2x2");
  ASSERT_NE(handle.lmp, nullptr);

  setup_partition_lj_system(handle.lmp);
  lammps_command(handle.lmp, "fix cp all pimd/langevin method nmpimd ensemble nvt integrator obabo "
                     "temp 0.8 sp 0.2 thermostat PILE_L 1234 tau 1.0 fixcom no");
  lammps_command(handle.lmp, "run 2 post no");

  EXPECT_EQ(fix_vector_size(handle.lmp, "cp"), pimd_langevin_test::kNuclearPrefixScalars);
  expect_all_finite(handle.lmp, "cp");
  EXPECT_GT(lammps_get_natoms(handle.lmp), 0.0);
}

TEST(FixPIMDLangevinMPI, CartesianIntegration)
{
  if (!has_exact_procs(2)) GTEST_SKIP() << "This test requires exactly 2 MPI ranks";
  for (const char *settings : {"nve integrator obabo", "nvt integrator baoab"}) {
    PartitionedLammpsHandle handle;
    handle.lmp = open_partitioned_lammps("2x1");
    ASSERT_NE(handle.lmp, nullptr);
    setup_partition_zero_pair_system(handle.lmp);
    lammps_command(handle.lmp,
                   (std::string("fix cp all pimd/langevin method pimd ensemble ") + settings +
                    " temp 1.0 thermostat PILE_L 1234 tau 1.0 fixcom no").c_str());
    lammps_command(handle.lmp, "run 2 post no");
    EXPECT_EQ(fix_vector_size(handle.lmp, "cp"), pimd_langevin_test::kNuclearPrefixScalars);
    expect_all_finite(handle.lmp, "cp");
    EXPECT_GT(std::abs(fix_value(handle.lmp, "cp", pimd_langevin_test::SE_BEAD)), 0.0);
  }
}

TEST(FixPIMDLangevinMPI, IsotropicBarostats)
{
  if (!has_exact_procs(4)) GTEST_SKIP() << "This test requires exactly 4 MPI ranks";
  for (const char *barostat : {"BZP", "MTTK"}) {
    PartitionedLammpsHandle handle;
    handle.lmp = open_partitioned_lammps("4x1");
    ASSERT_NE(handle.lmp, nullptr);
    setup_barostat_system(handle.lmp);
    lammps_command(handle.lmp,
                   (std::string("fix cp all pimd/langevin method nmpimd ensemble npt integrator ") +
                    "obabo temp 0.8 thermostat PILE_L 1234 tau 1.0 iso 1.0 barostat " +
                    barostat + " taup 1.0 fixcom no").c_str());
    lammps_command(handle.lmp, "run 2 post no");
    EXPECT_EQ(fix_vector_size(handle.lmp, "cp"), pimd_langevin_test::kIsoBarostatVectorSize);
    expect_all_finite(handle.lmp, "cp");
  }
}

TEST(FixPIMDLangevinMPI, AnisotropicBarostatStateMatchesAcrossPartitions)
{
  if (!has_exact_procs(4)) GTEST_SKIP() << "This test requires exactly 4 MPI ranks";
  for (const char *partition : {"4x1", "2x2"}) {
    for (const char *integrator : {"obabo", "baoab"}) {
      for (const char *ensemble : {"npt", "nph"}) {
        SCOPED_TRACE(std::string(partition) + " " + integrator + " " + ensemble);
        PartitionedLammpsHandle handle;
        handle.lmp = open_partitioned_lammps(partition);
        ASSERT_NE(handle.lmp, nullptr);
        setup_barostat_system(handle.lmp);
        lammps_command(handle.lmp,
                       (std::string("fix cp all pimd/langevin method nmpimd ensemble ") +
                        ensemble + " integrator " + integrator +
                        " temp 0.8 thermostat PILE_L 1234 tau 1.0 aniso 1.0 "
                        "barostat BZP taup 1.0 fixcom no").c_str());
        lammps_command(handle.lmp, "run 3 post no");

        EXPECT_EQ(fix_vector_size(handle.lmp, "cp"), pimd_langevin_test::kAnisoBarostatVectorSize);
        // Compare the public barostat diagnostics and the serialized state.
        // Nuclear bead energies are intentionally not compared across modes.
        double local[11];
        for (int i = 0; i < 3; ++i) local[i] = fix_value(handle.lmp, "cp", 10 + i);
        local[3] = fix_value(handle.lmp, "cp", 13);
        local[4] = fix_value(handle.lmp, "cp", 16);
        auto *fix = static_cast<FixPIMDLangevin *>(
            as_lammps(handle.lmp)->modify->get_fix_by_id("cp"));
        EXPECT_EQ(fix->pack_base_restart(local + 5), 6);
        double minimum[11], maximum[11];
        MPI_Allreduce(local, minimum, 11, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(local, maximum, 11, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        for (int i = 0; i < 11; ++i) {
          EXPECT_TRUE(std::isfinite(local[i]));
          EXPECT_NEAR(minimum[i], maximum[i], 1.0e-12) << "field " << i;
        }
      }
    }
  }
}

TEST(FixPIMDLangevinRestartTest, NVTStateRestored)
{
  if (!has_exact_procs(2)) GTEST_SKIP() << "This test requires exactly 2 MPI ranks";
  for (const char *method : {"nmpimd", "pimd"}) {
    PartitionedLammpsHandle handle;
    handle.lmp = open_partitioned_lammps("2x1");
    ASSERT_NE(handle.lmp, nullptr);
    RestartFileGuard restart_guard{std::string("langevin_nvt_") + method + ".restart", 2};
    setup_partition_zero_pair_system(handle.lmp);
    const std::string fix_command = std::string("fix cp all pimd/langevin method ") + method +
        " ensemble nvt integrator baoab temp 1.0 thermostat PILE_L 1234 tau 1.0 fixcom no";
    lammps_command(handle.lmp, "variable beadid world 01 02");
    lammps_command(handle.lmp, fix_command.c_str());
    lammps_command(handle.lmp, "run 2 post no");
    lammps_command(handle.lmp, fmt::format("write_restart {}.${{beadid}}", restart_guard.prefix).c_str());
    lammps_command(handle.lmp, "clear");
    lammps_command(handle.lmp, fmt::format("read_restart {}.${{beadid}}", restart_guard.prefix).c_str());
    lammps_command(handle.lmp, fix_command.c_str());
    lammps_command(handle.lmp, "run 1 post no");
    EXPECT_EQ(fix_vector_size(handle.lmp, "cp"), pimd_langevin_test::kNuclearPrefixScalars);
    expect_all_finite(handle.lmp, "cp");
    if (std::string(method) == "nmpimd")
      ASSERT_TRUE(check_langevin_state_collective(handle.lmp, "after_restart"));
  }
}

TEST(FixPIMDLangevinRestartTest, NMPIMDNPTBZPRestartContinuation)
{
  const auto iso = pimd_langevin_test::iso_barostat_indices();
  if (!has_exact_procs(4)) GTEST_SKIP() << "This test requires exactly 4 MPI ranks";

  PartitionedLammpsHandle handle;
  handle.lmp = open_partitioned_lammps("4x1");
  ASSERT_NE(handle.lmp, nullptr);
  RestartFileGuard restart_guard{"pimd_langevin_npt_bzp.restart", 4};

  setup_barostat_system(handle.lmp);
  lammps_command(handle.lmp, "variable beadid world 01 02 03 04");
  lammps_command(handle.lmp, "fix cp all pimd/langevin method nmpimd ensemble npt integrator obabo "
                     "temp 0.8 thermostat PILE_L 1234 tau 1.0 iso 1.0 barostat BZP taup 1.0 "
                     "fixcom no");
  lammps_command(handle.lmp, "run 2 post no");
  const double vw_before = fix_value(handle.lmp, "cp", iso.vw0);
  lammps_command(handle.lmp, fmt::format("write_restart {}.${{beadid}}", restart_guard.prefix).c_str());
  lammps_command(handle.lmp, "clear");
  lammps_command(handle.lmp, fmt::format("read_restart {}.${{beadid}}", restart_guard.prefix).c_str());
  lammps_command(handle.lmp, "fix cp all pimd/langevin method nmpimd ensemble npt integrator obabo "
                     "temp 0.8 thermostat PILE_L 1234 tau 1.0 iso 1.0 barostat BZP taup 1.0 "
                     "fixcom no");
  lammps_command(handle.lmp, "run 0 post no");

  EXPECT_EQ(fix_vector_size(handle.lmp, "cp"), pimd_langevin_test::kIsoBarostatVectorSize);
  EXPECT_NEAR(fix_value(handle.lmp, "cp", iso.vw0), vw_before, 1.0e-12);
  expect_all_finite(handle.lmp, "cp");
}

TEST(FixPIMDLangevinExpectedErrorTest, PIMDRejectsMultipleRanksPerBead)
{
  if (!has_exact_procs(4)) GTEST_SKIP() << "This test requires exactly 4 MPI ranks";
  void *lmp = open_partitioned_lammps("2x2");
  ASSERT_NE(lmp, nullptr);

  setup_partition_lj_system(lmp);
  lammps_command(lmp, "fix cp all pimd/langevin method pimd ensemble nvt integrator obabo "
                     "temp 0.8 thermostat PILE_L 1234 tau 1.0 fixcom no");
  lammps_command(lmp, "run 0 post no");

  EXPECT_EQ(lammps_has_error(lmp), 1);
  EXPECT_THAT(pimd_langevin_test::last_error(lmp),
              ContainsRegex(".*Method pimd only supports a single processor per bead.*"));
}

TEST(FixPIMDLangevinBosonicTest, CartesianPIMDSmoke)
{
  int nprocs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (nprocs != 4) GTEST_SKIP() << "This test requires exactly 4 MPI ranks";

  const char *args[] = {"LAMMPS_test", "-log", "none", "-partition", "4x1", "-echo",
                        "screen",      "-nocite",       "-in",        "none", nullptr};
  char **argv = (char **) args;
  const int argc = (sizeof(args) / sizeof(char *)) - 1;

  void *lmp = nullptr;
  EXPECT_NO_THROW(lmp = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr));
  ASSERT_NE(lmp, nullptr);

  auto command = [lmp](const char *line) { lammps_command(lmp, line); };

  command("units electron");
  command("dimension 3");
  command("boundary p p p");
  command("atom_style atomic");
  command("atom_modify map yes");
  command("pair_style none");
  command("timestep 0.5");
  command("variable Temp equal 17.4");
  command("variable k equal 1.2154614120000001e-08");
  command("variable ibead world 01 02 03 04");
  command("region box block -1500 1500 -1500 1500 -1500 1500");
  command("create_box 1 box");
  command("create_atoms 1 single -15.0 10.0 1.0");
  command("create_atoms 1 single -14.5 9.5 1.5");
  command("create_atoms 1 single -14.0 9.0 2.0");
  command("mass 1 0.00054858");
  command("velocity all create ${Temp} 18889${ibead} mom yes rot yes dist gaussian");
  command("fix harm all spring/self ${k}");
  command("fix_modify harm energy yes");
  command("fix cp all pimd/langevin/bosonic ensemble nvt temp ${Temp} thermostat PILE_L "
          "12345 tau 50 fixcom no");
  command("run 5 post no");

  EXPECT_EQ(pimd_langevin_test::fix_vector_size(lmp, "cp"), pimd_langevin_test::kBosonicVectorSize);
  for (int i = 0; i < pimd_langevin_test::kBosonicVectorSize; ++i) {
    EXPECT_TRUE(std::isfinite(pimd_langevin_test::fix_value(lmp, "cp", i))) << "index=" << i;
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
