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
#include "universe.h"

#include "../testing/core.h"
#include "gtest/gtest.h"
#include "pimd_test_utils.h"

#include <cmath>
#include <mpi.h>
#include <vector>

bool verbose = false;

namespace LAMMPS_NS {

class FixPIMDNVTSerialTest : public LAMMPSTest {
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

  void setup_nuclear_fix(const char *style, const char *method = "nmpimd")
  {
    command(std::string("fix cp all ") + style +
            " method " + method +
            " thermostat NHC temp 1.0 Tdamp 0.5 "
            "tchain 3 tloop 1");
  }

  double fix_value(const char *id, int index)
  {
    return pimd_test::fix_value(lmp, id, index);
  }
};

TEST_F(FixPIMDNVTSerialTest, NMPIMDStyleParsesAndRuns)
{
  setup_zero_pair_system();
  setup_nuclear_fix("pimd/nvt");
  command("run 0 post no");
  for (int i = 0; i < pimd_test::nuclear_vector_size(); ++i) {
    EXPECT_TRUE(std::isfinite(fix_value("cp", i))) << "index=" << i;
  }
}

TEST_F(FixPIMDNVTSerialTest, InternalThermostatIsNotRegisteredAsAnIntegrator)
{
  setup_zero_pair_system();
  const auto nfix = lmp->modify->get_fix_list().size();
  const auto ncompute = lmp->modify->get_compute_list().size();
  for (int repeat = 0; repeat < 2; ++repeat) {
    setup_nuclear_fix("pimd/nvt");
    EXPECT_EQ(lmp->modify->get_fix_list().size(), nfix + 1);
    EXPECT_EQ(lmp->modify->get_fix_by_id("cp_nhc"), nullptr);
    command("run 2 post no");
    command("unfix cp");
    EXPECT_EQ(lmp->modify->get_fix_list().size(), nfix);
    EXPECT_EQ(lmp->modify->get_compute_list().size(), ncompute);
  }
}

TEST_F(FixPIMDNVTSerialTest, PIMDStyleParsesAndRuns)
{
  setup_zero_pair_system();
  setup_nuclear_fix("pimd/nvt", "pimd");
  command("run 20 post no");
  for (int i = 0; i < pimd_test::nuclear_vector_size(); ++i) {
    EXPECT_TRUE(std::isfinite(fix_value("cp", i))) << "index=" << i;
  }
}

TEST_F(FixPIMDNVTSerialTest, CMDStyleParsesAndRuns)
{
  setup_zero_pair_system();
  setup_nuclear_fix("pimd/nvt", "cmd");
  command("run 20 post no");
  for (int i = 0; i < pimd_test::nuclear_vector_size(); ++i) {
    EXPECT_TRUE(std::isfinite(fix_value("cp", i))) << "index=" << i;
  }
}

TEST_F(FixPIMDNVTSerialTest, BosonicNVTStyleParsesAndRuns)
{
  setup_zero_pair_system();
  command("fix cp all pimd/nvt/bosonic method pimd temp 1.0 nhc 2");
  command("run 20 post no");
  for (int i = 0; i < pimd_test::nuclear_vector_size(2); ++i) {
    EXPECT_TRUE(std::isfinite(fix_value("cp", i))) << "index=" << i;
  }
}

TEST_F(FixPIMDNVTSerialTest, BosonicNVTRejectsNMPIMD)
{
  setup_zero_pair_system();
  EXPECT_ANY_THROW(command("fix cp all pimd/nvt/bosonic method nmpimd temp 1.0 nhc 2"));
}

TEST_F(FixPIMDNVTSerialTest, RejectsEnsembleKeyword)
{
  setup_zero_pair_system();
  EXPECT_ANY_THROW(command("fix cp all pimd/nvt method nmpimd ensemble nvt thermostat NHC "
                           "temp 1.0 Tdamp 0.5 tchain 3 tloop 1"));
}

TEST_F(FixPIMDNVTSerialTest, ChainArgumentsMatchFixNHBounds)
{
  for (const char *style : {"pimd/nvt", "pimd/nvt/bosonic", "pimd/uvt"}) {
    for (const char *options : {"tchain 0", "tchain -1", "tloop -1"}) {
      setup_zero_pair_system();
      try {
        command(std::string("fix cp all ") + style + " Tdamp 0.5 " + options);
        FAIL() << "Accepted " << options << " for " << style;
      } catch (const LAMMPSException &e) {
        EXPECT_NE(std::string(e.what()).find(std::string("Invalid fix ") + style),
                  std::string::npos);
      }
      command("clear");
    }
  }

  setup_zero_pair_system();
  EXPECT_NO_THROW(command("fix cp all pimd/nvt Tdamp 0.5 tchain 1 tloop 0"));
}

TEST_F(FixPIMDNVTSerialTest, RejectsUnknownMethods)
{
  const char *unsupported_methods[] = {"nmrpmd", "tprpmd", "tp-rpmd"};
  for (const char *method : unsupported_methods) {
    setup_zero_pair_system();
    EXPECT_ANY_THROW(command(std::string("fix cp all pimd/nvt method ") + method +
                             " thermostat NHC temp 1.0 Tdamp 0.5 tchain 3 tloop 1"))
        << "method=" << method;
    command("clear");
  }
}

TEST_F(FixPIMDNVTSerialTest, P1StandaloneRunProducesNuclearState)
{
  setup_zero_pair_system();
  setup_nuclear_fix("pimd/nvt", "nmpimd");
  command("run 20 post no");

  for (int i = 0; i < pimd_test::nuclear_vector_size(); ++i) {
    EXPECT_TRUE(std::isfinite(fix_value("cp", i))) << "index=" << i;
  }
  EXPECT_GT(fix_value("cp", 0), 0.0);
  EXPECT_TRUE(std::isfinite(fix_value("cp", 13)));
}

TEST_F(FixPIMDNVTSerialTest, ThermostatAdvancesFullElapsedTime)
{
  for (const std::string style : {"pimd/nvt", "pimd/uvt"})
    for (const std::string split : {"obabo", "baoab"})
      for (const std::string loops : {"1", "3"}) {
        SCOPED_TRACE(style + " " + split + " tloop=" + loops);
        command("clear");
        setup_zero_pair_system();
        command("velocity all set 0.0 0.0 0.0");
        command("timestep 0.001");
        std::string extra;
        if (style == "pimd/uvt") {
          command("variable dEdN equal 0.0");
          extra = " mu 0.0 Udamp 1.0 ne 1.0 ne_velocity 0.0 dedn v_dEdN";
        }
        command("fix cp all " + style + " method nmpimd temp 1.0 Tdamp 1.0 "
                "tchain 1 removecom no integrator " + split + " tloop " + loops + extra);
        command("run 7 post no");
        // With zero forces and velocities, K remains zero. For a single
        // thermostat and Tdamp=1, zeta_dot=-1, zeta=-t, eta=-t*t/2.
        const double elapsed = 0.007;
        EXPECT_NEAR(fix_value("cp", 10), -0.5 * elapsed * elapsed, 1.0e-13);
        EXPECT_NEAR(fix_value("cp", 11), -elapsed, 1.0e-13);
      }
}

TEST_F(FixPIMDNVTSerialTest, RestartRestoresNuclearState)
{
  std::vector<double> before(pimd_test::nuclear_vector_size());

  setup_zero_pair_system();
  setup_nuclear_fix("pimd/nvt", "nmpimd");
  command("run 20 post no");

  for (int i = 0; i < static_cast<int>(before.size()); ++i) before[i] = fix_value("cp", i);

  command("write_restart pimd_nvt.restart");
  command("clear");
  command("read_restart pimd_nvt.restart");
  setup_nuclear_fix("pimd/nvt", "nmpimd");
  command("run 0 post no");
  platform::unlink("pimd_nvt.restart");

  for (int i = 0; i < static_cast<int>(before.size()); ++i) {
    EXPECT_NEAR(fix_value("cp", i), before[i], 1.0e-10) << "index=" << i;
  }
}

TEST(FixPIMDNVTMPI, PartitionedRunExercisesBeadExpansion)
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
  command("fix cp all pimd/nvt method nmpimd thermostat NHC temp 0.8 "
          "Tdamp 0.2 tchain 3 tloop 1");
  command("run 1 post no");

  for (int i = 0; i < pimd_test::nuclear_vector_size(); ++i) {
    EXPECT_TRUE(std::isfinite(fix_value("cp", i))) << "index=" << i;
  }

  lammps_close(lmp);
}

TEST(FixPIMDNVTMPI, PIMDPartitionedRunSupportsOneRankPerBead)
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
  command("velocity all create 0.8 86420 mom yes rot no dist gaussian");
  command("fix cp all pimd/nvt method pimd thermostat NHC temp 0.8 "
          "Tdamp 0.2 tchain 3 tloop 1");
  command("run 2 post no");

  for (int i = 0; i < pimd_test::nuclear_vector_size(); ++i) {
    EXPECT_TRUE(std::isfinite(fix_value("cp", i))) << "index=" << i;
  }

  lammps_close(lmp);
}

TEST(FixPIMDNVTMPI, CMDLeavesCentroidChainInactive)
{
  int nprocs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (nprocs != 2) GTEST_SKIP() << "This test requires exactly 2 MPI ranks";

  for (const char *integrator : {"obabo", "baoab"}) {
    const char *args[] = {"LAMMPS_test", "-log", "none", "-partition", "2x1",
                          "-echo", "screen", "-nocite", "-in", "none"};
    void *handle = lammps_open(10, (char **) args, MPI_COMM_WORLD, nullptr);
    ASSERT_NE(handle, nullptr);
    auto command = [handle](const std::string &line) { lammps_command(handle, line.c_str()); };
    command("units lj");
    command("atom_style atomic");
    command("atom_modify map yes");
    command("region box block 0 4 0 4 0 4");
    command("create_box 1 box");
    command("create_atoms 1 single 1 1 1");
    command("create_atoms 1 single 2 2 2");
    command("mass 1 1.0");
    command("pair_style zero 2.5");
    command("pair_coeff * *");
    command("velocity all set 0.2 0.1 -0.05");
    command("timestep 0.001");
    command(std::string("fix cp all pimd/nvt method cmd temp 1.0 Tdamp 0.5 ") +
            "tchain 3 tloop 2 removecom no integrator " + integrator);
    command("run 5 post no");

    auto *lmp = static_cast<LAMMPS *>(handle);
    for (int i = 0; i < pimd_test::nuclear_vector_size(); ++i)
      EXPECT_TRUE(std::isfinite(pimd_test::fix_value(handle, "cp", i)));
    if (lmp->universe->iworld == 0) {
      for (int i = 10; i < 16; ++i)
        EXPECT_DOUBLE_EQ(pimd_test::fix_value(handle, "cp", i), 0.0);
    } else {
      EXPECT_NE(pimd_test::fix_value(handle, "cp", 13), 0.0);
    }
    lammps_close(handle);
  }
}

TEST(FixPIMDNVTMPI, MultiRankPerBeadRunProducesFiniteThermostatState)
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
  command("region box block 0 3 0 2 0 2");
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
  command("velocity all create 0.8 24680 mom yes rot no dist gaussian");
  command("fix cp all pimd/nvt method nmpimd thermostat NHC temp 0.8 "
          "Tdamp 0.2 tchain 3 tloop 1");
  command("run 2 post no");

  for (int i = 0; i < pimd_test::nuclear_vector_size(); ++i) {
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
