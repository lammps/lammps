/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with the U.S. Government retains
   certain rights in this software.
------------------------------------------------------------------------- */

#define LAMMPS_LIB_MPI 1

#include "atom.h"
#include "compute.h"
#include "fix.h"
#include "force.h"
#include "lammps.h"
#include "modify.h"
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

  void setup_uvt_fix(const char *style, double ne, double ne_velocity, double mu_value = 2.0,
                     double udamp = 0.5)
  {
    command("variable k_quad equal 5.0");
    command("variable N0_quad equal 1.0");
    command(pimd_test::quadratic_dedn_variable());
    command(std::string("fix cp all ") + style +
            " method nmpimd thermostat NHC temp 1.0 Tdamp 0.5 "
            "tchain 3 tloop 1 mu " +
            std::to_string(mu_value) + " Udamp " + std::to_string(udamp) + " ne " +
            std::to_string(ne) + " ne_velocity " +
            std::to_string(ne_velocity) + " dedn v_dEdN");
  }

  double fix_value(const char *id, int index)
  {
    return pimd_test::fix_value(lmp, id, index);
  }
};

TEST_F(FixPIMDUVTSerialTest, ComputeDerivativeStaysCurrentThroughVariable)
{
  std::vector<double> reference;
  for (const char *source : {"c_position", "v_derivative", "c_positions[1]", "v_derivatives[1]"}) {
    setup_quadratic_system();
    command("velocity all set 1.0 0.0 0.0");
    command("compute position all reduce sum x");
    command("compute positions all reduce sum x y");
    command("variable derivative equal c_position");
    command("variable derivatives vector c_positions");
    command(std::string("fix cp all pimd/uvt temp 1.0 Tdamp 0.5 mu 2.0 Udamp 0.5 ") +
            "ne 1.8 ne_velocity 0.0 removecom no dedn " + source);
    // Do not let thermo refresh the derivative's compute on intermediate steps.
    command("thermo 1000");
    command("thermo_style custom step");
    command("run 10 post no");

    double expected = 0.0;
    for (int i = 0; i < lmp->atom->nlocal; ++i) expected += lmp->atom->x[i][0];
    auto *fix = lmp->modify->get_fix_by_id("cp");
    int dim;
    const auto value = [&](const char *name) {
      return *static_cast<double *>(fix->extract(name, dim));
    };
    EXPECT_NEAR(value("dedn"), expected, 1.0e-10) << source;
    if (reference.empty()) {
      reference = {value("ne"), value("ne_dot")};
    } else {
      EXPECT_NEAR(value("ne"), reference[0], 1.0e-12) << source;
      EXPECT_NEAR(value("ne_dot"), reference[1], 1.0e-12) << source;
    }
    command("clear");
  }
}

TEST_F(FixPIMDUVTSerialTest, VariableDerivativeSchedulesPerAtomEnergy)
{
  setup_quadratic_system();
  command("pair_style lj/cut 2.5");
  command("pair_coeff * * 1.0 1.0");
  command("compute total all pe");
  command("compute energies all pe/atom");
  command("compute energy all reduce sum c_energies");
  command("variable derivative equal c_energy");
  command("fix cp all pimd/uvt temp 1.0 Tdamp 0.5 mu 2.0 Udamp 0.5 "
          "ne 1.8 ne_velocity 0.0 dedn v_derivative");
  command("thermo 1000");
  command("thermo_style custom step");
  command("run 10 post no");
  int dim;
  auto *fix = lmp->modify->get_fix_by_id("cp");
  const double derivative = *static_cast<double *>(fix->extract("dedn", dim));
  auto *total = lmp->modify->get_compute_by_id("total");
  EXPECT_NEAR(derivative, total->compute_scalar(), 1.0e-9);
}

TEST_F(FixPIMDUVTSerialTest, RestartRestoresElectronState)
{
  const auto uvt = pimd_test::uvt_vector_indices();
  double ne_before = 0.0;
  double nedot_before = 0.0;
  double dedn_before = 0.0;

  setup_quadratic_system();
  setup_uvt_fix("pimd/uvt", 1.8, 0.0);
  command("run 20 post no");

  ne_before = fix_value("cp", uvt.ne);
  nedot_before = fix_value("cp", uvt.ne_dot);
  dedn_before = fix_value("cp", uvt.dedn);

  command("write_restart pimd_uvt.restart");
  command("clear");
  command("read_restart pimd_uvt.restart");
  setup_uvt_fix("pimd/uvt", 1.8, 0.0);
  command("run 0 post no");
  platform::unlink("pimd_uvt.restart");

  EXPECT_NEAR(fix_value("cp", uvt.ne), ne_before, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", uvt.ne_dot), nedot_before, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", uvt.dedn), dedn_before, 1.0e-10);
}

TEST_F(FixPIMDUVTSerialTest, InputValidation)
{
  setup_quadratic_system();
  command("variable derivative equal 0.0");
  TEST_FAILURE("Missing mu keyword", command("fix cp all pimd/uvt temp 1 Tdamp 0.5 "
                                           "Udamp 0.5 ne 1 dedn v_derivative"););
  command("clear");
  setup_quadratic_system();
  command("variable derivative equal 0.0");
  EXPECT_ANY_THROW(command("fix cp all pimd/uvt temp 1 Tdamp 0.5 mu 2 2 0.5 "
                           "ne 1 dedn v_derivative"));
  for (const char *method : {"pimd", "cmd", "invalid"}) {
    EXPECT_ANY_THROW(command(std::string("fix cp all pimd/uvt method ") + method +
                             " temp 1 Tdamp 0.5 mu 0 Udamp 0.5 ne 1 dedn v_derivative"));
  }
  command("fix cp all pimd/uvt temp 1 Tdamp 0.5 mu 0 Udamp 0.5 ne 1 dedn v_derivative");
  command("run 0 post no");
  EXPECT_DOUBLE_EQ(fix_value("cp", pimd_test::uvt_vector_indices().ne), 1.0);
}

TEST_F(FixPIMDUVTSerialTest, ElectronicOutputLayoutOmitsTargetMu)
{
  for (int tchain : {1, 3}) {
    setup_quadratic_system();
    command("variable derivative equal 5.0*(f_cp[11]-1.0)");
    command("fix cp all pimd/uvt temp 1.0 Tdamp 0.5 tchain " + std::to_string(tchain) +
            " mu 2.0 Udamp 0.5 ne 1.8 ne_velocity 0.25 dedn v_derivative");
    command("run 0 post no");
    const auto indices = pimd_test::uvt_vector_indices();
    auto *fix = lmp->modify->get_fix_by_id("cp");
    int dim;
    const double mass = *static_cast<double *>(fix->extract("ne_mass", dim));
    EXPECT_EQ(fix->size_vector, 16 + 4 * tchain);
    EXPECT_DOUBLE_EQ(fix_value("cp", indices.ne), 1.8);
    EXPECT_NEAR(fix_value("cp", indices.ne_ke), 0.5 * mass * 0.25 * 0.25, 1.0e-12);
    EXPECT_NEAR(fix_value("cp", indices.electronic_potential), -3.6, 1.0e-12);
    const double nuclear_dof = 3.0 * 108 - 3.0;
    EXPECT_NEAR(fix_value("cp", 11), (nuclear_dof + mass * 0.25 * 0.25) /
                                      (nuclear_dof + 1.0), 1.0e-12);
    EXPECT_DOUBLE_EQ(fix_value("cp", 12), 0.25);
    EXPECT_NEAR(fix_value("cp", 13), 4.0, 1.0e-12);
    EXPECT_NEAR(fix_value("cp", 0), 0.5 * nuclear_dof, 1.0e-10);
    EXPECT_EQ(fix->get_thermo_colname(10), "f_cp:Ne");
    EXPECT_EQ(fix->get_thermo_colname(11), "f_cp:T_ins");
    EXPECT_EQ(fix->get_thermo_colname(0), "");
    for (int i = 0; i < fix->size_vector; ++i)
      EXPECT_EQ(fix->extlist[i], (i < 10 || i >= 16) ? 1 : 0);
    const char *states[] = {"ne", "ne_dot", "dedn"};
    const int state_indices[] = {indices.ne, indices.ne_dot, indices.dedn};
    for (int i = 0; i < 3; ++i) {
      auto *value = static_cast<double *>(fix->extract(states[i], dim));
      ASSERT_NE(value, nullptr);
      EXPECT_EQ(dim, 1);
      EXPECT_NEAR(*value, fix_value("cp", state_indices[i]), 1.0e-12);
    }
    command("clear");
  }
}

TEST_F(FixPIMDUVTSerialTest, P1PhysicalOutputsMatchClassicalUVT)
{
  if (!Info::has_package("EXTRA-FIX") || !Info::has_package("EXTRA-COMPUTE"))
    GTEST_SKIP() << "Classical UVT requires EXTRA-FIX and EXTRA-COMPUTE";
  std::vector<double> reference;
  for (const std::string style : {"uvt", "pimd/uvt"}) {
    setup_quadratic_system();
    command("variable derivative equal 4.0");
    command("fix cp all " + style +
            (style == "uvt" ? " temp 1.0 1.0 0.5" : " temp 1.0 Tdamp 0.5") +
            " mu 2.0 Udamp 0.5 ne 1.8 ne_velocity 0.25 dedn v_derivative");
    command("run 0 post no");
    const int offset = style == "uvt" ? 0 : 10;
    if (reference.empty()) {
      for (int i = 0; i < 6; ++i) reference.push_back(fix_value("cp", i));
    } else {
      for (int i = 0; i < 6; ++i)
        EXPECT_NEAR(fix_value("cp", offset + i), reference[i], 1.0e-12) << "field " << i;
    }
    if (style == "pimd/uvt") {
      command("run 5 post no");
      EXPECT_TRUE(std::isfinite(fix_value("cp", 16)));
      EXPECT_GT(std::fabs(fix_value("cp", 19)), 1.0e-12);
    }
    command("clear");
  }
}

TEST(FixPIMDUVTMPI, InitialVectorDerivativeUsesBeadCoordinates)
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
  command("compute positions all reduce sum x y");
  command("variable derivative vector c_positions");
  command("fix cp all pimd/uvt temp 0.8 Tdamp 0.2 mu 1.5 Udamp 0.2 "
          "ne 1.8 dedn v_derivative[1]");
  command("run 0 post no");

  auto *instance = static_cast<LAMMPS *>(lmp);
  double local = 0.0;
  for (int i = 0; i < instance->atom->nlocal; ++i) local += instance->atom->x[i][0];
  double total;
  MPI_Allreduce(&local, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  EXPECT_NEAR(fix_value("cp", uvt.dedn), local, 1.0e-10);
  int dim;
  auto *fix = instance->modify->get_fix_by_id("cp");
  EXPECT_NEAR(*static_cast<double *>(fix->extract("dedn", dim)), total / 2.0, 1.0e-10);

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
  command("fix cp all pimd/uvt method nmpimd thermostat NHC temp 0.8 Tdamp 0.2 "
          "tchain 3 tloop 1 mu 1.5 Udamp 0.2 ne 1.7 ne_velocity 0.0 dedn v_dEdN");
  command("run 0 post no");

  EXPECT_NEAR(fix_value("cp", uvt.ne), 1.7, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", uvt.ne_dot), 0.0, 1.0e-10);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  EXPECT_NEAR(fix_value("cp", uvt.dedn), rank == 0 ? 2.1 : 0.9, 1.0e-10);
  auto *fix = static_cast<LAMMPS *>(lmp)->modify->get_fix_by_id("cp");
  int dim;
  EXPECT_NEAR(*static_cast<double *>(fix->extract("dedn", dim)), 1.5, 1.0e-10);

  command("run 1 post no");
  EXPECT_NEAR(fix_value("cp", uvt.ne), 1.7, 1.0e-10);
  EXPECT_NEAR(fix_value("cp", uvt.ne_dot), 0.0, 1.0e-10);
  lammps_close(lmp);
}

TEST(FixPIMDUVTMPI, MultiRankPerBeadAveragesDerivativeOncePerBead)
{
  int nprocs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (nprocs != 4) GTEST_SKIP() << "This test requires exactly 4 MPI ranks";

  const auto uvt = pimd_test::uvt_vector_indices();

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
  command("variable bead_dedn world 1.0 2.0");
  command("variable dEdN equal ${bead_dedn}");
  command("fix cp all pimd/uvt method nmpimd thermostat NHC temp 0.8 Tdamp 0.2 "
          "tchain 3 tloop 1 mu 1.5 Udamp 0.2 ne 1.7 ne_velocity 0.25 dedn v_dEdN removecom no");
  command("run 0 post no");

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  EXPECT_NEAR(fix_value("cp", uvt.dedn), rank < 2 ? 1.0 : 2.0, 1.0e-12);
  auto *instance = static_cast<LAMMPS *>(lmp);
  auto *fix = instance->modify->get_fix_by_id("cp");
  int dim;
  EXPECT_NEAR(*static_cast<double *>(fix->extract("dedn", dim)), 1.5, 1.0e-12);
  const double mass = *static_cast<double *>(fix->extract("ne_mass", dim));
  double local_ke = 0.0;
  for (int i = 0; i < instance->atom->nlocal; ++i)
    for (int j = 0; j < 3; ++j)
      local_ke += instance->atom->v[i][j] * instance->atom->v[i][j];
  double total_ke;
  MPI_Allreduce(&local_ke, &total_ke, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  const double expected_temp = (total_ke * instance->force->mvv2e / 2.0 +
                                mass * 0.25 * 0.25) /
      ((2.0 * 3.0 * instance->atom->natoms + 1.0) * instance->force->boltz);
  EXPECT_NEAR(fix_value("cp", 11), expected_temp, 1.0e-12);
  command("run 1 post no");
  EXPECT_NEAR(fix_value("cp", uvt.dedn), rank < 2 ? 1.0 : 2.0, 1.0e-12);

  lammps_close(lmp);
}

TEST(FixPIMDUVTMPI, SharedElectronicMassAndEnergyAreConsistentAcrossBeads)
{
  int nprocs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (nprocs != 2 && nprocs != 4)
    GTEST_SKIP() << "This test requires 2 or 4 MPI ranks";

  const std::vector<std::string> layouts =
      nprocs == 2 ? std::vector<std::string>{"2x1"} :
                    std::vector<std::string>{"2x2", "4x1"};
  for (const auto &layout : layouts) {
    for (const std::string removecom : {"yes", "no"}) {
      for (const std::string split : {"obabo", "baoab"}) {
        SCOPED_TRACE(layout + " removecom=" + removecom +
                     " integrator=" + split);
        const char *args[] = {"LAMMPS_test", "-log", "none", "-partition", layout.c_str(),
                              "-screen", "none", "-nocite", "-in", "none"};
        void *lmp = lammps_open(sizeof(args) / sizeof(char *), (char **) args,
                                MPI_COMM_WORLD, nullptr);
        ASSERT_NE(lmp, nullptr);
        auto command = [lmp](const std::string &line) { lammps_command(lmp, line.c_str()); };
        command("units lj");
        command("atom_style atomic");
        command("atom_modify map yes");
        command("region box block 0 10 0 10 0 10");
        command("create_box 1 box");
        command("create_atoms 1 single 1 1 1");
        command("create_atoms 1 single 2 2 2");
        command("mass 1 1.0");
        command("pair_style zero 2.0");
        command("pair_coeff * *");
        command("velocity all create 1.0 97531 mom yes rot no dist gaussian");
        command("timestep 0.002");
        command("variable response equal 0.0");
        command("fix cp all pimd/uvt temp 1.0 Tdamp 0.5 mu 0.0 Udamp 0.5 "
                "ne 1.0 ne_velocity 1.0 dedn v_response integrator " + split +
                " removecom " + removecom);
        command("thermo_style custom step");
        command("thermo 1000");
        command("run 0 post no");

        // The global mass uses the centroid's nuclear DOF: three for two
        // atoms with COM removal, otherwise six. Internal-mode DOF differ.
        const double expected_mass = removecom == "no" ? 1.5 : 0.75;
        auto check_shared_state = [&]() {
          auto *fix = static_cast<LAMMPS *>(lmp)->modify->get_fix_by_id("cp");
          int dim;
          double state[] = {*static_cast<double *>(fix->extract("ne_mass", dim)),
                             *static_cast<double *>(fix->extract("ne", dim)),
                             *static_cast<double *>(fix->extract("ne_dot", dim)),
                             pimd_test::fix_value(lmp, "cp", 14)};
          double minimum[4], maximum[4];
          MPI_Allreduce(state, minimum, 4, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
          MPI_Allreduce(state, maximum, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
          for (int i = 0; i < 4; ++i)
            EXPECT_NEAR(minimum[i], maximum[i], 1.0e-12) << "state index=" << i;
          EXPECT_DOUBLE_EQ(state[0], expected_mass);
          EXPECT_NEAR(state[3], 0.5 * expected_mass * state[2] * state[2], 1.0e-12);
        };
        {
          SCOPED_TRACE("after setup");
          check_shared_state();
        }
        command("run 1 post no");
        {
          SCOPED_TRACE("after integration");
          check_shared_state();
        }
        lammps_close(lmp);
      }
    }
  }
}

TEST(FixPIMDUVTMPI, ElectronicVelocityUsesMeanChainFriction)
{
  int nprocs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (nprocs != 2 && nprocs != 4)
    GTEST_SKIP() << "This test requires 2 or 4 MPI ranks";

  for (const char *split : {"obabo", "baoab"}) {
    const char *args[] = {"LAMMPS_test", "-log", "none", "-partition",
                         nprocs == 2 ? "2x1" : "2x2", "-screen", "none",
                         "-nocite", "-in", "none"};
    void *lmp = lammps_open(sizeof(args) / sizeof(char *), (char **) args,
                            MPI_COMM_WORLD, nullptr);
    ASSERT_NE(lmp, nullptr);
    auto command = [lmp](const char *line) { lammps_command(lmp, line); };
    command("units lj");
    command("atom_style atomic");
    command("atom_modify map yes");
    command("lattice sc 0.7");
    command("region box block 0 2 0 2 0 2");
    command("create_box 1 box");
    command("create_atoms 1 box");
    command("mass 1 1.0");
    command("pair_style zero 2.5");
    command("pair_coeff * *");
    command("variable beadtemp world 0.6 1.2");
    command("velocity all create ${beadtemp} 97531 mom yes rot no dist gaussian");
    command("timestep 0.002");
    // Equal mu and dE/dN leave only thermostat scaling of the electronic velocity.
    command("variable dEdN equal 1.5");
    command((std::string("fix cp all pimd/uvt method nmpimd integrator ") + split +
             " temp 0.8 Tdamp 0.2 tchain 3 tloop 3 drag 0.1 removecom no "
             "mu 1.5 Udamp 0.2 ne 1.7 ne_velocity 0.1 dedn v_dEdN").c_str());
    command("run 20 post no");

    double eta0 = pimd_test::fix_value(lmp, "cp", pimd_test::kNuclearPrefixScalars + 6);
    double eta_sum = 0.0, eta_min = 0.0, eta_max = 0.0;
    MPI_Allreduce(&eta0, &eta_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&eta0, &eta_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&eta0, &eta_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    EXPECT_GT(eta_max - eta_min, 1.0e-6);
    const double nedot = pimd_test::fix_value(lmp, "cp", pimd_test::uvt_vector_indices().ne_dot);
    // Integral of mean friction equals the mean displacement of eta[0].
    EXPECT_NEAR(nedot, 0.1 * std::exp(-eta_sum / nprocs), 1.0e-12);
    lammps_close(lmp);
  }
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
