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

class FixTPRPMDCompatibilityTest : public LAMMPSTest {
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

  void setup_uvt_fix(const char *style, double ne = 1.4, double ne_velocity = 0.0)
  {
    command("variable k_quad equal 5.0");
    command("variable N0_quad equal 1.0");
    command(pimd_test::quadratic_dedn_variable().c_str());
    command(std::string("fix cp all ") + style +
            " method tp-rpmd ensemble uvt thermostat NHC temp 1.0 Tdamp 0.5 "
            "tchain 3 tloop 1 mu 2.0 2.0 0.5 ne " + std::to_string(ne) +
            " ne_velocity " + std::to_string(ne_velocity) + " dedn v_dEdN");
  }

  std::vector<double> sample_fix_vector(int size)
  {
    std::vector<double> values(size);
    for (int i = 0; i < size; ++i) values[i] = pimd_test::fix_value(lmp, "cp", i);
    return values;
  }
};

TEST_F(FixTPRPMDCompatibilityTest, AliasStyleParsesAndRuns)
{
  const auto uvt = pimd_test::uvt_vector_indices();

  setup_quadratic_system();
  setup_uvt_fix("tprpmd");
  command("run 0 post no");

  EXPECT_NEAR(pimd_test::fix_value(lmp, "cp", uvt.ne), 1.4, 1.0e-10);
  EXPECT_NEAR(pimd_test::fix_value(lmp, "cp", uvt.dedn), 2.0, 1.0e-10);
}

TEST_F(FixTPRPMDCompatibilityTest, AliasMatchesPIMDUVTKeyOutputs)
{
  const int vector_size = pimd_test::nuclear_vector_size() + 6;

  setup_quadratic_system();
  setup_uvt_fix("tprpmd");
  command("run 0 post no");
  auto alias_values = sample_fix_vector(vector_size);

  command("clear");
  setup_quadratic_system();
  setup_uvt_fix("pimd/uvt");
  command("run 0 post no");
  auto uvt_values = sample_fix_vector(vector_size);

  for (int i = 0; i < vector_size; ++i) {
    EXPECT_NEAR(alias_values[i], uvt_values[i], 1.0e-12) << "index=" << i;
  }
}

TEST_F(FixTPRPMDCompatibilityTest, AliasRestartRestoresIntoPIMDUVT)
{
  const auto uvt = pimd_test::uvt_vector_indices();

  setup_quadratic_system();
  setup_uvt_fix("tprpmd", 1.8, 0.0);
  command("run 20 post no");
  command("write_restart tprpmd_compat.restart");

  command("clear");
  command("read_restart tprpmd_compat.restart");
  setup_uvt_fix("pimd/uvt", 1.8, 0.0);
  command("run 0 post no");
  platform::unlink("tprpmd_compat.restart");

  EXPECT_TRUE(std::isfinite(pimd_test::fix_value(lmp, "cp", uvt.ne)));
  EXPECT_TRUE(std::isfinite(pimd_test::fix_value(lmp, "cp", uvt.ne_dot)));
  EXPECT_TRUE(std::isfinite(pimd_test::fix_value(lmp, "cp", uvt.dedn)));
  EXPECT_TRUE(std::isfinite(pimd_test::fix_value(lmp, "cp", uvt.mu)));
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
