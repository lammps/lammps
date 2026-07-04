/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with the U.S. Government retains
   certain rights in this software.
------------------------------------------------------------------------- */

#define LAMMPS_LIB_MPI 1

#include "lammps.h"
#include "utils.h"

#include "gtest/gtest.h"
#include "pimd_langevin_test_utils.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <mpi.h>
#include <vector>

bool verbose = false;

namespace LAMMPS_NS {

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
  ::testing::InitGoogleTest(&argc, argv);

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
