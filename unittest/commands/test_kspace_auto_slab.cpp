/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "../testing/core.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "gtest/gtest.h"
#include "utils.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <mpi.h>
#include <string>
#include <utility>
#include <vector>

bool verbose = false;

namespace LAMMPS_NS {

class KSpaceAutoSlabTest : public LAMMPSTest {
 protected:
  using ForceVector = std::vector<std::array<double, 3>>;

  void build_q2d_system()
  {
    command("clear");
    command("units real");
    command("atom_style charge");
    command("boundary p p f");
    command("region box block 0 40 0 40 0 40");
    command("create_box 1 box");

    uint64_t state = 246813579ULL;
    auto next_uniform = [&state]() {
      state = state * 6364136223846793005ULL + 1442695040888963407ULL;
      return static_cast<double>((state >> 11) & ((1ULL << 53) - 1)) /
          static_cast<double>(1ULL << 53);
    };

    for (int atom_id = 1; atom_id <= 100; ++atom_id) {
      const double x = 1.0 + 38.0 * next_uniform();
      const double y = 1.0 + 38.0 * next_uniform();
      const double z = 1.0 + 38.0 * next_uniform();
      command(fmt::format("create_atoms 1 single {:.15g} {:.15g} {:.15g}", x, y, z));
      const double charge = (atom_id % 2 == 1) ? 1.0 : -1.0;
      command(fmt::format("set atom {} charge {:.1f}", atom_id, charge));
    }

    command("mass 1 1.0");

    ASSERT_DOUBLE_EQ(lmp->domain->xprd, 40.0);
    ASSERT_DOUBLE_EQ(lmp->domain->yprd, 40.0);
    ASSERT_DOUBLE_EQ(lmp->domain->zprd, 40.0);
    ASSERT_EQ(lmp->atom->natoms, 100);
  }

  ForceVector run_total_ewald_forces(const std::string &accuracy, const std::string &slab_setting)
  {
    build_q2d_system();

    HIDE_OUTPUT([&] {
      command("pair_style coul/long 8.0");
      command("pair_coeff * *");
      command("pair_modify table 0");
      command("kspace_style ewald " + accuracy);
      command("kspace_modify slab " + slab_setting);
      command("run 0 post no");
    });

    return snapshot_forces();
  }

  ForceVector snapshot_forces() const
  {
    ForceVector forces(lmp->atom->natoms);
    auto **f = lmp->atom->f;
    tagint *tag = lmp->atom->tag;
    for (int i = 0; i < lmp->atom->nlocal; ++i) {
      const int idx = tag[i] - 1;
      forces[idx] = {f[i][0], f[i][1], f[i][2]};
    }
    return forces;
  }

  double relative_l2_error(const ForceVector &reference, const ForceVector &trial) const
  {
    double numerator = 0.0;
    double denominator = 0.0;
    for (size_t i = 0; i < reference.size(); ++i) {
      for (int dim = 0; dim < 3; ++dim) {
        const double delta = trial[i][dim] - reference[i][dim];
        numerator += delta * delta;
        denominator += reference[i][dim] * reference[i][dim];
      }
    }

    if (denominator == 0.0) return 0.0;
    return std::sqrt(numerator / denominator);
  }
};

TEST_F(KSpaceAutoSlabTest, EwaldAutoMatchesTightReferenceWithinTenXAccuracy)
{
  if (!info->has_style("kspace", "ewald")) GTEST_SKIP();
  if (!info->has_style("pair", "coul/long")) GTEST_SKIP();

  const auto reference = run_total_ewald_forces("1.0e-8", "10.0");
  const std::vector<std::pair<std::string, double>> cases = {
      {"1.0e-3", 1.0e-3},
      {"1.0e-4", 1.0e-4},
      {"1.0e-5", 1.0e-5},
  };

  for (const auto &[accuracy, tolerance] : cases) {
    const auto trial = run_total_ewald_forces(accuracy, "auto");
    const double error = relative_l2_error(reference, trial);
    EXPECT_LT(error, 10.0 * tolerance) << "accuracy=" << accuracy;
    if (verbose && (lmp->comm->me == 0)) {
      utils::print("accuracy={} auto_force_l2_rel_err={} bound={}\n", accuracy, error,
                 10.0 * tolerance);
    }
  }
}

}    // namespace LAMMPS_NS

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleMock(&argc, argv);

  // handle arguments passed via environment variable
  if (const char *var = getenv("TEST_ARGS")) {
    std::vector<std::string> env = LAMMPS_NS::utils::split_words(var);
    for (auto arg : env) {
      if (arg == "-v") {
        verbose = true;
      }
    }
  }

  if ((argc > 1) && (std::string(argv[1]) == "-v")) verbose = true;

  const int rv = RUN_ALL_TESTS();
  MPI_Finalize();
  return rv;
}
