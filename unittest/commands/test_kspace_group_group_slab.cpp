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

// Regression tests for "compute group/group ... kspace yes" with the
// Yeh-Berkowitz EW3DC slab correction on a triclinic (xy-tilted) box.  The
// group/group slab path (KSpace::slabcorr_groups) was only enabled for
// triclinic cells once the EW3DC correction was generalized to non-orthogonal
// boxes, and it is not exercised by the force-style YAML tests (which never
// invoke compute group/group).  Two checks:
//   (1) the triclinic xy-tilted result must match the orthogonal result, and
//   (2) the orthogonal result must match the analytic two-charged-sheet value.
//
// Physical setup (the reproducer from lammps/lammps#411): two atoms of charge
// +/-1 in a box that is periodic in x,y and non-periodic in z forms two
// infinite, oppositely charged sheets.  The group/group force between the
// sheets along z converges to 2*pi*q^2/(Lx*Ly) = pi and the interaction energy
// to 2*pi*q^2*d/(Lx*Ly) = 5*pi for Lx=1, Ly=2, q=1, d=|z2-z1|=5.

#include "../testing/core.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "gtest/gtest.h"
#include "modify.h"
#include "utils.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <mpi.h>
#include <string>
#include <vector>

bool verbose = false;

namespace LAMMPS_NS {

class KSpaceGroupGroupSlabTest : public LAMMPSTest {
 protected:
  // {energy, fx, fy, fz} of the group A <-> group B interaction
  using GroupGroup = std::array<double, 4>;

  void build_two_sheet_system()
  {
    command("clear");
    command("units lj");
    command("atom_style charge");
    command("boundary p p f");
    command("region box block 0.0 1.0 0.0 2.0 0.0 10.0 units box");
    command("create_box 1 box");
    command("create_atoms 1 single 0.0 0.0 0.0 units box");
    command("create_atoms 1 single 0.5 1.0 5.0 units box");
    command("set atom 1 charge  1.0");
    command("set atom 2 charge -1.0");
    command("mass 1 1.0");

    ASSERT_EQ(lmp->atom->natoms, 2);
  }

  // run a single group/group calculation; tilt the box into a triclinic cell
  // with an xy tilt when tilt_xy != 0.0 (xz = yz = 0 is required for EW3DC).
  GroupGroup run_group_group(const std::string &kspace, double tilt_xy)
  {
    build_two_sheet_system();

    HIDE_OUTPUT([&] {
      command("pair_style lj/cut/coul/long 4.5");
      command("pair_coeff * * 0.0 0.0");
      command("pair_modify table 0");
      command("kspace_style " + kspace + " 1.0e-10");
      command("kspace_modify slab 3.0");

      if (tilt_xy != 0.0) {
        command(fmt::format("change_box all triclinic xy final {:.15g} "
                            "xz final 0.0 yz final 0.0",
                            tilt_xy));
        // the kspace style must be re-created after switching to triclinic
        command("kspace_style " + kspace + " 1.0e-10");
        command("kspace_modify slab 3.0");
      }

      command("group pos id 1");
      command("group neg id 2");
      command("compute gg pos group/group neg kspace yes pair yes");
      command("run 0 post no");
    });

    // compute_scalar() runs the full pair+kspace group/group sum and fills
    // both scalar and vector[], so a single call yields energy and force
    auto *gg = lmp->modify->get_compute_by_id("gg");
    const double energy = gg->compute_scalar();
    return {energy, gg->vector[0], gg->vector[1], gg->vector[2]};
  }
};

TEST_F(KSpaceGroupGroupSlabTest, EwaldTriclinicSlabMatchesOrthogonal)
{
  if (!info->has_style("kspace", "ewald")) GTEST_SKIP();
  if (!info->has_style("pair", "coul/long")) GTEST_SKIP();
  if (!info->has_style("compute", "group/group")) GTEST_SKIP();

  const auto ortho = run_group_group("ewald", 0.0);
  const auto tri   = run_group_group("ewald", 0.3);

  if (verbose && (lmp->comm->me == 0))
    utils::print("ewald  ortho E={} fz={}   tri E={} fz={}\n", ortho[0], ortho[3], tri[0],
                 tri[3]);

  // analytic anchor for the orthogonal reference (two charged sheets)
  const double MY_PI = 3.14159265358979323846;
  EXPECT_NEAR(ortho[0], 5.0 * MY_PI, 1.0e-5);    // energy = 2*pi*q^2*d/(Lx*Ly)
  EXPECT_NEAR(ortho[3], 1.0 * MY_PI, 1.0e-5);    // z force = 2*pi*q^2/(Lx*Ly)
  EXPECT_NEAR(ortho[1], 0.0, 1.0e-10);           // x force
  EXPECT_NEAR(ortho[2], 0.0, 1.0e-10);           // y force

  // the triclinic xy-tilted result must reproduce the orthogonal one (Ewald is
  // exact, so this is bit-close, not just within slab accuracy)
  EXPECT_NEAR(tri[0], ortho[0], 1.0e-10);
  EXPECT_NEAR(tri[1], ortho[1], 1.0e-10);
  EXPECT_NEAR(tri[2], ortho[2], 1.0e-10);
  EXPECT_NEAR(tri[3], ortho[3], 1.0e-10);
}

TEST_F(KSpaceGroupGroupSlabTest, PPPMTriclinicSlabMatchesOrthogonal)
{
  if (!info->has_style("kspace", "pppm")) GTEST_SKIP();
  if (!info->has_style("pair", "coul/long")) GTEST_SKIP();
  if (!info->has_style("compute", "group/group")) GTEST_SKIP();

  const auto ortho = run_group_group("pppm", 0.0);
  const auto tri   = run_group_group("pppm", 0.3);

  if (verbose && (lmp->comm->me == 0))
    utils::print("pppm   ortho E={} fz={}   tri E={} fz={}\n", ortho[0], ortho[3], tri[0],
                 tri[3]);

  const double MY_PI = 3.14159265358979323846;
  // PPPM is mesh-based; allow a looser analytic tolerance than Ewald
  EXPECT_NEAR(ortho[0], 5.0 * MY_PI, 1.0e-3);    // energy = 2*pi*q^2*d/(Lx*Ly)
  EXPECT_NEAR(ortho[3], 1.0 * MY_PI, 1.0e-3);    // z force = 2*pi*q^2/(Lx*Ly)
  EXPECT_NEAR(ortho[1], 0.0, 1.0e-8);            // x force
  EXPECT_NEAR(ortho[2], 0.0, 1.0e-8);            // y force

  // orthogonal vs triclinic PPPM pick slightly different meshes/gewald, so
  // compare within a relative tolerance rather than bit-for-bit
  EXPECT_NEAR(tri[0], ortho[0], 1.0e-4 * std::abs(ortho[0]));
  EXPECT_NEAR(tri[3], ortho[3], 1.0e-4 * std::abs(ortho[3]));
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
      if (arg == "-v") verbose = true;
    }
  }

  if ((argc > 1) && (std::string(argv[1]) == "-v")) verbose = true;

  const int rv = RUN_ALL_TESTS();
  MPI_Finalize();
  return rv;
}
