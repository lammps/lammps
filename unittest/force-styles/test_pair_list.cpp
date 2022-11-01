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

#include "library.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

const char parms[] = "print \"\"\"\n"
                     "1 2 lj126 0.3 3.5 4.0\n"
                     "2 3 lj126 0.3 3.5 5.0\n"
                     "1 4 harmonic 10.0 7.0\n"
                     "2 4 harmonic 10.0 7.0\n"
                     "3 4 morse 10.0 0.3 6.5\n"
                     "\"\"\" file list.param\n";

const char first[] = "units           real\n"
                     "atom_style      bond\n"
                     "atom_modify map array\n"
                     "boundary        f f f\n"
                     "special_bonds   lj/coul 0.0 1.0 1.0\n"
                     "region          box block -5 5 -5 5 -5 5\n"
                     "create_box      1 box bond/types 2 extra/bond/per/atom 4\n"
                     "create_atoms 1 single -2.0  0.0 0.0\n"
                     "create_atoms 1 single  0.0  1.0 0.0\n"
                     "create_atoms 1 single  4.0  1.0 0.0\n"
                     "create_atoms 1 single  4.0 -4.0 -4.0\n"
                     "create_bonds single/bond 1 1 4\n"
                     "create_bonds single/bond 1 2 4\n"
                     "create_bonds single/bond 2 3 4\n"
                     "mass            1 10.0\n"
                     "velocity        all create 10.0 87287 loop geom\n";

const char second[] = "timestep        0.2\n"
                      "fix             1 all nve\n"
                      "run 2 post no\n";

static constexpr double EPSILON = 1.0e-10;

namespace LAMMPS_NS {

TEST(PairList, ListVsPairBond)
{
    if (!lammps_config_has_package("MOLECULE")) GTEST_SKIP();
    if (!lammps_config_has_package("MISC")) GTEST_SKIP();

    const char *lmpargv[] = {"melt", "-log", "none", "-nocite"};
    int lmpargc           = sizeof(lmpargv) / sizeof(const char *);

    ::testing::internal::CaptureStdout();
    void *ljmelt = lammps_open_no_mpi(lmpargc, (char **)lmpargv, nullptr);
    lmpargv[0]   = "plist";
    void *plist  = lammps_open_no_mpi(lmpargc, (char **)lmpargv, nullptr);

    lammps_commands_string(ljmelt, first);
    lammps_command(ljmelt, "pair_style lj/cut 5.0");
    lammps_command(ljmelt, "pair_coeff * * 0.3 3.5");
    lammps_command(ljmelt, "bond_style hybrid harmonic morse");
    lammps_command(ljmelt, "bond_coeff 1 harmonic 10.0 7.0");
    lammps_command(ljmelt, "bond_coeff 2 morse 10.0 0.3 6.5");

    lammps_commands_string(ljmelt, second);

    lammps_command(plist, parms);
    lammps_commands_string(plist, first);
    lammps_command(plist, "pair_style list list.param 10.0");
    lammps_command(plist, "pair_coeff * *");
    lammps_command(plist, "bond_style zero");
    lammps_command(plist, "bond_coeff * 2.0");
    lammps_commands_string(plist, second);
    ::testing::internal::GetCapturedStdout();

    double lj_pe = lammps_get_thermo(ljmelt, "pe");
    double ml_pe = lammps_get_thermo(plist, "pe");
    EXPECT_NEAR(lj_pe, ml_pe, EPSILON);
    double lj_ke = lammps_get_thermo(ljmelt, "ke");
    double ml_ke = lammps_get_thermo(plist, "ke");
    EXPECT_NEAR(lj_ke, ml_ke, EPSILON);
    double lj_press = lammps_get_thermo(ljmelt, "press");
    double ml_press = lammps_get_thermo(plist, "press");
    EXPECT_NEAR(lj_press, ml_press, EPSILON);

    ::testing::internal::CaptureStdout();
    lammps_command(plist, "shell rm list.param");
    lammps_close(ljmelt);
    lammps_close(plist);
    ::testing::internal::GetCapturedStdout();
}

} // namespace LAMMPS_NS
