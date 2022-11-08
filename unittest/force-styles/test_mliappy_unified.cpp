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

const char pickle[] = "python create_pickle here \"\"\"\n"
                      "import lammps\n"
                      "import lammps.mliap\n"
                      "from lammps.mliap.mliap_unified_lj import MLIAPUnifiedLJ\n"
                      "def create_pickle():\n"
                      "    unified = MLIAPUnifiedLJ(['Ar'])\n"
                      "    unified.pickle('mliap_unified_lj_Ar.pkl')\n"
                      "\"\"\"\n";

const char first[] = "units           lj\n"
                     "atom_style      atomic\n"
                     "lattice         fcc 0.8442\n"
                     "region          box block 0 2 0 2 0 2\n"
                     "create_box      1 box\n"
                     "create_atoms    1 box\n"
                     "mass            1 1.0\n"
                     "velocity        all create 3.0 87287 loop geom\n";

const char second[] = "neighbor        0.3 bin\n"
                      "neigh_modify    every 20 delay 0 check no\n"
                      "fix             1 all nve\n"
                      "run 2 post no\n";

namespace LAMMPS_NS {

TEST(MliapUnified, VersusLJMelt)
{
    const char *lmpargv[] = {"melt", "-log", "none", "-nocite"};
    int lmpargc           = sizeof(lmpargv) / sizeof(const char *);

    void *ljmelt = lammps_open_no_mpi(lmpargc, (char **)lmpargv, nullptr);
    void *mliap  = lammps_open_no_mpi(lmpargc, (char **)lmpargv, nullptr);

    lammps_commands_string(ljmelt, first);
    lammps_command(ljmelt, "pair_style lj/cut 2.5");
    lammps_command(ljmelt, "pair_coeff * * 1.0 1.0");
    lammps_commands_string(ljmelt, second);

    lammps_command(mliap, pickle);
    lammps_command(mliap, "python create_pickle invoke");

    lammps_commands_string(mliap, first);
    lammps_command(mliap, "pair_style mliap unified mliap_unified_lj_Ar.pkl 0");
    lammps_command(mliap, "pair_coeff * * Ar");
    lammps_commands_string(mliap, second);

    double lj_pe = lammps_get_thermo(ljmelt, "pe");
    double ml_pe = lammps_get_thermo(mliap, "pe");
    EXPECT_NEAR(lj_pe, ml_pe, 1.0e-14);
    double lj_ke = lammps_get_thermo(ljmelt, "ke");
    double ml_ke = lammps_get_thermo(mliap, "ke");
    EXPECT_NEAR(lj_ke, ml_ke, 1.0e-14);
    double lj_press = lammps_get_thermo(ljmelt, "press");
    double ml_press = lammps_get_thermo(mliap, "press");
    EXPECT_NEAR(lj_press, ml_press, 1.0e-14);

    lammps_command(mliap, "shell rm mliap_unified_lj_Ar.pkl");
    lammps_close(ljmelt);
    lammps_close(mliap);
}

TEST(MliapUnified, VersusLJMeltGhost)
{
    const char *lmpargv[] = {"melt", "-log", "none", "-nocite"};
    int lmpargc           = sizeof(lmpargv) / sizeof(const char *);

    void *ljmelt = lammps_open_no_mpi(lmpargc, (char **)lmpargv, nullptr);
    void *mliap  = lammps_open_no_mpi(lmpargc, (char **)lmpargv, nullptr);

    lammps_commands_string(ljmelt, first);
    lammps_command(ljmelt, "pair_style lj/cut 2.5");
    lammps_command(ljmelt, "pair_coeff * * 1.0 1.0");
    lammps_commands_string(ljmelt, second);

    lammps_command(mliap, pickle);
    lammps_command(mliap, "python create_pickle invoke");

    lammps_commands_string(mliap, first);
    lammps_command(mliap, "pair_style mliap unified mliap_unified_lj_Ar.pkl 1");
    lammps_command(mliap, "pair_coeff * * Ar");
    lammps_commands_string(mliap, second);

    double lj_pe = lammps_get_thermo(ljmelt, "pe");
    double ml_pe = lammps_get_thermo(mliap, "pe");
    EXPECT_NEAR(lj_pe, ml_pe, 1.0e-14);
    double lj_ke = lammps_get_thermo(ljmelt, "ke");
    double ml_ke = lammps_get_thermo(mliap, "ke");
    EXPECT_NEAR(lj_ke, ml_ke, 1.0e-14);
    double lj_press = lammps_get_thermo(ljmelt, "press");
    double ml_press = lammps_get_thermo(mliap, "press");
    EXPECT_NEAR(lj_press, ml_press, 1.0e-14);

    lammps_command(mliap, "shell rm mliap_unified_lj_Ar.pkl");
    lammps_close(ljmelt);
    lammps_close(mliap);
}

} // namespace LAMMPS_NS
