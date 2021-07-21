/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
#ifndef TEST_EXAMPLE_MELT__H
#define TEST_EXAMPLE_MELT__H

#include "../core.h"

class MeltTest : public LAMMPSTest {
protected:
    virtual void InitSystem() override
    {
        HIDE_OUTPUT([&] {
            command("units           lj");
            command("atom_style      atomic");
            command("atom_modify     map yes");

            command("lattice         fcc 0.8442");
            command("region          box block 0 2 0 2 0 2");
            command("create_box      1 box");
            command("create_atoms    1 box");
            command("mass            1 1.0");

            command("velocity        all create 3.0 87287");

            command("pair_style      lj/cut 2.5");
            command("pair_coeff      1 1 1.0 1.0 2.5");

            command("neighbor        0.3 bin");
            command("neigh_modify    every 20 delay 0 check no");
        });
    }
};

#endif
