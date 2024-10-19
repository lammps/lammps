/* -*- c++ -*--------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   Contributing author:  Mitch Murphy (alphataubio at gmail)
------------------------------------------------------------------------- */

#include "charmm2reaxff.h"

#include "atom.h"
#include "error.h"
#include "input.h"

#include <cmath>

using namespace LAMMPS_NS;

static bool is_mass_equal(double m1, double m2) { return std::abs(m1-m2)<.001; }

/* ---------------------------------------------------------------------- */

void Charmm2Reaxff::command(int narg, char **arg)
{

  lmp->atom->molecular = Atom::ATOMIC;
  lmp->atom->avec->molecular = Atom::ATOMIC;
  lmp->atom->q_flag = 1;
  lmp->atom->deallocate_topology();
  lmp->atom->nbondtypes = 0;
  lmp->atom->nangletypes = 0;
  lmp->atom->ndihedraltypes = 0;
  lmp->atom->nimpropertypes = 0;
  lmp->atom->nbonds = 0;
  lmp->atom->nangles = 0;
  lmp->atom->ndihedrals = 0;
  lmp->atom->nimpropers = 0;

  if (narg < 1) utils::missing_cmd_args(FLERR, "charmm2reaxff", error);

  if (atom->tag_enable == 0) error->all(FLERR, "Must have atom IDs for charmm2reaxff command");

  int *type = lmp->atom->type;
  double *mass = lmp->atom->mass;

  for ( int i = 0; i < lmp->atom->nlocal; i++) {

    double mass_i = mass[type[i]];

    if(is_mass_equal(mass_i,12.011))      type[i] = 1;  // C
    else if(is_mass_equal(mass_i,1.008))  type[i] = 2;  // H
    else if(is_mass_equal(mass_i,15.999)) type[i] = 3;  // O
    else if(is_mass_equal(mass_i,14.007)) type[i] = 4;  // N
    else if(is_mass_equal(mass_i,32.060)) type[i] = 5;  // S
    else if(is_mass_equal(mass_i,24.305)) type[i] = 6;  // Mg
    else if(is_mass_equal(mass_i,22.990)) type[i] = 8;  // Na
    else if(is_mass_equal(mass_i,35.450)) type[i] = 10;  // Cl
    else
      error->all(FLERR, "charmm2reaxff: atom {} with unknown mass", i);

  }

  atom->ntypes = 13;

  mass[1] = 12.011;  // C
  mass[2] = 1.008;   // H
  mass[3] = 15.999;  // O
  mass[4] = 14.007;  // N
  mass[5] = 32.060;  // S
  mass[6] = 24.305;  // Mg
  mass[7] = 30.974;  // P
  mass[8] = 22.990;  // Na
  mass[9] = 47.867;  // Ti
  mass[10] = 35.450; // Cl
  mass[11] = 18.998; // F
  mass[12] = 196.97; // Au
  mass[13] = 0.0000; // X

  std::string pair_coeff_cmd = fmt::format("pair_coeff * * {} C H O N S Mg P Na Ti Cl F Au X", arg[0]);

  //input->one("newton on");
  //input->one("labelmap atom 1 C 2 H 3 O 4 N 5 S 6 Mg 7 P 8 Na 9 Ti 10 Cl 11 F 12 Au 13 X");
  input->one("bond_style none");
  input->one("angle_style none");
  input->one("dihedral_style none");
  input->one("improper_style none");
  input->one("unfix cmap");
  input->one("kspace_style none");
  input->one("pair_style reaxff NULL");
  input->one(pair_coeff_cmd);

}

