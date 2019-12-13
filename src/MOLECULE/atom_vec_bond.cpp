/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom_vec_bond.h"
#include "atom.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecBond::AtomVecBond(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 1;
  bonds_allow = 1;
  mass_type = 1;

  atom->molecule_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = (char *) 
    "molecule num_bond bond_type bond_atom nspecial special";
  fields_copy = (char *)
    "molecule num_bond bond_type bond_atom nspecial special";
  fields_comm = (char *) "";
  fields_comm_vel = (char *) "";
  fields_reverse = (char *) "";
  fields_border = (char *) "molecule";
  fields_border_vel = (char *) "molecule";
  fields_exchange = (char *)
    "molecule num_bond bond_type bond_atom nspecial special";
  fields_restart = (char *) "molecule num_bond bond_type bond_atom";
  fields_create = (char *) "molecule num_bond nspecial";
  fields_data_atom = (char *) "id molecule type x";
  fields_data_vel = (char *) "id v";

  setup_fields();

  bond_per_atom = 0;
  bond_negative = NULL;
}

/* ---------------------------------------------------------------------- */

AtomVecBond::~AtomVecBond()
{
  delete [] bond_negative;
}

/* ----------------------------------------------------------------------
   grow atom arrays
   must set local copy of body ptr
   needed in replicate when 2 atom classes exist and pack_restart() is called
------------------------------------------------------------------------- */

void AtomVecBond::grow(int n)
{
  AtomVec::grow(n);
  num_bond = atom->num_bond;
  bond_type = atom->bond_type;
}


/* ----------------------------------------------------------------------
   modify values for AtomVec::pack_restart() to pack
------------------------------------------------------------------------- */

void AtomVecBond::pack_restart_pre(int ilocal)
{
  // insure bond_negative vector is needed length

  if (bond_per_atom < atom->bond_per_atom) {
    delete [] bond_negative;
    bond_per_atom = atom->bond_per_atom;
    bond_negative = new int[bond_per_atom];
  }

  // flip any negative types to positive and flag which ones

  //int *num_bond = atom->num_bond;
  //int **bond_type = atom->bond_type;

  any_bond_negative = 0;
  for (int m = 0; m < num_bond[ilocal]; m++) {
    if (bond_type[ilocal][m] < 0) {
      bond_negative[m] = 1;
      bond_type[ilocal][m] = -bond_type[ilocal][m];
      any_bond_negative = 1;
    } else bond_negative[m] = 0;
  }
}

/* ----------------------------------------------------------------------
   unmodify values packed by AtomVec::pack_restart()
------------------------------------------------------------------------- */

void AtomVecBond::pack_restart_post(int ilocal)
{
  // restore the flagged types to their negative values

  if (any_bond_negative) {
    //int *num_bond = atom->num_bond;
    //int **bond_type = atom->bond_type;
    for (int m = 0; m < num_bond[ilocal]; m++)
      if (bond_negative[m]) bond_type[ilocal][m] = -bond_type[ilocal][m];
  }
}

/* ----------------------------------------------------------------------
   initialize other atom quantities after AtomVec::unpack_restart()
------------------------------------------------------------------------- */

void AtomVecBond::unpack_restart_init(int ilocal)
{
  atom->nspecial[ilocal][0] = 0;
  atom->nspecial[ilocal][1] = 0;
  atom->nspecial[ilocal][2] = 0;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecBond::data_atom_post(int ilocal)
{
  atom->num_bond[ilocal] = 0;
  atom->nspecial[ilocal][0] = 0;
  atom->nspecial[ilocal][1] = 0;
  atom->nspecial[ilocal][2] = 0;
}
