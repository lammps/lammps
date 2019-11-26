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

#include "atom_vec_angle.h"
#include "atom.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecAngle::AtomVecAngle(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 1;
  bonds_allow = angles_allow = 1;
  mass_type = 1;

  atom->molecule_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in the string does not matter
  //   except fields_data_atom and fields_data_vel which must match data file

  fields_grow = (char *) 
    "molecule num_bond bond_type bond_atom "
    "num_angle angle_type angle_atom1 angle_atom2 angle_atom3 nspecial special";
  fields_copy = (char *)
    "molecule num_bond bond_type bond_atom "
    "num_angle angle_type angle_atom1 angle_atom2 angle_atom3 nspecial special";
  fields_comm = NULL;
  fields_comm_vel = NULL;
  fields_reverse = NULL;
  fields_border = (char *) "molecule";
  fields_border_vel = (char *) "molecule";
  fields_exchange = (char *)
    "molecule num_bond bond_type bond_atom "
    "num_angle angle_type angle_atom1 angle_atom2 angle_atom3 nspecial special";
  fields_restart = (char *) 
    "molecule num_bond bond_type bond_atom "
    "num_angle angle_type angle_atom1 angle_atom2 angle_atom3";
  fields_create = (char *) "molecule num_bond num_angle nspecial";
  fields_data_atom = (char *) "id molecule type x";
  fields_data_vel = NULL;

  setup_fields();

  bond_per_atom = angle_per_atom = 0;
  bond_negative = angle_negative = NULL;
}

/* ---------------------------------------------------------------------- */

AtomVecAngle::~AtomVecAngle()
{
  delete [] bond_negative;
  delete [] angle_negative;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file
   modify/unmodify values for default AtomVec::pack_restart() to pack
------------------------------------------------------------------------- */

int AtomVecAngle::pack_restart(int i, double *buf)
{
  // insure negative vectors are needed length

  if (bond_per_atom < atom->bond_per_atom) {
    delete [] bond_negative;
    bond_per_atom = atom->bond_per_atom;
    bond_negative = new int[bond_per_atom];
  }
  if (angle_per_atom < atom->angle_per_atom) {
    delete [] angle_negative;
    angle_per_atom = atom->angle_per_atom;
    angle_negative = new int[angle_per_atom];
  }

  // flip any negative types to positive and flag which ones

  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  int *num_angle = atom->num_angle;
  int **angle_type = atom->angle_type;

  int any_bond_negative = 0;
  for (int m = 0; m < num_bond[i]; m++) {
    if (bond_type[i][m] < 0) {
      bond_negative[m] = 1;
      bond_type[i][m] = -bond_type[i][m];
      any_bond_negative = 1;
    } else bond_negative[m] = 0;
  }

  int any_angle_negative = 0;
  for (int m = 0; m < num_angle[i]; m++) {
    if (angle_type[i][m] < 0) {
      angle_negative[m] = 1;
      angle_type[i][m] = -angle_type[i][m];
      any_angle_negative = 1;
    } else angle_negative[m] = 0;
  }

  // perform the pack with adjusted values

  int n = AtomVec::pack_restart(i,buf);

  // restore the flagged types to their negative values

  if (any_bond_negative) {
    for (int m = 0; m < num_bond[i]; m++)
      if (bond_negative[m]) bond_type[i][m] = -bond_type[i][m];
  }
  if (any_angle_negative) {
    for (int m = 0; m < num_angle[i]; m++)
      if (angle_negative[m]) angle_type[i][m] = -angle_type[i][m];
  }

  return n;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
   initialize other atom quantities
------------------------------------------------------------------------- */

int AtomVecAngle::unpack_restart(double *buf)
{
  AtomVec::unpack_restart(buf);
  int ilocal = atom->nlocal-1;

  atom->nspecial[ilocal][0] = 0;
  atom->nspecial[ilocal][1] = 0;
  atom->nspecial[ilocal][2] = 0;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   modify what default AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecAngle::data_atom(double *coord, imageint imagetmp, char **values)
{
  AtomVec::data_atom(coord,imagetmp,values);
  int ilocal = atom->nlocal-1;

  atom->num_bond[ilocal] = 0;
  atom->num_angle[ilocal] = 0;
  atom->nspecial[ilocal][0] = 0;
  atom->nspecial[ilocal][1] = 0;
  atom->nspecial[ilocal][2] = 0;
}
