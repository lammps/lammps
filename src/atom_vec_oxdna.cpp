// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom_vec_oxdna.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
AtomVecOxdna::AtomVecOxdna(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::MOLECULAR;
  bonds_allow = 1;
  mass_type = PER_TYPE;

  atom->molecule_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = (char *) "id5p";
  fields_copy = (char *) "id5p";
  fields_comm = (char *) "";
  fields_comm_vel = (char *) "";
  fields_reverse = (char *) "";
  fields_border = (char *) "";
  fields_border_vel = (char *) "";
  fields_exchange = (char *) "id5p";
  fields_restart = (char *) "id5p";
  fields_create = (char *) "";
  fields_data_atom = (char *) "id type x";
  fields_data_vel = (char *) "id v";

  setup_fields();

}

/* ---------------------------------------------------------------------- */
AtomVecOxdna::~AtomVecOxdna()
{

}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecOxdna::grow_pointers()
{
  id5p = atom->id5p;

}

/* ----------------------------------------------------------------------
   process bond information as per data file
   store 5' partner to inform 3'->5' bond directionality 
------------------------------------------------------------------------- */

void AtomVecOxdna::data_bonds_post()
{

printf("CALLED FROM ATOM_VEC_OXDNA.CPP\n");

}
