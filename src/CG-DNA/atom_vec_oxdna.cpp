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

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"

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
  fields_border = (char *) "id5p";
  fields_border_vel = (char *) "";
  fields_exchange = (char *) "id5p";
  fields_restart = (char *) "id5p";
  fields_create = (char *) "";
  fields_data_atom = (char *) "id type x";
  fields_data_vel = (char *) "id v";

  setup_fields();

  if(!force->newton_bond) error->warning(FLERR,"Write_data command requires newton on to preserve 3'->5' bond polarity");
}

/* ---------------------------------------------------------------------- */
AtomVecOxdna::~AtomVecOxdna() {}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecOxdna::grow_pointers()
{
  id5p = atom->id5p;
}

/* ----------------------------------------------------------------------
   initialize atom quantity 5' partner
------------------------------------------------------------------------- */

void AtomVecOxdna::data_atom_post(int ilocal)
{
  tagint *id5p = atom->id5p;
  id5p[ilocal] = -1;
}

/* ----------------------------------------------------------------------
   process bond information as per data file
   store 5' partner to inform 3'->5' bond directionality
------------------------------------------------------------------------- */

void AtomVecOxdna::data_bonds_post(int /*m*/, int /*num_bond*/, tagint atom1, tagint atom2,
                                   tagint id_offset)
{
  int n;
  tagint *id5p = atom->id5p;

  if (id_offset) {
    atom1 += id_offset;
    atom2 += id_offset;
  }

  if ((n = atom->map(atom1)) >= 0) { id5p[n] = atom2; }
}
