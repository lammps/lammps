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

#include "atom_vec_tdpd.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecTDPD::AtomVecTDPD(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;
  mass_type = 1;
  forceclearflag = 1;

  atom->tdpd_flag = 1;
  atom->vest_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = (char *) "cc cc_flux vest";
  fields_copy = (char *) "cc vest";
  fields_comm = (char *) "cc vest";
  fields_comm_vel = (char *) "cc vest";
  fields_reverse = (char *) "cc_flux";
  fields_border = (char *) "cc vest";
  fields_border_vel = (char *) "cc vest";
  fields_exchange = (char *) "cc vest";
  fields_restart = (char * ) "cc vest";
  fields_create = (char *) "cc vest";
  fields_data_atom = (char *) "id type x cc";
  fields_data_vel = (char *) "id v";
}

/* ----------------------------------------------------------------------
   process additional args
   single arg = number of cc_species
------------------------------------------------------------------------- */

void AtomVecTDPD::process_args(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Invalid atom_style tdpd command");

  atom->cc_species = utils::inumeric(FLERR,arg[0],false,lmp);
  cc_species = atom->cc_species;

  atom->add_peratom_change_columns("cc",cc_species);
  atom->add_peratom_change_columns("cc_flux",cc_species);

  // delay setting up of fields until now

  setup_fields();
}

/* ---------------------------------------------------------------------- */

void AtomVecTDPD::init()
{
  AtomVec::init();

  if (strcmp(update->unit_style,"lj") != 0)
    error->all(FLERR,"Atom style tdpd requires lj units");
}

/* ----------------------------------------------------------------------
   clear extra forces starting at atom N
   nbytes = # of bytes to clear for a per-atom vector
------------------------------------------------------------------------- */

void AtomVecTDPD::force_clear(int n, size_t nbytes)
{
  memset(&atom->cc_flux[n][0],0,cc_species*nbytes);
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecTDPD::data_atom_post(int ilocal)
{
  atom->vest[ilocal][0] = 0.0;
  atom->vest[ilocal][1] = 0.0;
  atom->vest[ilocal][2] = 0.0;
}
