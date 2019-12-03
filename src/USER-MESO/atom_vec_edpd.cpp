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

#include "atom_vec_edpd.h"
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "update.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecEDPD::AtomVecEDPD(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;
  mass_type = 1;
  forceclearflag = 1;

  atom->edpd_flag = 1;
  atom->vest_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = (char *) "edpd_cv edpd_temp edpd_flux vest";
  fields_copy = (char *) "edpd_cv edpd_temp edpd_flux vest";
  fields_comm = (char *) "edpd_temp vest";
  fields_comm_vel = (char *) "edpd_temp vest";
  fields_reverse = (char *) "edpd_flux";
  fields_border = (char *) "edpd_cv edpd_temp vest";
  fields_border_vel = (char *) "edpd_cv edpd_temp vest";
  fields_exchange = (char *) "edpd_cv edpd_temp vest";
  fields_restart = (char * ) "edpd_cv edpd_temp vest";
  fields_create = (char *) "edpd_cv edpd_temp edpd_flux vest";
  fields_data_atom = (char *) "id type edpd_temp edpd_cv x";
  fields_data_vel = (char *) "id v";

  setup_fields();
}

/* ---------------------------------------------------------------------- */

void AtomVecEDPD::init()
{
  AtomVec::init();

  if (strcmp(update->unit_style,"lj") != 0)
    error->all(FLERR,"Atom style edpd requires lj units");
}

/* ----------------------------------------------------------------------
   clear extra forces starting at atom N
   nbytes = # of bytes to clear for a per-atom vector
------------------------------------------------------------------------- */

void AtomVecEDPD::force_clear(int n, size_t nbytes)
{
  memset(&atom->edpd_flux[n],0,nbytes);
}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecEDPD::create_atom_post(int ilocal)
{
  atom->edpd_temp[ilocal] = 1.0;
  atom->edpd_cv[ilocal]= 1.0e5;
  atom->vest[ilocal][3] = atom->edpd_temp[ilocal];
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecEDPD::data_atom_post(int ilocal)
{
  atom->edpd_flux[ilocal] = 0.0;
  atom->vest[ilocal][0] = 0.0;
  atom->vest[ilocal][1] = 0.0;
  atom->vest[ilocal][2] = 0.0;
  atom->vest[ilocal][3] = atom->edpd_temp[ilocal];
}
