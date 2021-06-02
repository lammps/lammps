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

#include "atom_vec_edpd.h"

#include "atom.h"
#include "error.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecEDPD::AtomVecEDPD(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::ATOMIC;
  mass_type = PER_TYPE;
  forceclearflag = 1;

  atom->edpd_flag = 1;
  atom->vest_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = (char *) "edpd_cv edpd_temp edpd_flux vest vest_temp";
  fields_copy = (char *) "edpd_cv edpd_temp edpd_flux vest vest_temp";
  fields_comm = (char *) "edpd_temp vest vest_temp";
  fields_comm_vel = (char *) "edpd_temp vest vest_temp";
  fields_reverse = (char *) "edpd_flux";
  fields_border = (char *) "edpd_cv edpd_temp vest vest_temp";
  fields_border_vel = (char *) "edpd_cv edpd_temp vest vest_temp";
  fields_exchange = (char *) "edpd_cv edpd_temp vest vest_temp";
  fields_restart = (char * ) "edpd_cv edpd_temp vest vest_temp";
  fields_create = (char *) "edpd_cv edpd_temp edpd_flux vest vest_temp";
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
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecEDPD::grow_pointers()
{
  edpd_cv = atom->edpd_cv;
  edpd_temp = atom->edpd_temp;
  edpd_flux = atom->edpd_flux;
  vest = atom->vest;
  vest_temp = atom->vest_temp;
}

/* ----------------------------------------------------------------------
   clear extra forces starting at atom N
   nbytes = # of bytes to clear for a per-atom vector
------------------------------------------------------------------------- */

void AtomVecEDPD::force_clear(int n, size_t nbytes)
{
  memset(&edpd_flux[n],0,nbytes);
}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecEDPD::create_atom_post(int ilocal)
{
  edpd_temp[ilocal] = 1.0;
  edpd_cv[ilocal]= 1.0e5;
  vest_temp[ilocal] = edpd_temp[ilocal];
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecEDPD::data_atom_post(int ilocal)
{
  edpd_flux[ilocal] = 0.0;
  vest[ilocal][0] = 0.0;
  vest[ilocal][1] = 0.0;
  vest[ilocal][2] = 0.0;
  vest_temp[ilocal] = edpd_temp[ilocal];
}
