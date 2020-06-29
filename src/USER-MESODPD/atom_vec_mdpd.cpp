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

#include "atom_vec_mdpd.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecMDPD::AtomVecMDPD(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;
  mass_type = 1;
  forceclearflag = 1;

  atom->rho_flag = 1;
  atom->vest_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = (char *) "rho drho vest";
  fields_copy = (char *) "rho drho vest";
  fields_comm = (char *) "rho vest";
  fields_comm_vel = (char *) "rho vest";
  fields_reverse = (char *) "drho";
  fields_border = (char *) "rho vest";
  fields_border_vel = (char *) "rho vest";
  fields_exchange = (char *) "rho vest";
  fields_restart = (char * ) "rho vest";
  fields_create = (char *) "rho drho vest";
  fields_data_atom = (char *) "id type rho x";
  fields_data_vel = (char *) "id v";

  setup_fields();
}

/* ---------------------------------------------------------------------- */

void AtomVecMDPD::init()
{
  AtomVec::init();

  if (strcmp(update->unit_style,"lj") != 0)
    error->all(FLERR,"Atom style mdpd requires lj units");
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecMDPD::grow_pointers()
{
  rho = atom->rho;
  drho = atom->drho;
  vest = atom->vest;
}

/* ----------------------------------------------------------------------
   clear extra forces starting at atom N
   nbytes = # of bytes to clear for a per-atom vector
------------------------------------------------------------------------- */

void AtomVecMDPD::force_clear(int n, size_t nbytes)
{
  memset(&drho[n],0,nbytes);
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecMDPD::data_atom_post(int ilocal)
{
  drho[ilocal] = 0.0;
  vest[ilocal][0] = 0.0;
  vest[ilocal][1] = 0.0;
  vest[ilocal][2] = 0.0;
}

/* ----------------------------------------------------------------------
   assign an index to named atom property and return index
   return -1 if name is unknown to this atom style
------------------------------------------------------------------------- */

int AtomVecMDPD::property_atom(char *name)
{
  if (strcmp(name,"rho") == 0) return 0;
  if (strcmp(name,"drho") == 0) return 1;
  return -1;
}

/* ----------------------------------------------------------------------
   pack per-atom data into buf for ComputePropertyAtom
   index maps to data specific to this atom style
------------------------------------------------------------------------- */

void AtomVecMDPD::pack_property_atom(int index, double *buf,
                                     int nvalues, int groupbit)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int n = 0;
  if (index == 0) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = rho[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 1) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = drho[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  }
}
