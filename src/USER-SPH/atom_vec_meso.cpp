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

#include "atom_vec_meso.h"
#include <cstring>
#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecMeso::AtomVecMeso(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;
  mass_type = 1;
  forceclearflag = 1;

  atom->e_flag = 1;
  atom->rho_flag = 1;
  atom->cv_flag = 1;
  atom->vest_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = (char *) "rho drho e de cv vest";
  fields_copy = (char *) "rho drho e de cv vest";
  fields_comm = (char *) "rho e vest";
  fields_comm_vel = (char *) "rho e vest";
  fields_reverse = (char *) "drho de";
  fields_border = (char *) "rho e cv vest";
  fields_border_vel = (char *) "rho e cv vest";
  fields_exchange = (char *) "rho e cv vest";
  fields_restart = (char * ) "rho e cv vest";
  fields_create = (char *) "rho e cv vest de drho";
  fields_data_atom = (char *) "id type rho e cv x";
  fields_data_vel = (char *) "id v";

  setup_fields();
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecMeso::grow_pointers()
{
  rho = atom->rho;
  drho = atom->drho;
  e = atom->e;
  de = atom->de;
  cv = atom->cv;
  vest = atom->vest;
}

/* ----------------------------------------------------------------------
   clear extra forces starting at atom N
   nbytes = # of bytes to clear for a per-atom vector
------------------------------------------------------------------------- */

void AtomVecMeso::force_clear(int n, size_t nbytes)
{
  memset(&de[n],0,nbytes);
  memset(&drho[n],0,nbytes);
}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecMeso::create_atom_post(int ilocal)
{
  cv[ilocal] = 1.0;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecMeso::data_atom_post(int ilocal)
{
  vest[ilocal][0] = 0.0;
  vest[ilocal][1] = 0.0;
  vest[ilocal][2] = 0.0;
  de[ilocal] = 0.0;
  drho[ilocal] = 0.0;
}

/* ----------------------------------------------------------------------
   assign an index to named atom property and return index
   return -1 if name is unknown to this atom style
------------------------------------------------------------------------- */

int AtomVecMeso::property_atom(char *name)
{
  if (strcmp(name,"rho") == 0) return 0;
  if (strcmp(name,"drho") == 0) return 1;
  if (strcmp(name,"e") == 0) return 2;
  if (strcmp(name,"de") == 0) return 3;
  if (strcmp(name,"cv") == 0) return 4;
  return -1;
}

/* ----------------------------------------------------------------------
   pack per-atom data into buf for ComputePropertyAtom
   index maps to data specific to this atom style
------------------------------------------------------------------------- */

void AtomVecMeso::pack_property_atom(int index, double *buf,
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
  } else if (index == 2) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = e[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 3) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = de[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 4) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = cv[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  }
}
