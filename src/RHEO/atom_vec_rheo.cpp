/* ----------------------------------------------------------------------
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
   Contributing authors:
   Joel Clemmer (SNL)
----------------------------------------------------------------------- */

#include "atom_vec_rheo.h"

#include "atom.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecRHEO::AtomVecRHEO(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::ATOMIC;
  mass_type = PER_TYPE;
  forceclearflag = 1;

  atom->rheo_status_flag = 1;
  atom->pressure_flag = 1;
  atom->rho_flag = 1;
  atom->viscosity_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = {"rheo_status", "rho", "drho", "pressure", "viscosity"};
  fields_copy = {"rheo_status", "rho", "drho", "pressure", "viscosity"};
  fields_comm = {"rheo_status", "rho"};
  fields_comm_vel = {"rheo_status", "rho"};
  fields_reverse = {"drho"};
  fields_border = {"rheo_status", "rho"};
  fields_border_vel = {"rheo_status", "rho"};
  fields_exchange = {"rheo_status", "rho"};
  fields_restart = {"rheo_status", "rho"};
  fields_create = {"rheo_status", "rho", "drho", "pressure", "viscosity"};
  fields_data_atom = {"id", "type", "rheo_status", "rho", "x"};
  fields_data_vel = {"id", "v"};

  setup_fields();
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecRHEO::grow_pointers()
{
  rheo_status = atom->rheo_status;
  pressure = atom->pressure;
  rho = atom->rho;
  drho = atom->drho;
  viscosity = atom->viscosity;
}

/* ----------------------------------------------------------------------
   clear extra forces starting at atom N
   nbytes = # of bytes to clear for a per-atom vector
------------------------------------------------------------------------- */

void AtomVecRHEO::force_clear(int n, size_t nbytes)
{
  memset(&drho[n], 0, nbytes);
}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecRHEO::create_atom_post(int ilocal)
{
  rho[ilocal] = 1.0;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecRHEO::data_atom_post(int ilocal)
{
  drho[ilocal] = 0.0;
  pressure[ilocal] = 0.0;
  viscosity[ilocal] = 0.0;
}

/* ----------------------------------------------------------------------
   assign an index to named atom property and return index
   return -1 if name is unknown to this atom style
------------------------------------------------------------------------- */

int AtomVecRHEO::property_atom(const std::string &name)
{
  if (name == "rheo_status") return 0;
  if (name == "pressure") return 1;
  if (name == "rho") return 2;
  if (name == "drho") return 3;
  if (name == "viscosity") return 4;
  return -1;
}

/* ----------------------------------------------------------------------
   pack per-atom data into buf for ComputePropertyAtom
   index maps to data specific to this atom style
------------------------------------------------------------------------- */

void AtomVecRHEO::pack_property_atom(int index, double *buf, int nvalues, int groupbit)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int n = 0;

  if (index == 0) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = rheo_status[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 1) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = pressure[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 2) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = rho[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 3) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = drho[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 4) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = viscosity[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  }
}
