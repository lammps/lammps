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

/* ----------------------------------------------------------------------
   Contributing author: Andres Jaramillo-Botero (Caltech)
------------------------------------------------------------------------- */

#include "atom_vec_electron.h"

#include "atom.h"
#include "citeme.h"

#include <cstring>

using namespace LAMMPS_NS;

static const char cite_user_eff_package[] =
  "EFF package:\n\n"
  "@Article{Jaramillo-Botero11,\n"
  " author = {A. Jaramillo-Botero, J. Su, A. Qi, W. A. Goddard III},\n"
  " title = {Large-Scale, Long-Term Nonadiabatic Electron Molecular Dynamics for Describing Material Properties and Phenomena in Extreme Environments},\n"
  " journal = {J.~Comp.~Chem.},\n"
  " year =    2011,\n"
  " volume =  32,\n"
  " pages =   {497--512}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

AtomVecElectron::AtomVecElectron(LAMMPS *lmp) : AtomVec(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_user_eff_package);

  mass_type = PER_TYPE;
  molecular = Atom::ATOMIC;
  forceclearflag = 1;

  atom->electron_flag = 1;
  atom->q_flag = atom->spin_flag = atom->eradius_flag =
    atom->ervel_flag = atom->erforce_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = (char *) "q spin eradius ervel erforce";
  fields_copy = (char *) "q spin eradius ervel";
  fields_comm = (char *) "eradius";
  fields_comm_vel = (char *) "eradius";
  fields_reverse = (char *) "erforce";
  fields_border = (char *) "q spin eradius";
  fields_border_vel = (char *) "q spin eradius";
  fields_exchange = (char *) "q spin eradius ervel";
  fields_restart = (char *) "q spin eradius ervel";
  fields_create = (char *) "q spin eradius ervel";
  fields_data_atom = (char *) "id type q spin eradius x";
  fields_data_vel = (char *) "id v ervel";

  setup_fields();
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecElectron::grow_pointers()
{
  spin = atom->spin;
  eradius = atom->eradius;
  ervel = atom->ervel;
  erforce = atom->erforce;
}

/* ----------------------------------------------------------------------
   clear extra forces starting at atom N
   nbytes = # of bytes to clear for a per-atom vector
------------------------------------------------------------------------- */

void AtomVecElectron::force_clear(int n, size_t nbytes)
{
  memset(&erforce[n],0,nbytes);
}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecElectron::create_atom_post(int ilocal)
{
  spin[ilocal] = 1;
  eradius[ilocal] = 1.0;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecElectron::data_atom_post(int ilocal)
{
  ervel[ilocal] = 0.0;
}

/* ----------------------------------------------------------------------
   assign an index to named atom property and return index
   return -1 if name is unknown to this atom style
------------------------------------------------------------------------- */

int AtomVecElectron::property_atom(char *name)
{
  if (strcmp(name,"spin") == 0) return 0;
  if (strcmp(name,"eradius") == 0) return 1;
  if (strcmp(name,"ervel") == 0) return 2;
  if (strcmp(name,"erforce") == 0) return 3;
  return -1;
}

/* ----------------------------------------------------------------------
   pack per-atom data into buf for ComputePropertyAtom
   index maps to data specific to this atom style
------------------------------------------------------------------------- */

void AtomVecElectron::pack_property_atom(int index, double *buf,
                                         int nvalues, int groupbit)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int n = 0;

  if (index == 0) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = spin[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 1) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = eradius[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 2) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = ervel[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 3) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = erforce[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  }
}
