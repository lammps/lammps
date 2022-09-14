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

using namespace LAMMPS_NS;

static const char cite_user_eff_package[] =
    "EFF package: doi:10.1002/jcc.21637\n\n"
    "@Article{Jaramillo-Botero11,\n"
    " author = {A. Jaramillo-Botero and J. Su and A. Qi and Goddard, III, W. A.},\n"
    " title = {Large-Scale, Long-Term Nonadiabatic Electron Molecular Dynamics for Describing "
    "Material Properties and Phenomena in Extreme Environments},\n"
    " journal = {J.~Comp.\\ Chem.},\n"
    " year =    2011,\n"
    " volume =  32,\n"
    " number =  3,\n"
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
  atom->q_flag = atom->spin_flag = atom->eradius_flag = atom->ervel_flag = atom->erforce_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = {"q", "espin", "eradius", "ervel", "erforce"};
  fields_copy = {"q", "espin", "eradius", "ervel"};
  fields_comm = {"eradius"};
  fields_comm_vel = {"eradius"};
  fields_reverse = {"erforce"};
  fields_border = {"q", "espin", "eradius"};
  fields_border_vel = {"q", "espin", "eradius"};
  fields_exchange = {"q", "espin", "eradius", "ervel"};
  fields_restart = {"q", "espin", "eradius", "ervel"};
  fields_create = {"q", "espin", "eradius", "ervel"};
  fields_data_atom = {"id", "type", "q", "espin", "eradius", "x"};
  fields_data_vel = {"id", "v", "ervel"};

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
  memset(&erforce[n], 0, nbytes);
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

int AtomVecElectron::property_atom(const std::string &name)
{
  if (name == "espin") return 0;
  if (name == "spin") return 0;    // backward compatibility
  if (name == "eradius") return 1;
  if (name == "ervel") return 2;
  if (name == "erforce") return 3;
  return -1;
}

/* ----------------------------------------------------------------------
   pack per-atom data into buf for ComputePropertyAtom
   index maps to data specific to this atom style
------------------------------------------------------------------------- */

void AtomVecElectron::pack_property_atom(int index, double *buf, int nvalues, int groupbit)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int n = 0;

  if (index == 0) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = spin[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 1) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = eradius[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 2) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = ervel[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 3) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = erforce[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  }
}
