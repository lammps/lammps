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

#include "atom_vec_dipole.h"
#include <cmath>
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecDipole::AtomVecDipole(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;
  mass_type = 1;

  atom->q_flag = atom->mu_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = (char *) "q mu";
  fields_copy = (char *) "q mu";
  fields_comm = (char *) "mu3";
  fields_comm_vel = (char *) "mu3";
  fields_reverse = (char *) "";
  fields_border = (char *) "q mu";
  fields_border_vel = (char *) "q mu";
  fields_exchange = (char *) "q mu";
  fields_restart = (char *) "q mu";
  fields_create = (char *) "q mu";
  fields_data_atom = (char *) "id type q x mu3";
  fields_data_vel = (char *) "id v";

  setup_fields();
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecDipole::grow_pointers()
{
  mu = atom->mu;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecDipole::data_atom_post(int ilocal)
{
  double *mu_one = mu[ilocal];
  mu_one[3] = 
    sqrt(mu_one[0]*mu_one[0] + mu_one[1]*mu_one[1] + mu_one[2]*mu_one[2]);
}
