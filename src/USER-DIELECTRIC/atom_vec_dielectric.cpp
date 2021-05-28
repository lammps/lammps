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

#include <cmath>
#include <cstdlib>
#include "atom_vec_dielectric.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecDielectric::AtomVecDielectric(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::ATOMIC;
  mass_type = PER_TYPE;
  
  atom->q_flag = atom->mu_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = (char *) "q mu3 area ed em epsilon curvature";
  fields_copy = (char *) "q mu3 area ed em epsilon curvature";
  fields_comm = (char *) "q mu3 area ed em epsilon curvature";
  fields_comm_vel = (char *) "q mu3 area ed em epsilon curvature";
  fields_reverse = (char *) "";
  fields_border = (char *) "q mu3 area ed em epsilon curvature";
  fields_border_vel = (char *) "q mu3 area ed em epsilon curvature";
  fields_exchange = (char *) "q mu3 area ed em epsilon curvature";
  fields_restart = (char * ) "q mu3 area ed em epsilon curvature";
  fields_create = (char *) "q mu3 area ed em epsilon curvature";
  fields_data_atom = (char *) "id molecule type q x mu3 area ed em epsilon curvature";
  fields_data_vel = (char *) "id v";

  setup_fields();
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecDielectric::grow_pointers()
{
  mu = atom->mu;
  area = atom->area;
  ed = atom->ed;
  em = atom->em;
  epsilon = atom->epsilon;
  curvature = atom->curvature;
  q_unscaled = atom->q_unscaled;
}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecDielectric::create_atom_post(int ilocal)
{
  area[ilocal] = 1.0;
  em[ilocal] = 1.0;
  epsilon[ilocal] = 1.0;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecDielectric::data_atom_post(int ilocal)
{
  double* q = atom->q;
  q_unscaled[ilocal] = q[ilocal];
  q[ilocal] /= epsilon[ilocal];
}

/* ----------------------------------------------------------------------
   unmodify values packed by AtomVec::pack_data()
------------------------------------------------------------------------- */

void AtomVecDielectric::pack_data_post(int ilocal)
{
  double* q = atom->q;
  q_unscaled[ilocal] = q[ilocal]/epsilon[ilocal];
}




