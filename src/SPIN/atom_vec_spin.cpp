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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)

   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018).
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics
   and molecular dynamics. Journal of Computational Physics.
------------------------------------------------------------------------- */

#include "atom_vec_spin.h"
#include <cmath>
#include <cstring>
#include "atom.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecSpin::AtomVecSpin(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;
  mass_type = 1;
  forceclearflag = 1;

  atom->sp_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = (char *) "sp fm fm_long";
  fields_copy = (char *) "sp";
  fields_comm = (char *) "sp";
  fields_comm_vel = (char *) "sp";
  fields_reverse = (char *) "fm fm_long";
  fields_border = (char *) "sp";
  fields_border_vel = (char *) "sp";
  fields_exchange = (char *) "sp";
  fields_restart = (char *) "sp";
  fields_create = (char *) "sp";
  fields_data_atom = (char *) "id type x sp";
  fields_data_vel = (char *) "id v omega";

  setup_fields();
}

/* ----------------------------------------------------------------------
   clear all forces (mechanical and magnetic)
------------------------------------------------------------------------- */

void AtomVecSpin::force_clear(int /*n*/, size_t nbytes)
{
  memset(&atom->f[0][0],0,3*nbytes);
  memset(&atom->fm[0][0],0,3*nbytes);
  memset(&atom->fm_long[0][0],0,3*nbytes);
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecSpin::data_atom_post(int ilocal)
{
  double *sp = atom->sp[ilocal];
  double inorm = 1.0/sqrt(sp[0]*sp[0] + sp[1]*sp[1] + sp[2]*sp[2]);
  sp[0] *= inorm;
  sp[1] *= inorm;
  sp[2] *= inorm;
}
