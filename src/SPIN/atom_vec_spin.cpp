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
  molecular = Atom::ATOMIC;
  mass_type = PER_TYPE;
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
  fields_data_vel = (char *) "id v";

  setup_fields();
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecSpin::grow_pointers()
{
  sp = atom->sp;
  fm = atom->fm;
  fm_long = atom->fm_long;
}

/* ----------------------------------------------------------------------
   clear extra forces starting at atom N
   nbytes = # of bytes to clear for a per-atom vector
   include f b/c this is invoked from within SPIN pair styles
------------------------------------------------------------------------- */

void AtomVecSpin::force_clear(int n, size_t nbytes)
{
  memset(&fm[n][0],0,3*nbytes);
  memset(&fm_long[n][0],0,3*nbytes);
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecSpin::data_atom_post(int ilocal)
{
  double *sp_one = sp[ilocal];
  double norm =
    1.0/sqrt(sp_one[0]*sp_one[0] + sp_one[1]*sp_one[1] + sp_one[2]*sp_one[2]);
  sp_one[0] *= norm;
  sp_one[1] *= norm;
  sp_one[2] *= norm;
}
