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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)

   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018).
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics
   and molecular dynamics. Journal of Computational Physics.
------------------------------------------------------------------------- */

#include "atom_vec_spin.h"

#include "atom.h"
#include "domain.h"
#include "memory.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecSpin::AtomVecSpin(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::ATOMIC;
  mass_type = PER_TYPE;
  forceclearflag = 1;

  atom->sp_flag = 1;

  sp_hold = nullptr;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = {"sp", "fm", "fm_long"};
  fields_copy = {"sp"};
  fields_comm = {"sp"};
  fields_comm_vel = {"sp"};
  fields_reverse = {"fm", "fm_long"};
  fields_border = {"sp"};
  fields_border_vel = {"sp"};
  fields_exchange = {"sp"};
  fields_restart = {"sp"};
  fields_create = {"sp"};
  fields_data_atom = {"id", "type", "x", "sp"};
  fields_data_vel = {"id", "v"};

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
  memset(&fm[n][0], 0, 3 * nbytes);
  memset(&fm_long[n][0], 0, 3 * nbytes);
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecSpin::data_atom_post(int ilocal)
{
  double *sp_one = sp[ilocal];
  double norm = 1.0 / sqrt(sp_one[0] * sp_one[0] + sp_one[1] * sp_one[1] + sp_one[2] * sp_one[2]);
  sp_one[0] *= norm;
  sp_one[1] *= norm;
  sp_one[2] *= norm;
}

/* ----------------------------------------------------------------------
   convert read_data file info from general to restricted triclinic
   parent class operates on data from Velocities section of data file
   child class operates on spin vector sp
------------------------------------------------------------------------- */

void AtomVecSpin::read_data_general_to_restricted(int nlocal_previous, int nlocal)
{
  AtomVec::read_data_general_to_restricted(nlocal_previous, nlocal);

  for (int i = nlocal_previous; i < nlocal; i++)
    domain->general_to_restricted_vector(sp[i]);
}

/* ----------------------------------------------------------------------
   convert info output by write_data from restricted to general triclinic
   parent class operates on x and data from Velocities section of data file
   child class operates on spin vector sp which has 4 values per atom
------------------------------------------------------------------------- */

void AtomVecSpin::write_data_restricted_to_general()
{
  AtomVec::write_data_restricted_to_general();

  int nlocal = atom->nlocal;
  memory->create(sp_hold,nlocal,3,"atomvec:sp_hold");
  for (int i = 0; i < nlocal; i++) {
    memcpy(&sp_hold[i],&sp[i],3*sizeof(double));
    domain->restricted_to_general_vector(sp[i]);
  }
}

/* ----------------------------------------------------------------------
   restore info output by write_data to restricted triclinic
   original data is in "hold" arrays
   parent class operates on x and data from Velocities section of data file
   child class operates on spin vector sp which has 4 values per atom
------------------------------------------------------------------------- */

void AtomVecSpin::write_data_restore_restricted()
{
  AtomVec::write_data_restore_restricted();

  if (!sp_hold) return;

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++)
    memcpy(&sp[i],&sp_hold[i],3*sizeof(double));
  memory->destroy(sp_hold);
  sp_hold = nullptr;
}
