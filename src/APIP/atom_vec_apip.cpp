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

#include "atom_vec_apip.h"

#include "atom.h"

#include <cstring>
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecApip::AtomVecApip(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::ATOMIC;
  mass_type = PER_TYPE;
  forceclearflag = 1;

  atom->apip_lambda_flag = 1;
  atom->apip_lambda_input_flag = 1;
  atom->apip_lambda_input_ta_flag = 1;
  atom->apip_lambda_const_flag = 1;
  atom->apip_lambda_required_flag = 1;
  atom->apip_e_fast_flag = 1;
  atom->apip_e_precise_flag = 1;
  atom->apip_f_const_lambda_flag = 1;
  atom->apip_f_dyn_lambda_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  // The full list of fields is in atom_vec.cpp
  fields_copy = {"apip_lambda", "apip_lambda_required", "apip_lambda_input", "apip_lambda_input_ta",
                 "apip_lambda_const"};
  fields_comm = {"apip_lambda", "apip_lambda_required", "apip_lambda_input_ta",
                 "apip_lambda_const"};
  fields_comm_vel = {};
  fields_border = {"apip_lambda", "apip_lambda_required", "apip_lambda_input_ta",
                   "apip_lambda_const"};
  fields_border_vel = {};
  fields_exchange = {"apip_lambda", "apip_lambda_required", "apip_lambda_input_ta",
                     "apip_lambda_const"};
  fields_restart = {"apip_lambda", "apip_lambda_required", "apip_lambda_input",
                    "apip_lambda_input_ta", "apip_lambda_const"};
  fields_create = {};
  fields_grow = {
      "apip_lambda",          "apip_lambda_required", "apip_lambda_input",
      "apip_lambda_input_ta", "apip_lambda_const",    "apip_e_fast",
      "apip_e_precise",       "apip_f_const_lambda",  "apip_f_dyn_lambda"};    // allocates memory
  fields_reverse = {"apip_f_const_lambda",
                    "apip_f_dyn_lambda"};    // communication of force after calculation
  fields_data_atom = {"id", "type", "x"};
  fields_data_vel = {"id", "v"};

  setup_fields();
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecApip::grow_pointers()
{
  apip_lambda = atom->apip_lambda;
  apip_lambda_required = atom->apip_lambda_required;
  apip_lambda_input = atom->apip_lambda_input;
  apip_lambda_input_ta = atom->apip_lambda_input_ta;
  apip_lambda_const = atom->apip_lambda_const;
  apip_e_fast = atom->apip_e_fast;
  apip_e_precise = atom->apip_e_precise;
  apip_f_const_lambda = atom->apip_f_const_lambda;
  apip_f_dyn_lambda = atom->apip_f_dyn_lambda;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecApip::data_atom_post(int ilocal)
{
  apip_lambda[ilocal] = 0;
  apip_lambda_const[ilocal] = 0;
  apip_lambda_required[ilocal] = ApipLambdaRequired::UNKNOWN;
  apip_lambda_input[ilocal] = 0;
  apip_lambda_input_ta[ilocal] = 0;
  apip_e_fast[ilocal] = 0;
  apip_e_precise[ilocal] = 0;
  apip_f_const_lambda[ilocal][0] = 0;
  apip_f_const_lambda[ilocal][1] = 0;
  apip_f_const_lambda[ilocal][2] = 0;
  apip_f_dyn_lambda[ilocal][0] = 0;
  apip_f_dyn_lambda[ilocal][1] = 0;
  apip_f_dyn_lambda[ilocal][2] = 0;
}

/* ----------------------------------------------------------------------
   clear extra forces starting at atom n
   natoms = # of atoms to clear
   nbytes = natoms * sizeof(double)
   requires forceclearflag = 1 to be called
------------------------------------------------------------------------- */

void AtomVecApip::force_clear(int n, size_t nbytes)
{
  memset(&apip_f_const_lambda[n][0], 0, 3 * nbytes);
  memset(&apip_f_dyn_lambda[n][0], 0, 3 * nbytes);
  memset(&apip_lambda_required[n], 0, nbytes / sizeof(double) * sizeof(int));
}
