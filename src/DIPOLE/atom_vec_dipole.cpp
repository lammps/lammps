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

#include "atom_vec_dipole.h"

#include "atom.h"
#include "domain.h"
#include "memory.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecDipole::AtomVecDipole(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::ATOMIC;
  mass_type = PER_TYPE;

  atom->q_flag = atom->mu_flag = 1;

  mu_hold = nullptr;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = {"q", "mu"};
  fields_copy = {"q", "mu"};
  fields_comm = {"mu3"};
  fields_comm_vel = {"mu3"};
  fields_border = {"q", "mu"};
  fields_border_vel = {"q", "mu"};
  fields_exchange = {"q", "mu"};
  fields_restart = {"q", "mu"};
  fields_create = {"q", "mu"};
  fields_data_atom = {"id", "type", "q", "x", "mu3"};
  fields_data_vel = {"id", "v"};

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
  mu_one[3] = sqrt(mu_one[0] * mu_one[0] + mu_one[1] * mu_one[1] + mu_one[2] * mu_one[2]);
}

/* ----------------------------------------------------------------------
   convert read_data file info from general to restricted triclinic
   parent class operates on data from Velocities section of data file
   child class operates on dipole moment mu
------------------------------------------------------------------------- */

void AtomVecDipole::read_data_general_to_restricted(int nlocal_previous, int nlocal)
{
  AtomVec::read_data_general_to_restricted(nlocal_previous, nlocal);

  for (int i = nlocal_previous; i < nlocal; i++)
    domain->general_to_restricted_vector(mu[i]);
}

/* ----------------------------------------------------------------------
   convert info output by write_data from restricted to general triclinic
   parent class operates on x and data from Velocities section of data file
   child class operates on dipole momemt mu which has 4 values per atom
------------------------------------------------------------------------- */

void AtomVecDipole::write_data_restricted_to_general()
{
  AtomVec::write_data_restricted_to_general();

  int nlocal = atom->nlocal;
  memory->create(mu_hold,nlocal,3,"atomvec:mu_hold");
  for (int i = 0; i < nlocal; i++) {
    memcpy(&mu_hold[i],&mu[i],3*sizeof(double));
    domain->restricted_to_general_vector(mu[i]);
  }
}

/* ----------------------------------------------------------------------
   restore info output by write_data to restricted triclinic
   original data is in "hold" arrays
   parent class operates on x and data from Velocities section of data file
   child class operates on dipole moment mu which has 4 values per atom
------------------------------------------------------------------------- */

void AtomVecDipole::write_data_restore_restricted()
{
  AtomVec::write_data_restore_restricted();

  if (!mu_hold) return;

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++)
    memcpy(&mu[i],&mu_hold[i],3*sizeof(double));
  memory->destroy(mu_hold);
  mu_hold = nullptr;
}
