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

#include <cstring>
#include "fix_external.h"
#include "atom.h"
#include "update.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{PF_CALLBACK,PF_ARRAY};

/* ---------------------------------------------------------------------- */

FixExternal::FixExternal(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  fexternal(NULL), caller_vector(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal fix external command");

  scalar_flag = 1;
  global_freq = 1;
  virial_flag = 1;
  thermo_virial = 1;
  extscalar = 1;

  if (strcmp(arg[3],"pf/callback") == 0) {
    if (narg != 6) error->all(FLERR,"Illegal fix external command");
    mode = PF_CALLBACK;
    ncall = force->inumeric(FLERR,arg[4]);
    napply = force->inumeric(FLERR,arg[5]);
    if (ncall <= 0 || napply <= 0)
      error->all(FLERR,"Illegal fix external command");
  } else if (strcmp(arg[3],"pf/array") == 0) {
    if (narg != 5) error->all(FLERR,"Illegal fix external command");
    mode = PF_ARRAY;
    napply = force->inumeric(FLERR,arg[4]);
    if (napply <= 0) error->all(FLERR,"Illegal fix external command");
  } else error->all(FLERR,"Illegal fix external command");

  callback = NULL;

  // perform initial allocation of atom-based array
  // register with Atom class

  grow_arrays(atom->nmax);
  atom->add_callback(0);

  user_energy = 0.0;

  // optional vector of values provided by caller
  // vector_flag and size_vector are setup via set_vector_length()

  caller_vector = NULL;
}

/* ---------------------------------------------------------------------- */

FixExternal::~FixExternal()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  memory->destroy(fexternal);
  delete [] caller_vector;
}

/* ---------------------------------------------------------------------- */

int FixExternal::setmask()
{
  int mask = 0;
  if (mode == PF_CALLBACK || mode == PF_ARRAY) {
    mask |= PRE_REVERSE;
    mask |= POST_FORCE;
    mask |= THERMO_ENERGY;
    mask |= MIN_POST_FORCE;
  }
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixExternal::init()
{
  if (mode == PF_CALLBACK && callback == NULL)
    error->all(FLERR,"Fix external callback function not set");
}

/* ---------------------------------------------------------------------- */

void FixExternal::setup(int vflag)
{
  post_force(vflag);
}

/* --------------------------------------------------------------------- */

void FixExternal::setup_pre_reverse(int eflag, int vflag)
{
  pre_reverse(eflag,vflag);
}

/* ---------------------------------------------------------------------- */

void FixExternal::min_setup(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   store eflag, so can use it in post_force to tally per-atom energies
------------------------------------------------------------------------- */

void FixExternal::pre_reverse(int eflag, int /*vflag*/)
{
  eflag_caller = eflag;
}

/* ---------------------------------------------------------------------- */

void FixExternal::post_force(int vflag)
{
  bigint ntimestep = update->ntimestep;

  int eflag = eflag_caller;
  ev_init(eflag,vflag);

  // invoke the callback in driver program
  // it will fill fexternal with forces

  if (mode == PF_CALLBACK && ntimestep % ncall == 0)
    (this->callback)(ptr_caller,update->ntimestep,
                     atom->nlocal,atom->tag,atom->x,fexternal);

  // add forces from fexternal to atoms in group

  if (ntimestep % napply == 0) {
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        f[i][0] += fexternal[i][0];
        f[i][1] += fexternal[i][1];
        f[i][2] += fexternal[i][2];
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixExternal::min_post_force(int vflag)
{
  post_force(vflag);
}

// ----------------------------------------------------------------------
// "set" methods caller can invoke directly
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   caller invokes this method to set its contribution to global energy
   unlike other energy/virial set methods:
     do not just return if eflag_global is not set
     b/c input script could access this quantity via compute_scalar()
     even if eflag is not set on a particular timestep
------------------------------------------------------------------------- */

void FixExternal::set_energy_global(double caller_energy)
{
  user_energy = caller_energy;
}

/* ----------------------------------------------------------------------
   caller invokes this method to set its contribution to global virial
------------------------------------------------------------------------- */

void FixExternal::set_virial_global(double *caller_virial)
{
  if (!evflag) return;
  if (!vflag_global) return;

  for (int i = 0; i < 6; i++)
    virial[i] = caller_virial[i];
}

/* ----------------------------------------------------------------------
   caller invokes this method to set its contribution to peratom energy
------------------------------------------------------------------------- */

void FixExternal::set_energy_peratom(double *caller_energy)
{
  if (!eflag_atom) return;

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++)
    eatom[i] = caller_energy[i];
}

/* ----------------------------------------------------------------------
   caller invokes this method to set its contribution to peratom virial
------------------------------------------------------------------------- */

void FixExternal::set_virial_peratom(double **caller_virial)
{
  int i,j;

  if (!evflag) return;
  if (!vflag_atom) return;

  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++)
    for (j = 0; j < 6; j++)
      vatom[i][j] = caller_virial[i][j];
}

/* ----------------------------------------------------------------------
   caller invokes this method to set length of vector of values
   assume all vector values are extensive, could make this an option
------------------------------------------------------------------------- */

void FixExternal::set_vector_length(int n)
{
  delete [] caller_vector;

  vector_flag = 1;
  size_vector = n;
  extvector = 1;

  caller_vector = new double[n];
}

/* ----------------------------------------------------------------------
   caller invokes this method to set Index value in vector
   index ranges from 1 to N inclusive
------------------------------------------------------------------------- */

void FixExternal::set_vector(int index, double value)
{
  if (index >= size_vector)
    error->all(FLERR,"Invalid set_vector index in fix external");
  caller_vector[index-1] = value;
}

/* ----------------------------------------------------------------------
   potential energy of added force
   up to user to set it via set_energy()
------------------------------------------------------------------------- */

double FixExternal::compute_scalar()
{
  return user_energy;
}

/* ----------------------------------------------------------------------
   arbitrary value computed by caller
   up to user to set it via set_vector()
------------------------------------------------------------------------- */

double FixExternal::compute_vector(int n)
{
  return caller_vector[n];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixExternal::memory_usage()
{
  double bytes = 3*atom->nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixExternal::grow_arrays(int nmax)
{
  memory->grow(fexternal,nmax,3,"external:fexternal");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixExternal::copy_arrays(int i, int j, int /*delflag*/)
{
  fexternal[j][0] = fexternal[i][0];
  fexternal[j][1] = fexternal[i][1];
  fexternal[j][2] = fexternal[i][2];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixExternal::pack_exchange(int i, double *buf)
{
  buf[0] = fexternal[i][0];
  buf[1] = fexternal[i][1];
  buf[2] = fexternal[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixExternal::unpack_exchange(int nlocal, double *buf)
{
  fexternal[nlocal][0] = buf[0];
  fexternal[nlocal][1] = buf[1];
  fexternal[nlocal][2] = buf[2];
  return 3;
}

/* ----------------------------------------------------------------------
   external caller sets a callback function to invoke in post_force()
------------------------------------------------------------------------- */

void FixExternal::set_callback(FnPtr caller_callback, void *caller_ptr)
{
  callback = caller_callback;
  ptr_caller = caller_ptr;
}
