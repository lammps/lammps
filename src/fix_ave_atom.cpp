// clang-format off
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

#include "fix_ave_atom.h"

#include "arg_info.h"
#include "atom.h"
#include "compute.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAveAtom::FixAveAtom(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg), array(nullptr)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix ave/atom", error);

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  nrepeat = utils::inumeric(FLERR, arg[4], false, lmp);
  peratom_freq = utils::inumeric(FLERR, arg[5], false, lmp);
  time_depend = 1;

  // expand args if any have wildcard character "*"
  // this can reset nvalues

  int expand = 0;
  char **earg;
  int nvalues = utils::expand_args(FLERR, narg - 6, &arg[6], 1, earg, lmp);

  if (earg != &arg[6]) expand = 1;
  arg = earg;

  // parse values

  values.clear();
  for (int i = 0; i < nvalues; i++) {

    value_t val;
    val.id = "";
    val.val.c = nullptr;

    if (strcmp(arg[i], "x") == 0) {
      val.which = ArgInfo::X;
      val.argindex = 0;
    } else if (strcmp(arg[i], "y") == 0) {
      val.which = ArgInfo::X;
      val.argindex = 1;
    } else if (strcmp(arg[i], "z") == 0) {
      val.which = ArgInfo::X;
      val.argindex = 2;

    } else if (strcmp(arg[i], "vx") == 0) {
      val.which = ArgInfo::V;
      val.argindex = 0;
    } else if (strcmp(arg[i], "vy") == 0) {
      val.which = ArgInfo::V;
      val.argindex = 1;
    } else if (strcmp(arg[i], "vz") == 0) {
      val.which = ArgInfo::V;
      val.argindex = 2;

    } else if (strcmp(arg[i], "fx") == 0) {
      val.which = ArgInfo::F;
      val.argindex = 0;
    } else if (strcmp(arg[i], "fy") == 0) {
      val.which = ArgInfo::F;
      val.argindex = 1;
    } else if (strcmp(arg[i], "fz") == 0) {
      val.which = ArgInfo::F;
      val.argindex = 2;

    } else {
      ArgInfo argi(arg[i]);

      val.which = argi.get_type();
      val.argindex = argi.get_index1();
      val.id = argi.get_name();

      if ((val.which == ArgInfo::UNKNOWN) || (val.which == ArgInfo::NONE) || (argi.get_dim() > 1))
        error->all(FLERR, "Invalid fix ave/atom argument: {}", arg[i]);
    }
    values.push_back(val);
  }

  // if wildcard expansion occurred, free earg memory from exapnd_args()

  if (expand) {
    for (int i = 0; i < nvalues; i++) delete[] earg[i];
    memory->sfree(earg);
  }

  // setup and error check
  // for fix inputs, check that fix frequency is acceptable

  if (nevery <= 0) error->all(FLERR,"Illegal fix ave/atom nevery value: {}", nevery);
  if (nrepeat <= 0) error->all(FLERR,"Illegal fix ave/atom nrepeat value: {}", nrepeat);
  if (peratom_freq <= 0) error->all(FLERR,"Illegal fix ave/atom nfreq value: {}", peratom_freq);
  if (peratom_freq % nevery || nrepeat*nevery > peratom_freq)
    error->all(FLERR,"Inconsistent fix ave/atom nevery/nrepeat/nfreq values");

  for (auto &val : values) {

    if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c) error->all(FLERR,"Compute ID {} for fix ave/atom does not exist", val.id);
      if (val.val.c->peratom_flag == 0)
        error->all(FLERR, "Fix ave/atom compute {} does not calculate per-atom values", val.id);
      if (val.argindex == 0 && val.val.c->size_peratom_cols != 0)
        error->all(FLERR,"Fix ave/atom compute {} does not calculate a per-atom vector", val.id);
      if (val.argindex && val.val.c->size_peratom_cols == 0)
        error->all(FLERR,"Fix ave/atom compute {} does not calculate a per-atom array", val.id);
      if (val.argindex && val.argindex > val.val.c->size_peratom_cols)
        error->all(FLERR,"Fix ave/atom compute {} array is accessed out-of-range", val.id);

    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f) error->all(FLERR, "Fix ID {} for fix ave/atom does not exist", val.id);
      if (val.val.f->peratom_flag == 0)
        error->all(FLERR, "Fix ave/atom fix {} does not calculate per-atom values", val.id);
      if (val.argindex == 0 && val.val.f->size_peratom_cols != 0)
        error->all(FLERR, "Fix ave/atom fix {} does not calculate a per-atom vector", val.id);
      if (val.argindex && val.val.f->size_peratom_cols == 0)
        error->all(FLERR, "Fix ave/atom fix {} does not calculate a per-atom array", val.id);
      if (val.argindex && val.argindex > val.val.f->size_peratom_cols)
        error->all(FLERR,"Fix ave/atom fix {} array is accessed out-of-range", val.id);
      if (nevery % val.val.f->peratom_freq)
        error->all(FLERR, "Fix {} for fix ave/atom not computed at compatible time", val.id);

    } else if (val.which == ArgInfo::VARIABLE) {
      val.val.v = input->variable->find(val.id.c_str());
      if (val.val.v < 0)
        error->all(FLERR,"Variable name {} for fix ave/atom does not exist", val.id);
      if (input->variable->atomstyle(val.val.v) == 0)
        error->all(FLERR,"Fix ave/atom variable {} is not atom-style variable", val.id);
    }
  }

  // this fix produces either a per-atom vector or array

  peratom_flag = 1;
  if (values.size() == 1) size_peratom_cols = 0;
  else size_peratom_cols = values.size();

  // perform initial allocation of atom-based array
  // register with Atom class

  FixAveAtom::grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);

  // zero the array since dump may access it on timestep 0
  // zero the array since a variable may access it before first run

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++)
    for (std::size_t m = 0; m < values.size(); m++)
      array[i][m] = 0.0;

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  irepeat = 0;
  nvalid_last = -1;
  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveAtom::~FixAveAtom()
{
  // unregister callback to this fix from Atom class

  atom->delete_callback(id,Atom::GROW);
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

int FixAveAtom::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveAtom::init()
{
  // set indices and check validity of all computes,fixes,variables

  for (auto &val : values) {
    if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c) error->all(FLERR, "Compute ID {} for fix ave/atom does not exist", val.id);

    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f) error->all(FLERR, "Fix ID {} for fix ave/atom does not exist", val.id);

    } else if (val.which == ArgInfo::VARIABLE) {
      val.val.v = input->variable->find(val.id.c_str());
      if (val.val.v < 0)
        error->all(FLERR,"Variable name {} for fix ave/atom does not exist", val.id);
    }
  }

  // need to reset nvalid if nvalid < ntimestep b/c minimize was performed

  if (nvalid < update->ntimestep) {
    irepeat = 0;
    nvalid = nextvalid();
    modify->addstep_compute_all(nvalid);
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveAtom::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveAtom::end_of_step()
{
  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;
  nvalid_last = nvalid;

  // zero if first step

  int nlocal = atom->nlocal;

  if (irepeat == 0)
    for (int i = 0; i < nlocal; i++)
      for (std::size_t m = 0; m < values.size(); m++)
        array[i][m] = 0.0;

  // accumulate results of attributes,computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  int *mask = atom->mask;

  int i, j, m = 0;
  for (auto &val : values) {
    j = val.argindex;

    if (val.which == ArgInfo::X) {
      double **x = atom->x;
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) array[i][m] += x[i][j];

    } else if (val.which == ArgInfo::V) {
      double **v = atom->v;
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) array[i][m] += v[i][j];

    } else if (val.which == ArgInfo::F) {
      double **f = atom->f;
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) array[i][m] += f[i][j];

    // invoke compute if not previously invoked

    } else if (val.which == ArgInfo::COMPUTE) {
      if (!(val.val.c->invoked_flag & Compute::INVOKED_PERATOM)) {
        val.val.c->compute_peratom();
        val.val.c->invoked_flag |= Compute::INVOKED_PERATOM;
      }

      if (j == 0) {
        double *compute_vector = val.val.c->vector_atom;
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit) array[i][m] += compute_vector[i];
      } else {
        int jm1 = j - 1;
        double **compute_array = val.val.c->array_atom;
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit) array[i][m] += compute_array[i][jm1];
      }

    // access fix fields, guaranteed to be ready

    } else if (val.which == ArgInfo::FIX) {
      if (j == 0) {
        double *fix_vector = val.val.f->vector_atom;
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit) array[i][m] += fix_vector[i];
      } else {
        int jm1 = j - 1;
        double **fix_array = val.val.f->array_atom;
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit) array[i][m] += fix_array[i][jm1];
      }

    // evaluate atom-style variable
    // final argument = 1 sums result to array

    } else if (val.which == ArgInfo::VARIABLE) {
      if (array) input->variable->compute_atom(val.val.v,igroup,&array[0][m],values.size(),1);
      else input->variable->compute_atom(val.val.v,igroup,nullptr,values.size(),1);
    }
    ++m;
  }

  // done if irepeat < nrepeat
  // else reset irepeat and nvalid

  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    modify->addstep_compute(nvalid);
    return;
  }

  irepeat = 0;
  nvalid = ntimestep+peratom_freq - ((bigint)nrepeat-1)*nevery;
  modify->addstep_compute(nvalid);

  if (array == nullptr) return;

  // average the final result for the Nfreq timestep

  double repeat = nrepeat;
  for (i = 0; i < nlocal; i++)
    for (m = 0; m < (int)values.size(); m++)
      array[i][m] /= repeat;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixAveAtom::memory_usage()
{
  double bytes;
  bytes = (double)atom->nmax*values.size() * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixAveAtom::grow_arrays(int nmax)
{
  memory->grow(array,nmax,values.size(),"fix_ave/atom:array");
  array_atom = array;
  if (array) vector_atom = array[0];
  else vector_atom = nullptr;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixAveAtom::copy_arrays(int i, int j, int /*delflag*/)
{
  for (std::size_t m = 0; m < values.size(); m++)
    array[j][m] = array[i][m];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixAveAtom::pack_exchange(int i, double *buf)
{
  for (std::size_t m = 0; m < values.size(); m++) buf[m] = array[i][m];
  return values.size();
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixAveAtom::unpack_exchange(int nlocal, double *buf)
{
  for (std::size_t m = 0; m < values.size(); m++) array[nlocal][m] = buf[m];
  return values.size();
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
------------------------------------------------------------------------- */

bigint FixAveAtom::nextvalid()
{
  bigint nvalid = (update->ntimestep/peratom_freq)*peratom_freq + peratom_freq;
  if (nvalid-peratom_freq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= ((bigint)nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += peratom_freq;
  return nvalid;
}
