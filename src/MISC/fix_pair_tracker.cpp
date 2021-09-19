/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_pair_tracker.h"

#include "atom.h"
#include "error.h"
#include "memory.h"
#include "tokenizer.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

#define DELTA 1000

/* ---------------------------------------------------------------------- */

FixPairTracker::FixPairTracker(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), nvalues(0), vector(nullptr), array(nullptr), pack_choice(nullptr)
{
  if (narg < 3) error->all(FLERR, "Illegal fix pair/tracker command");
  local_flag = 1;

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR, "Illegal fix pair/tracker command");
  local_freq = nevery;

  // If optional arguments included, this will be oversized
  nvalues = narg - 4;
  pack_choice = new FnPtrPack[nvalues];

  tmin = -1;
  type_filter = nullptr;
  int iarg = 4;
  nvalues = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "id1") == 0) {
      pack_choice[nvalues++] = &FixPairTracker::pack_id1;
    } else if (strcmp(arg[iarg], "id2") == 0) {
      pack_choice[nvalues++] = &FixPairTracker::pack_id2;

    } else if (strcmp(arg[iarg], "time/created") == 0) {
      pack_choice[nvalues++] = &FixPairTracker::pack_time_created;
    } else if (strcmp(arg[iarg], "time/broken") == 0) {
      pack_choice[nvalues++] = &FixPairTracker::pack_time_broken;
    } else if (strcmp(arg[iarg], "time/total") == 0) {
      pack_choice[nvalues++] = &FixPairTracker::pack_time_total;

    } else if (strcmp(arg[iarg], "x") == 0) {
      pack_choice[nvalues++] = &FixPairTracker::pack_x;
    } else if (strcmp(arg[iarg], "y") == 0) {
      pack_choice[nvalues++] = &FixPairTracker::pack_y;
    } else if (strcmp(arg[iarg], "z") == 0) {
      pack_choice[nvalues++] = &FixPairTracker::pack_z;

    } else if (strcmp(arg[iarg], "r/min") == 0) {
      pack_choice[nvalues++] = &FixPairTracker::pack_rmin;
    } else if (strcmp(arg[iarg], "r/ave") == 0) {
      pack_choice[nvalues++] = &FixPairTracker::pack_rave;

    } else if (strcmp(arg[iarg], "time/min") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "Invalid keyword in fix pair/tracker command");
      tmin = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg++;

    } else if (strcmp(arg[iarg], "type/include") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "Invalid keyword in fix pair/tracker command");
      int ntypes = atom->ntypes;

      int i, j, itype, jtype, in, jn, infield, jnfield;
      int inlo, inhi, jnlo, jnhi;
      char *istr, *jstr;
      if (!type_filter) {
        memory->create(type_filter, ntypes + 1, ntypes + 1, "fix/pair/tracker:type_filter");

        for (i = 0; i <= ntypes; i++) {
          for (j = 0; j <= ntypes; j++) { type_filter[i][j] = 0; }
        }
      }

      auto iwords = Tokenizer(arg[iarg + 1], ",").as_vector();
      auto jwords = Tokenizer(arg[iarg + 2], ",").as_vector();

      for (const auto &ifield : iwords) {
        utils::bounds(FLERR, ifield, 1, ntypes, inlo, inhi, error);
        for (const auto &jfield : jwords) {
          utils::bounds(FLERR, jfield, 1, ntypes, jnlo, jnhi, error);

          for (itype = inlo; itype <= inhi; itype++) {
            for (jtype = jnlo; jtype <= jnhi; jtype++) {
              type_filter[itype][jtype] = 1;
              type_filter[jtype][itype] = 1;
            }
          }
        }
      }
      iarg += 2;

    } else
      error->all(FLERR, "Invalid keyword in fix pair/tracker command");

    iarg++;
  }

  if (nvalues == 1)
    size_local_cols = 0;
  else
    size_local_cols = nvalues;

  nmax = 0;
  ncount = 0;
  vector = nullptr;
  array = nullptr;
}

/* ---------------------------------------------------------------------- */

FixPairTracker::~FixPairTracker()
{
  delete[] pack_choice;

  memory->destroy(vector);
  memory->destroy(array);
  memory->destroy(type_filter);
}

/* ---------------------------------------------------------------------- */

int FixPairTracker::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::init()
{
  // Set size of array/vector
  ncount = 0;

  if (ncount > nmax) reallocate(ncount);

  size_local_rows = ncount;
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::lost_contact(int i, int j, double time_tmp, double nstep_tmp, double rsum_tmp,
                                  double rmin_tmp)
{

  double time = update->atime + (update->ntimestep - update->atimestep) * update->dt;
  if ((time - time_tmp) < tmin) return;

  if (type_filter) {
    int *type = atom->type;
    if (type_filter[type[i]][type[j]] == 0) return;
  }

  int *mask = atom->mask;
  if (!(mask[i] & groupbit)) return;
  if (!(mask[j] & groupbit)) return;

  if (ncount == nmax) reallocate(ncount);

  index_i = i;
  index_j = j;

  rmin = rmin_tmp;
  rsum = rsum_tmp;
  time_initial = time_tmp;
  nstep_initial = nstep_tmp;

  // fill vector or array with local values
  if (nvalues == 1) {
    (this->*pack_choice[0])(0);
  } else {
    for (int k = 0; k < nvalues; k++) { (this->*pack_choice[k])(k); }
  }

  ncount += 1;
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::post_force(int /*vflag*/)
{
  if (update->ntimestep % nevery == 0) {
    size_local_rows = ncount;
    ncount = 0;
  }
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::reallocate(int n)
{
  // grow vector or array
  while (nmax <= n) nmax += DELTA;

  if (nvalues == 1) {
    memory->grow(vector, nmax, "fix_pair_tracker:vector");
    vector_local = vector;
  } else {
    memory->grow(array, nmax, nvalues, "fix_pair_tracker:array");
    array_local = array;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double FixPairTracker::memory_usage()
{
  double bytes = nmax * (double) nvalues * sizeof(double);
  bytes += nmax * 2 * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   one method for every keyword fix pair/tracker can output
   the atom property is packed into a local vector or array
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void FixPairTracker::pack_time_created(int n)
{
  if (nvalues == 1)
    vector[ncount] = time_initial;
  else
    array[ncount][n] = time_initial;
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::pack_time_broken(int n)
{
  double time = update->atime + (update->ntimestep - update->atimestep) * update->dt;
  if (nvalues == 1)
    vector[ncount] = time;
  else
    array[ncount][n] = time;
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::pack_time_total(int n)
{
  double time = update->atime + (update->ntimestep - update->atimestep) * update->dt;
  if (nvalues == 1)
    vector[ncount] = time - time_initial;
  else
    array[ncount][n] = time - time_initial;
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::pack_id1(int n)
{
  tagint *tag = atom->tag;

  if (nvalues == 1)
    vector[ncount] = tag[index_i];
  else
    array[ncount][n] = tag[index_i];
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::pack_id2(int n)
{
  tagint *tag = atom->tag;

  if (nvalues == 1)
    vector[ncount] = tag[index_j];
  else
    array[ncount][n] = tag[index_j];
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::pack_x(int n)
{
  double **x = atom->x;

  if (nvalues == 1)
    vector[ncount] = (x[index_i][0] + x[index_j][0]) / 2;
  else
    array[ncount][n] = (x[index_i][0] + x[index_j][0]) / 2;
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::pack_y(int n)
{
  double **x = atom->x;

  if (nvalues == 1)
    vector[ncount] = (x[index_i][1] + x[index_j][1]) / 2;
  else
    array[ncount][n] = (x[index_i][1] + x[index_j][1]) / 2;
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::pack_z(int n)
{
  double **x = atom->x;

  if (nvalues == 1)
    vector[ncount] = (x[index_i][2] + x[index_j][2]) / 2;
  else
    array[ncount][n] = (x[index_i][2] + x[index_j][2]) / 2;
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::pack_rmin(int n)
{
  if (nvalues == 1)
    vector[ncount] = rmin;
  else
    array[ncount][n] = rmin;
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::pack_rave(int n)
{
  if (nvalues == 1)
    vector[ncount] = rsum / (update->ntimestep - nstep_initial);
  else
    array[ncount][n] = rsum / (update->ntimestep - nstep_initial);
}
