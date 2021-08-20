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

#include "pair_tracker.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "fix_dummy.h"
#include "fix_neigh_history.h"
#include "fix_pair_tracker.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairTracker::PairTracker(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 1;
  no_virial_fdotr_compute = 1;

  neighprev = 0;
  history = 1;
  size_history = 4;
  nondefault_history_transfer = 1;

  finitecutflag = 0;

  // create dummy fix as placeholder for FixNeighHistory
  // this is so final order of Modify:fix will conform to input script

  fix_history = nullptr;
  modify->add_fix("NEIGH_HISTORY_TRACK_DUMMY all DUMMY");
  fix_dummy = (FixDummy *) modify->fix[modify->nfix - 1];
}

/* ---------------------------------------------------------------------- */

PairTracker::~PairTracker()
{
  if (!fix_history)
    modify->delete_fix("NEIGH_HISTORY_TRACK_DUMMY");
  else
    modify->delete_fix("NEIGH_HISTORY_TRACK");

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);

    delete[] onerad_dynamic;
    delete[] onerad_frozen;
    delete[] maxrad_dynamic;
    delete[] maxrad_frozen;
  }
}

/* ---------------------------------------------------------------------- */

void PairTracker::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, time;
  double radi, radj, radsum, rsq, r;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int *touch, **firsttouch;
  double *data, *alldata, **firstdata;

  int updateflag = 1;
  if (update->setupflag) updateflag = 0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = fix_history->firstflag;
  firstdata = fix_history->firstvalue;

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    if (finitecutflag) radi = radius[i];
    itype = type[i];
    touch = firsttouch[i];
    alldata = firstdata[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      r = sqrt(rsq);

      if (finitecutflag) {
        radj = radius[j];
        radsum = radi + radj;

        if (rsq >= radsum * radsum) {

          data = &alldata[size_history * jj];
          if (touch[jj] == 1) {
            fix_pair_tracker->lost_contact(i, j, data[0], data[1], data[2], data[3]);
          }
          touch[jj] = 0;
          data[0] = 0.0;    // initial time
          data[1] = 0.0;    // initial timestep
          data[2] = 0.0;    // sum of r, may overflow
          data[3] = 0.0;    // min of r

        } else {

          data = &alldata[size_history * jj];
          if (touch[jj] == 0) {
            time = update->atime + (update->ntimestep - update->atimestep) * update->dt;
            data[0] = time;
            data[1] = (double) update->ntimestep;
            data[2] = r;
            data[3] = r;
          } else if (updateflag) {
            data[2] += r;
            if (data[3] > r) data[3] = r;
          }
          touch[jj] = 1;
        }
      } else {
        jtype = type[j];
        if (rsq >= cutsq[itype][jtype]) {

          data = &alldata[size_history * jj];
          if (touch[jj] == 1) {
            fix_pair_tracker->lost_contact(i, j, data[0], data[1], data[2], data[3]);
          }

          touch[jj] = 0;
          data[0] = 0.0;    // initial time
          data[1] = 0.0;    // initial timestep
          data[2] = 0.0;    // sum of r, may overflow
          data[3] = 0.0;    // min of r

        } else {

          data = &alldata[size_history * jj];
          if (touch[jj] == 0) {
            time = update->atime + (update->ntimestep - update->atimestep) * update->dt;
            data[0] = time;
            data[1] = (double) update->ntimestep;
            data[2] = r;
            data[3] = r;
          } else if (updateflag) {
            data[2] += r;
            if (data[3] > r) data[3] = r;
          }
          touch[jj] = 1;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairTracker::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");

  onerad_dynamic = new double[n + 1];
  onerad_frozen = new double[n + 1];
  maxrad_dynamic = new double[n + 1];
  maxrad_frozen = new double[n + 1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairTracker::settings(int narg, char **arg)
{
  if (narg != 0 && narg != 1) error->all(FLERR, "Illegal pair_style command");

  if (narg == 1) {
    if (strcmp(arg[0], "finite") == 0)
      finitecutflag = 1;
    else
      error->all(FLERR, "Illegal pair_style command");
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairTracker::coeff(int narg, char **arg)
{
  if (narg > 2 && finitecutflag) error->all(FLERR, "Incorrect args for pair coefficients");
  if (narg != 3 && !finitecutflag) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double cut_one = 0.0;
  if (!finitecutflag) cut_one = utils::numeric(FLERR, arg[2], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      setflag[i][j] = 1;
      cut[i][j] = cut_one;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairTracker::init_style()
{
  int i;

  // error and warning checks

  if (!atom->radius_flag && finitecutflag)
    error->all(FLERR, "Pair tracker requires atom attribute radius for finite cutoffs");

  // need a history neigh list

  int irequest = neighbor->request(this, instance_me);
  if (finitecutflag) {
    neighbor->requests[irequest]->size = 1;
    neighbor->requests[irequest]->history = 1;
    // history flag won't affect results, but match granular pairstyles
    // so neighborlist can be copied to reduce overhead
  }

  // if history is stored and first init, create Fix to store history
  // it replaces FixDummy, created in the constructor
  // this is so its order in the fix list is preserved

  if (fix_history == nullptr) {
    modify->replace_fix("NEIGH_HISTORY_TRACK_DUMMY",
                        fmt::format("NEIGH_HISTORY_TRACK all NEIGH_HISTORY {}", size_history), 1);
    int ifix = modify->find_fix("NEIGH_HISTORY_TRACK");
    fix_history = (FixNeighHistory *) modify->fix[ifix];
    fix_history->pair = this;
    fix_history->use_bit_flag = 0;
  }

  if (finitecutflag) {

    if (force->pair->beyond_contact)
      error->all(FLERR,
                 "Pair tracker incompatible with granular pairstyles that extend beyond contact");
    // check for FixPour and FixDeposit so can extract particle radii

    int ipour;
    for (ipour = 0; ipour < modify->nfix; ipour++)
      if (strcmp(modify->fix[ipour]->style, "pour") == 0) break;
    if (ipour == modify->nfix) ipour = -1;

    int idep;
    for (idep = 0; idep < modify->nfix; idep++)
      if (strcmp(modify->fix[idep]->style, "deposit") == 0) break;
    if (idep == modify->nfix) idep = -1;

    // set maxrad_dynamic and maxrad_frozen for each type
    // include future FixPour and FixDeposit particles as dynamic

    int itype;
    for (i = 1; i <= atom->ntypes; i++) {
      onerad_dynamic[i] = onerad_frozen[i] = 0.0;
      if (ipour >= 0) {
        itype = i;
        onerad_dynamic[i] = *((double *) modify->fix[ipour]->extract("radius", itype));
      }
      if (idep >= 0) {
        itype = i;
        onerad_dynamic[i] = *((double *) modify->fix[idep]->extract("radius", itype));
      }
    }

    double *radius = atom->radius;
    int *mask = atom->mask;
    int *type = atom->type;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++)
      if (mask[i] & freeze_group_bit)
        onerad_frozen[type[i]] = MAX(onerad_frozen[type[i]], radius[i]);
      else
        onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]], radius[i]);

    MPI_Allreduce(&onerad_dynamic[1], &maxrad_dynamic[1], atom->ntypes, MPI_DOUBLE, MPI_MAX, world);
    MPI_Allreduce(&onerad_frozen[1], &maxrad_frozen[1], atom->ntypes, MPI_DOUBLE, MPI_MAX, world);
  }

  int ifix = modify->find_fix("NEIGH_HISTORY_TRACK");
  if (ifix < 0) error->all(FLERR, "Could not find pair fix neigh history ID");
  fix_history = (FixNeighHistory *) modify->fix[ifix];

  ifix = modify->find_fix_by_style("pair/tracker");
  if (ifix < 0) error->all(FLERR, "Cannot use pair tracker without fix pair/tracker");
  fix_pair_tracker = (FixPairTracker *) modify->fix[ifix];
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairTracker::init_one(int i, int j)
{
  if (!allocated) allocate();

  // always mix prefactors geometrically

  if (setflag[i][j] == 0) { cut[i][j] = mix_distance(cut[i][i], cut[j][j]); }

  cut[j][i] = cut[i][j];

  // if finite, cutoff = sum of max I,J radii for
  // dynamic/dynamic & dynamic/frozen interactions, but not frozen/frozen
  double cutoff;
  if (finitecutflag) {
    cutoff = maxrad_dynamic[i] + maxrad_dynamic[j];
    cutoff = MAX(cutoff, maxrad_frozen[i] + maxrad_dynamic[j]);
    cutoff = MAX(cutoff, maxrad_dynamic[i] + maxrad_frozen[j]);
  } else {
    cutoff = cut[i][j];
  }
  return cutoff;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairTracker::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) { fwrite(&cut[i][j], sizeof(double), 1, fp); }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairTracker::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) { utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error); }
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairTracker::write_restart_settings(FILE *fp)
{
  fwrite(&mix_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairTracker::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) { utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error); }
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

double PairTracker::single(int /*i*/, int /*j*/, int /*itype*/, int /*jtype*/, double /*rsq*/,
                           double /*factor_coul*/, double /*factor_lj*/, double &/*fforce*/)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   transfer history during fix/neigh/history exchange
   only needed if any history entries i-j are not just negative of j-i entries
------------------------------------------------------------------------- */

void PairTracker::transfer_history(double *source, double *target)
{
  for (int i = 0; i < size_history; i++) target[i] = source[i];
}

/* ----------------------------------------------------------------------
   self-interaction range of particle if finite particles
------------------------------------------------------------------------- */

double PairTracker::atom2cut(int i)
{
  double cut = atom->radius[i] * 2;
  return cut;
}

/* ----------------------------------------------------------------------
   maximum interaction range for two finite particles
------------------------------------------------------------------------- */

double PairTracker::radii2cut(double r1, double r2)
{
  double cut = r1 + r2;
  return cut;
}
