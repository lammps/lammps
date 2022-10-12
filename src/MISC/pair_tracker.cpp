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
#include "fix_store_local.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "tokenizer.h"
#include "update.h"
#include "utils.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairTracker::PairTracker(LAMMPS *lmp) :
    Pair(lmp), onerad_dynamic(nullptr), onerad_frozen(nullptr), maxrad_dynamic(nullptr),
    maxrad_frozen(nullptr), id_fix_store_local(nullptr), fix_dummy(nullptr), fix_history(nullptr),
    fix_store_local(nullptr), type_filter(nullptr), output_data(nullptr), pack_choice(nullptr)
{
  single_enable = 1;
  no_virial_fdotr_compute = 1;

  neighprev = 0;
  history = 1;
  size_history = 3;
  nondefault_history_transfer = 1;

  finitecutflag = 0;
  tmin = -1;

  // create dummy fix as placeholder for FixNeighHistory
  // this is so final order of Modify:fix will conform to input script
  fix_dummy = dynamic_cast<FixDummy *>(modify->add_fix("NEIGH_HISTORY_TRACK_DUMMY all DUMMY"));
}

/* ---------------------------------------------------------------------- */

PairTracker::~PairTracker()
{
  if (!fix_history)
    modify->delete_fix("NEIGH_HISTORY_TRACK_DUMMY");
  else
    modify->delete_fix("NEIGH_HISTORY_TRACK");
  if (id_fix_store_local) modify->delete_fix(id_fix_store_local);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);

    delete[] onerad_dynamic;
    delete[] onerad_frozen;
    delete[] maxrad_dynamic;
    delete[] maxrad_frozen;
  }

  delete[] pack_choice;
  delete[] id_fix_store_local;

  memory->destroy(output_data);
  memory->destroy(type_filter);
}

/* ---------------------------------------------------------------------- */

void PairTracker::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
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
          if (touch[jj] == 1) process_data(i, j, data);

          touch[jj] = 0;
          data[0] = 0.0;    // initial timestep
          data[1] = 0.0;    // sum of r, may be inaccurate over long times
          data[2] = 0.0;    // min of r

        } else {

          data = &alldata[size_history * jj];
          if (touch[jj] == 0) {
            data[0] = (double) update->ntimestep;
            data[1] = r;
            data[2] = r;
          } else if (updateflag) {
            data[1] += r;
            if (data[2] > r) data[2] = r;
          }
          touch[jj] = 1;
        }
      } else {
        jtype = type[j];
        if (rsq >= cutsq[itype][jtype]) {

          data = &alldata[size_history * jj];
          if (touch[jj] == 1) process_data(i, j, data);

          touch[jj] = 0;
          data[0] = 0.0;    // initial timestep
          data[1] = 0.0;    // sum of r, may be inaccurate over long times
          data[2] = 0.0;    // min of r

        } else {

          data = &alldata[size_history * jj];
          if (touch[jj] == 0) {
            data[0] = (double) update->ntimestep;
            data[1] = r;
            data[2] = r;
          } else if (updateflag) {
            data[1] += r;
            if (data[2] > r) data[2] = r;
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
  if (narg < 2) error->all(FLERR, "Illegal pair_style command");

  id_fix_store_local = utils::strdup(arg[0]);
  store_local_freq = utils::inumeric(FLERR, arg[1], false, lmp);

  // If optional arguments included, this will be oversized
  pack_choice = new FnPtrPack[narg - 1];

  nvalues = 0;
  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "finite") == 0) {
      finitecutflag = 1;
    } else if (strcmp(arg[iarg], "id1") == 0) {
      pack_choice[nvalues++] = &PairTracker::pack_id1;
    } else if (strcmp(arg[iarg], "id2") == 0) {
      pack_choice[nvalues++] = &PairTracker::pack_id2;
    } else if (strcmp(arg[iarg], "time/created") == 0) {
      pack_choice[nvalues++] = &PairTracker::pack_time_created;
    } else if (strcmp(arg[iarg], "time/broken") == 0) {
      pack_choice[nvalues++] = &PairTracker::pack_time_broken;
    } else if (strcmp(arg[iarg], "time/total") == 0) {
      pack_choice[nvalues++] = &PairTracker::pack_time_total;
    } else if (strcmp(arg[iarg], "x") == 0) {
      pack_choice[nvalues++] = &PairTracker::pack_x;
    } else if (strcmp(arg[iarg], "y") == 0) {
      pack_choice[nvalues++] = &PairTracker::pack_y;
    } else if (strcmp(arg[iarg], "z") == 0) {
      pack_choice[nvalues++] = &PairTracker::pack_z;
    } else if (strcmp(arg[iarg], "r/min") == 0) {
      pack_choice[nvalues++] = &PairTracker::pack_rmin;
    } else if (strcmp(arg[iarg], "r/ave") == 0) {
      pack_choice[nvalues++] = &PairTracker::pack_rave;
    } else if (strcmp(arg[iarg], "time/min") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "Invalid keyword in pair tracker command");
      tmin = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg++;

    } else if (strcmp(arg[iarg], "type/include") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "Invalid keyword in pair tracker command");
      int ntypes = atom->ntypes;
      int i, j, itype, jtype;
      int inlo, inhi, jnlo, jnhi;

      if (!type_filter) {
        memory->create(type_filter, ntypes + 1, ntypes + 1, "pair/tracker:type_filter");

        for (i = 0; i <= ntypes; i++) {
          for (j = 0; j <= ntypes; j++) type_filter[i][j] = 0;
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

    } else {
      error->all(FLERR, "Illegal pair_style command");
    }
    iarg++;
  }

  if (nvalues == 0) error->all(FLERR, "Must request at least one value to output");
  memory->create(output_data, nvalues, "pair/tracker:output_data");

  fix_store_local = dynamic_cast<FixStoreLocal *>(modify->get_fix_by_id(id_fix_store_local));
  if (!fix_store_local)
    fix_store_local = dynamic_cast<FixStoreLocal *>(modify->add_fix(
        fmt::format("{} all STORE/LOCAL {} {}", id_fix_store_local, store_local_freq, nvalues)));
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
  // error and warning checks

  if (!atom->radius_flag && finitecutflag)
    error->all(FLERR, "Pair tracker requires atom attribute radius for finite cutoffs");

  int neigh_flags = NeighConst::REQ_DEFAULT;
  // history flag won't affect results, but match granular pairstyles
  // so neighborlist can be copied to reduce overhead
  if (finitecutflag) neigh_flags |= NeighConst::REQ_SIZE | NeighConst::REQ_HISTORY;
  neighbor->add_request(this, neigh_flags);

  // if history is stored and first init, create Fix to store history
  // it replaces FixDummy, created in the constructor
  // this is so its order in the fix list is preserved

  if (fix_history == nullptr) {
    modify->replace_fix("NEIGH_HISTORY_TRACK_DUMMY",
                        fmt::format("NEIGH_HISTORY_TRACK all NEIGH_HISTORY {}", size_history), 1);
    fix_history = dynamic_cast<FixNeighHistory *>(modify->get_fix_by_id("NEIGH_HISTORY_TRACK"));
    fix_history->pair = this;
    fix_history->use_bit_flag = 0;
  } else {
    fix_history = dynamic_cast<FixNeighHistory *>(modify->get_fix_by_id("NEIGH_HISTORY_TRACK"));
    if (!fix_history) error->all(FLERR, "Could not find pair fix neigh history ID");
  }

  if (finitecutflag) {
    if (force->pair->beyond_contact)
      error->all(FLERR,
                 "Pair tracker incompatible with granular pairstyles that extend beyond contact");

    // check for FixFreeze and set freeze_group_bit

    auto fixlist = modify->get_fix_by_style("^freeze");
    if (fixlist.size() == 0)
      freeze_group_bit = 0;
    else if (fixlist.size() > 1)
      error->all(FLERR, "Only one fix freeze command at a time allowed");
    else
      freeze_group_bit = fixlist.front()->groupbit;

    // check for FixPour and FixDeposit so can extract particle radii

    auto pours = modify->get_fix_by_style("^pour");
    auto deps = modify->get_fix_by_style("^deposit");

    // set maxrad_dynamic and maxrad_frozen for each type
    // include future FixPour and FixDeposit particles as dynamic

    int itype = 0;
    for (int i = 1; i <= atom->ntypes; i++) {
      onerad_dynamic[i] = onerad_frozen[i] = 0.0;
      for (auto &ipour : pours) {
        itype = i;
        double maxrad = *((double *) ipour->extract("radius", itype));
        if (maxrad > 0.0) onerad_dynamic[i] = maxrad;
      }
      for (auto &idep : deps) {
        itype = i;
        double maxrad = *((double *) idep->extract("radius", itype));
        if (maxrad > 0.0) onerad_dynamic[i] = maxrad;
      }
    }

    double *radius = atom->radius;
    int *mask = atom->mask;
    int *type = atom->type;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & freeze_group_bit)
        onerad_frozen[type[i]] = MAX(onerad_frozen[type[i]], radius[i]);
      else
        onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]], radius[i]);

    MPI_Allreduce(&onerad_dynamic[1], &maxrad_dynamic[1], atom->ntypes, MPI_DOUBLE, MPI_MAX, world);
    MPI_Allreduce(&onerad_frozen[1], &maxrad_frozen[1], atom->ntypes, MPI_DOUBLE, MPI_MAX, world);
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairTracker::init_one(int i, int j)
{
  if (!allocated) allocate();

  // always mix prefactors geometrically
  if (setflag[i][j] == 0) cut[i][j] = mix_distance(cut[i][i], cut[j][j]);

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
      if (setflag[i][j]) fwrite(&cut[i][j], sizeof(double), 1, fp);
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
        if (me == 0) utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
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
                           double /*factor_coul*/, double /*factor_lj*/, double & /*fforce*/)
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

/* ---------------------------------------------------------------------- */

void PairTracker::process_data(int i, int j, double *input_data)
{
  if ((update->ntimestep - input_data[0]) < tmin) return;

  if (type_filter) {
    int *type = atom->type;
    if (type_filter[type[i]][type[j]] == 0) return;
  }

  for (int k = 0; k < nvalues; k++) (this->*pack_choice[k])(k, i, j, input_data);
  fix_store_local->add_data(output_data, i, j);
}

/* ----------------------------------------------------------------------
   one method for every keyword fix pair/tracker can output
   the atom property is packed into a local vector or array
------------------------------------------------------------------------- */

void PairTracker::pack_time_created(int n, int /*i*/, int /*j*/, double *data)
{
  output_data[n] = data[0];
}

/* ---------------------------------------------------------------------- */

void PairTracker::pack_time_broken(int n, int /*i*/, int /*j*/, double * /*data*/)
{
  output_data[n] = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void PairTracker::pack_time_total(int n, int /*i*/, int /*j*/, double *data)
{
  output_data[n] = update->ntimestep - data[0];
}

/* ---------------------------------------------------------------------- */

void PairTracker::pack_id1(int n, int i, int /*j*/, double * /*data*/)
{
  output_data[n] = atom->tag[i];
}

/* ---------------------------------------------------------------------- */

void PairTracker::pack_id2(int n, int /*i*/, int j, double * /*data*/)
{
  output_data[n] = atom->tag[j];
}

/* ---------------------------------------------------------------------- */

void PairTracker::pack_x(int n, int i, int j, double * /*data*/)
{
  output_data[n] = (atom->x[i][0] + atom->x[j][0]) * 0.5;
}

/* ---------------------------------------------------------------------- */

void PairTracker::pack_y(int n, int i, int j, double * /*data*/)
{
  output_data[n] = (atom->x[i][1] + atom->x[j][1]) * 0.5;
}

/* ---------------------------------------------------------------------- */

void PairTracker::pack_z(int n, int i, int j, double * /*data*/)
{
  output_data[n] = (atom->x[i][2] + atom->x[j][2]) * 0.5;
}

/* ---------------------------------------------------------------------- */

void PairTracker::pack_rmin(int n, int /*i*/, int /*j*/, double *data)
{
  output_data[n] = data[2];
}

/* ---------------------------------------------------------------------- */

void PairTracker::pack_rave(int n, int /*i*/, int /*j*/, double *data)
{
  output_data[n] = data[1] / (update->ntimestep - data[0]);
}
