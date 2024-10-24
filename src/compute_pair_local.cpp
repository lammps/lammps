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

#include "compute_pair_local.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

static constexpr int DELTA = 10000;

enum { DIST, ENG, FORCE, FX, FY, FZ, PN, DX, DY, DZ };
enum { TYPE, RADIUS };

/* ---------------------------------------------------------------------- */

ComputePairLocal::ComputePairLocal(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), pstyle(nullptr), pindex(nullptr), vlocal(nullptr), alocal(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "compute pair/local", error);

  local_flag = 1;
  nvalues = narg - 3;
  pstyle = new int[nvalues];
  pindex = new int[nvalues];

  nvalues = 0;
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "dist") == 0)
      pstyle[nvalues++] = DIST;
    else if (strcmp(arg[iarg], "eng") == 0)
      pstyle[nvalues++] = ENG;
    else if (strcmp(arg[iarg], "force") == 0)
      pstyle[nvalues++] = FORCE;
    else if (strcmp(arg[iarg], "fx") == 0)
      pstyle[nvalues++] = FX;
    else if (strcmp(arg[iarg], "fy") == 0)
      pstyle[nvalues++] = FY;
    else if (strcmp(arg[iarg], "fz") == 0)
      pstyle[nvalues++] = FZ;
    else if (strcmp(arg[iarg], "dx") == 0)
      pstyle[nvalues++] = DX;
    else if (strcmp(arg[iarg], "dy") == 0)
      pstyle[nvalues++] = DY;
    else if (strcmp(arg[iarg], "dz") == 0)
      pstyle[nvalues++] = DZ;
    else if (utils::strmatch(arg[iarg], "^p\\d+$")) {    // p1, p2, p3, ... pN
      int n = std::stoi(&arg[iarg][1]);
      if (n <= 0) error->all(FLERR, "Invalid keyword {} in compute pair/local command", arg[iarg]);
      pstyle[nvalues] = PN;
      pindex[nvalues++] = n - 1;
    } else
      break;
    iarg++;
  }

  // optional args

  cutstyle = TYPE;

  // group2-ID defaults to group-ID of compute if group2 option not used
  igroup2 = igroup;
  groupbit2 = groupbit;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "cutoff") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "compute pair/local cutoff", error);
      if (strcmp(arg[iarg + 1], "type") == 0)
        cutstyle = TYPE;
      else if (strcmp(arg[iarg + 1], "radius") == 0)
        cutstyle = RADIUS;
      else
        error->all(FLERR, "Unknown compute pair/local cutoff keyword: {}", arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "group2") == 0) {
      if (iarg + 1 > narg) error->all(FLERR, "Illegal compute property/local command");
      igroup2 = group->find(arg[iarg + 1]);
      if (igroup2 == -1) error->all(FLERR,"Could not find compute group2-ID {}", arg[iarg + 1] );
      groupbit2 = group->bitmask[igroup2];
      iarg += 2;
    } else
      error->all(FLERR, "Unknown compute pair/local keyword: {}", arg[iarg]);
  }

  // error check

  if (cutstyle == RADIUS && !atom->radius_flag)
    error->all(FLERR, "Compute pair/local with ID {} requires atom attribute radius", id);

  // set singleflag if need to call pair->single()

  singleflag = 0;
  for (int i = 0; i < nvalues; i++)
    if (pstyle[i] != DIST && pstyle[i] != DX && pstyle[i] != DY && pstyle[i] != DZ) singleflag = 1;

  if (nvalues == 1)
    size_local_cols = 0;
  else
    size_local_cols = nvalues;

  nmax = 0;
  vlocal = nullptr;
  alocal = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputePairLocal::~ComputePairLocal()
{
  memory->destroy(vlocal);
  memory->destroy(alocal);
  delete[] pstyle;
  delete[] pindex;
}

/* ---------------------------------------------------------------------- */

void ComputePairLocal::init()
{
  if (singleflag && force->pair == nullptr)
    error->all(FLERR, "No pair style is defined for compute pair/local");
  if (singleflag && force->pair->single_enable == 0)
    error->all(FLERR, "Pair style does not support compute pair/local");

  for (int i = 0; i < nvalues; i++)
    if (pstyle[i] == PN && pindex[i] >= force->pair->single_extra)
      error->all(FLERR, "Pair style does not have extra field requested by compute pair/local");

  // need an occasional half neighbor list
  // set size to same value as request made by force->pair
  // this should enable it to always be a copy list (e.g. for granular pstyle)

  int neighflags = NeighConst::REQ_OCCASIONAL;
  auto pairrequest = neighbor->find_request(force->pair);
  if (pairrequest && pairrequest->get_size()) neighflags |= NeighConst::REQ_SIZE;
  neighbor->add_request(this, neighflags);
}

/* ---------------------------------------------------------------------- */

void ComputePairLocal::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputePairLocal::compute_local()
{
  invoked_local = update->ntimestep;

  // count local entries and compute pair info

  ncount = compute_pairs(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
  compute_pairs(1);
}

/* ----------------------------------------------------------------------
   count pairs and compute pair info on this proc
   only count pair once if newton_pair is off
   both atom I,J must be in group
   if flag is set, compute requested info about pair
------------------------------------------------------------------------- */

int ComputePairLocal::compute_pairs(int flag)
{
  int i, j, m, n, ii, jj, inum, jnum, itype, jtype;
  tagint itag, jtag;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq, radsum, eng, fpair, factor_coul, factor_lj;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double *ptr;

  double **x = atom->x;
  double *radius = atom->radius;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  // invoke half neighbor list (will copy or build if necessary)

  if (flag == 0) neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  // skip if I or J are not in group
  // for newton = 0 and J = ghost atom,
  //   need to ensure I,J pair is only output by one proc
  //   use same itag,jtag logic as in Neighbor::neigh_half_nsq()
  // for flag = 0, just count pair interactions within force cutoff
  // for flag = 1, calculate requested output fields

  Pair *pair = force->pair;
  double **cutsq = force->pair->cutsq;

  m = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if ( !(mask[i] & groupbit) && !(mask[i] & groupbit2) ) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itag = tag[i];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;
      jtag = tag[j];

      if ( (mask[i] & groupbit) && !(mask[j] & groupbit2)) continue;
      if ( (mask[i] & groupbit2) && !(mask[j] & groupbit)) continue;

      // itag = jtag is possible for long cutoffs that include images of self

      if (newton_pair == 0 && j >= nlocal) {
        if (itag > jtag) {
          if ((itag + jtag) % 2 == 0) continue;
        } else if (itag < jtag) {
          if ((itag + jtag) % 2 == 1) continue;
        } else {
          if (x[j][2] < ztmp) continue;
          if (x[j][2] == ztmp) {
            if (x[j][1] < ytmp) continue;
            if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
          }
        }
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (cutstyle == TYPE) {
        if (rsq >= cutsq[itype][jtype]) continue;
      } else {
        radsum = radius[i] + radius[j];
        if (rsq >= radsum * radsum) continue;
      }

      if (flag) {
        if (singleflag)
          eng = pair->single(i, j, itype, jtype, rsq, factor_coul, factor_lj, fpair);
        else
          eng = fpair = 0.0;

        if (nvalues == 1)
          ptr = &vlocal[m];
        else
          ptr = alocal[m];

        // to make sure dx, dy and dz are always from the lower to the higher id
        double directionCorrection = itag > jtag ? -1.0 : 1.0;

        for (n = 0; n < nvalues; n++) {
          switch (pstyle[n]) {
            case DIST:
              ptr[n] = sqrt(rsq);
              break;
            case DX:
              ptr[n] = delx * directionCorrection;
              break;
            case DY:
              ptr[n] = dely * directionCorrection;
              break;
            case DZ:
              ptr[n] = delz * directionCorrection;
              break;
            case ENG:
              ptr[n] = eng;
              break;
            case FORCE:
              ptr[n] = sqrt(rsq) * fpair;
              break;
            case FX:
              ptr[n] = delx * fpair;
              break;
            case FY:
              ptr[n] = dely * fpair;
              break;
            case FZ:
              ptr[n] = delz * fpair;
              break;
            case PN:
              ptr[n] = pair->svector[pindex[n]];
              break;
          }
        }
      }

      m++;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputePairLocal::reallocate(int n)
{
  // grow vector_local or array_local

  while (nmax < n) nmax += DELTA;

  if (nvalues == 1) {
    memory->destroy(vlocal);
    memory->create(vlocal, nmax, "pair/local:vector_local");
    vector_local = vlocal;
  } else {
    memory->destroy(alocal);
    memory->create(alocal, nmax, nvalues, "pair/local:array_local");
    array_local = alocal;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputePairLocal::memory_usage()
{
  double bytes = (double) nmax * nvalues * sizeof(double);
  return bytes;
}
