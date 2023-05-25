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

/* ----------------------------------------------------------------------
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "compute_local_composition_atom.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ComputeLocalCompositionAtom::ComputeLocalCompositionAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), result(nullptr)
{
  if (narg < 3 || narg > 5) error->all(FLERR, "Illegal compute local_composition/atom command");

  // get nelements and ntypes

  int ntypes = atom->ntypes;

  // memory->create(map, ntypes + 1, "compute_sna_grid:map");
  // nelements = utils::inumeric(FLERR, 0, false, lmp); // !!!!!!!! what is 2nd arg in inumeric?
  // for (int i = 0; i < ntypes; i++) {
  //   int jelem = utils::inumeric(FLERR, i, false, lmp);
  //   if (jelem < 0 || jelem >= nelements) error->all(FLERR, "Illegal compute {} command", style);
  //   map[i + 1] = jelem;
  //   printf("mapp[x] %d jelem %d \n", map[i + 1], jelem);
  // }

  // printf("ntypes %d \n", ntypes);
  // printf("elements %d \n", nelements);

  // process optional args

  cutoff = 0.0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "cutoff") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute local_composition/atom command");
      cutoff = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (cutoff <= 0.0) error->all(FLERR, "Illegal compute local_composition/atom command");
      iarg += 2;
    } else
      error->all(FLERR, "Illegal compute local_composition/atom command");
  }

  peratom_flag = 1;

  size_peratom_cols = 1 + ntypes;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeLocalCompositionAtom::~ComputeLocalCompositionAtom()
{
  if (copymode) return;

  memory->destroy(result);
  // memory->destroy(map);
}

/* ---------------------------------------------------------------------- */

void ComputeLocalCompositionAtom::init()
{
  if (!force->pair && cutoff == 0.0)
    error->all(FLERR,
               "Compute local_composition/atom requires a cutoff be specified "
               "or a pair style be defined");

  double skin = neighbor->skin;
  if (cutoff != 0.0) {
    double cutghost;    // as computed by Neighbor and Comm
    if (force->pair)
      cutghost = MAX(force->pair->cutforce + skin, comm->cutghostuser);
    else
      cutghost = comm->cutghostuser;

    if (cutoff > cutghost)
      error->all(FLERR,
                 "Compute local_composition/atom cutoff exceeds ghost atom range - "
                 "use comm_modify cutoff command");
  }

  int cutflag = 1;
  if (force->pair) {
    if (cutoff == 0.0) { cutoff = force->pair->cutforce; }
    if (cutoff <= force->pair->cutforce + skin) cutflag = 0;
  }

  cutsq = cutoff * cutoff;
  if (domain->dimension == 3)
    volume = 4.0 / 3.0 * MY_PI * cutsq * cutoff;
  else
    volume = MY_PI * cutsq;

  // need an occasional full neighbor list

  auto req = neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);
  if (cutflag) req->set_cutoff(cutoff);
}

/* ---------------------------------------------------------------------- */

void ComputeLocalCompositionAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeLocalCompositionAtom::compute_peratom()
{
  int i, j, ii, jj, inum, jnum;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int count;

  invoked_peratom = update->ntimestep;

  // grow result array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(result);
    nmax = atom->nmax;
    memory->create(result, nmax, size_peratom_cols, "local_composition/atom:result");
    array_atom = result;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // compute properties for each atom in group
  // use full neighbor list to count atoms less than cutoff

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;

  int typeone_i, typeone_j;

  // TODO continue to implement map (nelements instead of ntypes)
  
  int ntypes = atom->ntypes;
  
  double lcomp[ntypes];

  // get per-atom local compositions

  for (ii = 0; ii < inum; ii++) {

    for (int i = 0; i < ntypes; i++) {
      lcomp[i] = 0;
    }

    i = ilist[ii];

    if (mask[i] & groupbit) {

      typeone_i = type[i];

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      // i atom contribution

      count = 1;

      typeone_i = type[i];
      lcomp[typeone_i - 1]++;

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        typeone_j = type[j];

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx * delx + dely * dely + delz * delz;
        if (rsq < cutsq) {
          count++;
          lcomp[typeone_j-1]++;
        }
      }

      // total count of atoms found in sampled radius range

      result[i][0] = count;

      // local composition fractions per element
      for (int n = 0; n < ntypes; n++) {
        result[i][n + 1] = lcomp[n] / count;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeLocalCompositionAtom::memory_usage()
{
  double bytes = (double) 2 * nmax * sizeof(double);
  return bytes;
}
