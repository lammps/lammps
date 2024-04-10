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
   Contributing author: Megan McCarthy (SNL)
------------------------------------------------------------------------- */

#include "compute_composition_atom.h"

#include "atom.h"
#include "comm.h"
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

ComputeCompositionAtom::ComputeCompositionAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), list(nullptr), result(nullptr)
{
  if (narg < 3 || narg > 5) error->all(FLERR, "Illegal compute composition/atom command");

  cutsq = cutoff = 0.0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "cutoff") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute composition/atom command");
      cutoff = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (cutoff <= 0.0) error->all(FLERR, "Illegal compute composition/atom command");
      iarg += 2;
    } else
      error->all(FLERR, "Illegal compute composition/atom command");
  }

  peratom_flag = 1;

  ntypes = atom->ntypes;
  size_peratom_cols = 1 + ntypes;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeCompositionAtom::~ComputeCompositionAtom()
{
  if (copymode) return;

  memory->destroy(result);
}

/* ---------------------------------------------------------------------- */

void ComputeCompositionAtom::init()
{
  if (!force->pair && cutoff == 0.0)
    error->all(FLERR,
               "Compute composition/atom requires a cutoff be specified "
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
                 "Compute composition/atom cutoff exceeds ghost atom range - "
                 "use comm_modify cutoff command");
  }

  int cutflag = 1;
  if (force->pair) {
    if (cutoff == 0.0) { cutoff = force->pair->cutforce; }
    if (cutoff <= force->pair->cutforce + skin) cutflag = 0;
  }

  cutsq = cutoff * cutoff;

  if ((neighbor->style == Neighbor::MULTI) || (neighbor->style == Neighbor::MULTI_OLD))
    error->all(FLERR, "Compute composition/atom requires neighbor style 'bin' or 'nsq'");

  // need an occasional full neighbor list

  auto req = neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);
  if (cutflag) req->set_cutoff(cutoff);
}

/* ---------------------------------------------------------------------- */

void ComputeCompositionAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeCompositionAtom::compute_peratom()
{
  int i, j, ii, jj, inum, jnum;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int count, itype, jtype;

  invoked_peratom = update->ntimestep;

  // grow result array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(result);
    nmax = atom->nmax;
    memory->create(result, nmax, size_peratom_cols, "composition/atom:result");
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

  // get per-atom local compositions

  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii];

    if (mask[i] & groupbit) {

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      // i atom contribution

      count = 1;

      itype = type[i];
      result[i][itype]++;

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        jtype = type[j];

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx * delx + dely * dely + delz * delz;
        if (rsq < cutsq) {
          count++;
          result[i][jtype]++;
        }
      }

      // total count of atoms found in sampled radius range

      result[i][0] = count;

      // local comp fractions per element

      double lfac = 1.0 / count;
      for (int n = 1; n < size_peratom_cols; n++) result[i][n + 1] *= lfac;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeCompositionAtom::memory_usage()
{
  double bytes = (double) 2 * nmax * sizeof(double);
  return bytes;
}
