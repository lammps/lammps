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

#include "compute_coord_atom.h"

#include "atom.h"
#include "comm.h"
#include "compute_orientorder_atom.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCoordAtom::ComputeCoordAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), typelo(nullptr), typehi(nullptr), cvec(nullptr), carray(nullptr),
    group2(nullptr), id_orientorder(nullptr), normv(nullptr)
{
  if (narg < 5) error->all(FLERR, "Illegal compute coord/atom command");

  jgroup = group->find("all");
  jgroupbit = group->bitmask[jgroup];
  cstyle = NONE;

  if (strcmp(arg[3], "cutoff") == 0) {
    cstyle = CUTOFF;
    double cutoff = utils::numeric(FLERR, arg[4], false, lmp);
    cutsq = cutoff * cutoff;

    int iarg = 5;
    if ((narg > 6) && (strcmp(arg[5], "group") == 0)) {
      group2 = utils::strdup(arg[6]);
      iarg += 2;
      jgroup = group->find(group2);
      if (jgroup == -1) error->all(FLERR, "Compute coord/atom group2 ID does not exist");
      jgroupbit = group->bitmask[jgroup];
    }

    ncol = narg - iarg + 1;
    int ntypes = atom->ntypes;
    typelo = new int[ncol];
    typehi = new int[ncol];

    if (narg == iarg) {
      ncol = 1;
      typelo[0] = 1;
      typehi[0] = ntypes;
    } else {
      ncol = 0;
      while (iarg < narg) {
        utils::bounds(FLERR, arg[iarg], 1, ntypes, typelo[ncol], typehi[ncol], error);
        if (typelo[ncol] > typehi[ncol]) error->all(FLERR, "Illegal compute coord/atom command");
        ncol++;
        iarg++;
      }
    }

  } else if (strcmp(arg[3], "orientorder") == 0) {
    cstyle = ORIENT;
    if (narg != 6) error->all(FLERR, "Illegal compute coord/atom command");

    id_orientorder = utils::strdup(arg[4]);

    int iorientorder = modify->find_compute(id_orientorder);
    if (iorientorder < 0) error->all(FLERR, "Could not find compute coord/atom compute ID");
    if (!utils::strmatch(modify->compute[iorientorder]->style, "^orientorder/atom"))
      error->all(FLERR, "Compute coord/atom compute ID is not orientorder/atom");

    threshold = utils::numeric(FLERR, arg[5], false, lmp);
    if (threshold <= -1.0 || threshold >= 1.0)
      error->all(FLERR, "Compute coord/atom threshold not between -1 and 1");

    ncol = 1;
    typelo = new int[ncol];
    typehi = new int[ncol];
    typelo[0] = 1;
    typehi[0] = atom->ntypes;

  } else
    error->all(FLERR, "Invalid cstyle in compute coord/atom");

  peratom_flag = 1;
  if (ncol == 1)
    size_peratom_cols = 0;
  else
    size_peratom_cols = ncol;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeCoordAtom::~ComputeCoordAtom()
{
  if (copymode) return;

  delete[] group2;
  delete[] typelo;
  delete[] typehi;
  memory->destroy(cvec);
  memory->destroy(carray);
  delete[] id_orientorder;
}

/* ---------------------------------------------------------------------- */

void ComputeCoordAtom::init()
{
  if (cstyle == ORIENT) {
    int iorientorder = modify->find_compute(id_orientorder);
    c_orientorder = dynamic_cast<ComputeOrientOrderAtom *>(modify->compute[iorientorder]);
    cutsq = c_orientorder->cutsq;
    l = c_orientorder->qlcomp;
    //  communicate real and imaginary 2*l+1 components of the normalized vector
    comm_forward = 2 * (2 * l + 1);
    if (!(c_orientorder->qlcompflag))
      error->all(FLERR,
                 "Compute coord/atom requires components option in compute orientorder/atom");
  }

  if (force->pair == nullptr)
    error->all(FLERR, "Compute coord/atom requires a pair style be defined");
  if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR, "Compute coord/atom cutoff is longer than pairwise cutoff");

  // need an occasional full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);
}

/* ---------------------------------------------------------------------- */

void ComputeCoordAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeCoordAtom::compute_peratom()
{
  int i, j, m, ii, jj, inum, jnum, jtype, n;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double *count;

  invoked_peratom = update->ntimestep;

  // grow coordination array if necessary

  if (atom->nmax > nmax) {
    if (ncol == 1) {
      memory->destroy(cvec);
      nmax = atom->nmax;
      memory->create(cvec, nmax, "coord/atom:cvec");
      vector_atom = cvec;
    } else {
      memory->destroy(carray);
      nmax = atom->nmax;
      memory->create(carray, nmax, ncol, "coord/atom:carray");
      array_atom = carray;
    }
  }

  if (cstyle == ORIENT) {
    if (!(c_orientorder->invoked_flag & Compute::INVOKED_PERATOM)) {
      c_orientorder->compute_peratom();
      c_orientorder->invoked_flag |= Compute::INVOKED_PERATOM;
    }
    nqlist = c_orientorder->nqlist;
    normv = c_orientorder->array_atom;
    comm->forward_comm(this);
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // compute coordination number(s) for each atom in group
  // use full neighbor list to count atoms less than cutoff

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;

  if (cstyle == CUTOFF) {

    if (ncol == 1) {

      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        if (mask[i] & groupbit) {
          xtmp = x[i][0];
          ytmp = x[i][1];
          ztmp = x[i][2];
          jlist = firstneigh[i];
          jnum = numneigh[i];

          n = 0;
          for (jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            j &= NEIGHMASK;

            if (mask[j] & jgroupbit) {
              jtype = type[j];
              delx = xtmp - x[j][0];
              dely = ytmp - x[j][1];
              delz = ztmp - x[j][2];
              rsq = delx * delx + dely * dely + delz * delz;
              if (rsq < cutsq && jtype >= typelo[0] && jtype <= typehi[0]) n++;
            }
          }

          cvec[i] = n;
        } else
          cvec[i] = 0.0;
      }

    } else {
      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        count = carray[i];
        for (m = 0; m < ncol; m++) count[m] = 0.0;

        if (mask[i] & groupbit) {
          xtmp = x[i][0];
          ytmp = x[i][1];
          ztmp = x[i][2];
          jlist = firstneigh[i];
          jnum = numneigh[i];

          for (jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            j &= NEIGHMASK;

            if (mask[j] & jgroupbit) {
              jtype = type[j];
              delx = xtmp - x[j][0];
              dely = ytmp - x[j][1];
              delz = ztmp - x[j][2];
              rsq = delx * delx + dely * dely + delz * delz;
              if (rsq < cutsq) {
                for (m = 0; m < ncol; m++)
                  if (jtype >= typelo[m] && jtype <= typehi[m]) count[m] += 1.0;
              }
            }
          }
        }
      }
    }

  } else if (cstyle == ORIENT) {

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      if (mask[i] & groupbit) {
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];

        n = 0;
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx * delx + dely * dely + delz * delz;
          if (rsq < cutsq) {
            double dot_product = 0.0;
            for (m = 0; m < 2 * (2 * l + 1); m++) {
              dot_product += normv[i][nqlist + m] * normv[j][nqlist + m];
            }
            if (dot_product > threshold) n++;
          }
        }
        cvec[i] = n;
      } else
        cvec[i] = 0.0;
    }
  }
}

/* ---------------------------------------------------------------------- */

int ComputeCoordAtom::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                        int * /*pbc*/)
{
  int i, m = 0, j;
  for (i = 0; i < n; ++i) {
    for (j = nqlist; j < nqlist + 2 * (2 * l + 1); ++j) { buf[m++] = normv[list[i]][j]; }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeCoordAtom::unpack_forward_comm(int n, int first, double *buf)
{
  int i, last, m = 0, j;
  last = first + n;
  for (i = first; i < last; ++i) {
    for (j = nqlist; j < nqlist + 2 * (2 * l + 1); ++j) { normv[i][j] = buf[m++]; }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeCoordAtom::memory_usage()
{
  double bytes = (double) ncol * nmax * sizeof(double);
  return bytes;
}
