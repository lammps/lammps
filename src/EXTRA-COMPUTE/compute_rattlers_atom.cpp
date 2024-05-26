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

/* ----------------------------------------------------------------------
   Contributing authors: Joel Clemmer (SNL), Ishan Srivastava (LBNL)
------------------------------------------------------------------------- */

#include "compute_rattlers_atom.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

enum { TYPE, RADIUS };

/* ---------------------------------------------------------------------- */

ComputeRattlersAtom::ComputeRattlersAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), ncontacts(nullptr), rattler(nullptr)
{
  if (narg != 6) error->all(FLERR, "Illegal compute rattlers/atom command");

  if (strcmp(arg[3], "type") == 0)
    cutstyle = TYPE;
  else if (strcmp(arg[3], "radius") == 0)
    cutstyle = RADIUS;
  else
    error->all(FLERR, "Illegal compute rattlers/atom command");

  if (cutstyle == RADIUS && !atom->radius_flag)
    error->all(FLERR, "Compute rattlers/atom radius style requires atom attribute radius");

  ncontacts_rattler = utils::inumeric(FLERR, arg[4], false, lmp);
  max_tries = utils::inumeric(FLERR, arg[5], false, lmp);

  nmax = 0;
  invoked_peratom = -1;

  scalar_flag = 1;
  extscalar = 1;
  peratom_flag = 1;
  size_peratom_cols = 0;
  comm_forward = 1;
  comm_reverse = 1;
}

/* ---------------------------------------------------------------------- */

ComputeRattlersAtom::~ComputeRattlersAtom()
{
  memory->destroy(ncontacts);
  memory->destroy(rattler);
}

/* ---------------------------------------------------------------------- */

void ComputeRattlersAtom::init()
{
  if (force->pair == nullptr) error->all(FLERR, "No pair style is defined for compute rattlers");

  // Cannot calculate distance from radii for JKR/DMT
  if (force->pair->beyond_contact)
    error->all(FLERR, "Compute rattlers does not currently support pair styles that extend beyond contact");

  // need an occasional half neighbor list
  // set size to same value as request made by force->pair
  // this should enable it to always be a copy list (e.g. for granular pstyle)

  auto pairrequest = neighbor->find_request(force->pair);
  if (pairrequest && pairrequest->get_size())
    neighbor->add_request(this, NeighConst::REQ_SIZE | NeighConst::REQ_OCCASIONAL);
  else
    neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);
}

/* ---------------------------------------------------------------------- */

void ComputeRattlersAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeRattlersAtom::compute_peratom()
{
  if (invoked_peratom == update->ntimestep) return;
  invoked_peratom = update->ntimestep;

  int i, j, ii, jj, inum, jnum, itype, jtype, tmp_flag;
  tagint itag, jtag;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq, radsum;

  if (nmax < atom->nmax) {
    nmax = atom->nmax;
    memory->destroy(ncontacts);
    memory->destroy(rattler);
    memory->create(ncontacts, nmax, "rattlers:ncontacts");
    memory->create(rattler, nmax, "rattlers:rattler");
    vector_atom = rattler;
  }

  for (i = 0; i < nmax; i++) rattler[i] = 0;

  int *ilist, *jlist, *numneigh, **firstneigh;

  double **x = atom->x;
  double *radius = atom->radius;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  // invoke half neighbor list (will copy or build if necessary)
  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double **cutsq = force->pair->cutsq;

  int change_flag = 1;
  int ntry = 0;
  while (ntry < max_tries) {
    change_flag = 0;

    for (i = 0; i < nmax; i++) ncontacts[i] = 0;

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      if (!(mask[i] & groupbit)) continue;
      if (rattler[i] == 1) continue;

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      itag = tag[i];
      itype = type[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        if (!(mask[j] & groupbit)) continue;
        if (rattler[j] == 1) continue;

        // itag = jtag is possible for long cutoffs that include images of self

        if (newton_pair == 0 && j >= nlocal) {
          jtag = tag[j];
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

        jtype = type[j];

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx * delx + dely * dely + delz * delz;

        if (cutstyle == TYPE) {
          if (rsq >= cutsq[itype][jtype]) continue;
        } else {
          radsum = radius[i] + radius[j];
          if (rsq >= radsum * radsum) continue;
        }
        ncontacts[i] += 1;
        if (newton_pair || j < nlocal)
          ncontacts[j] += 1;
      }
    }

    // add contributions from ghosts
    if (force->newton_pair) comm->reverse_comm(this);

    // Set flags for rattlers
    for (i = 0; i < atom->nlocal; i++) {
      if (ncontacts[i] < ncontacts_rattler && rattler[i] == 0) {
        rattler[i] = 1;
        change_flag = 1;
      }
    }

    comm->forward_comm(this);

    MPI_Allreduce(&change_flag, &tmp_flag, 1, MPI_INT, MPI_MAX, world);
    change_flag = tmp_flag;
    if (change_flag == 0) break;

    ntry += 1;
  }

  if (change_flag == 1)
    error->warning(FLERR, "Rattler calculation failed to converge within max tries");
}

/* ---------------------------------------------------------------------- */

double ComputeRattlersAtom::compute_scalar()
{
  if (invoked_peratom != update->ntimestep)
    compute_peratom();

  invoked_scalar = update->ntimestep;

  double total_rattlers = 0;
  for (int i = 0; i < atom->nlocal; i++) {
    if (rattler[i] == 1) {
      total_rattlers += 1;
    }
  }

  //Total across processors
  MPI_Allreduce(&total_rattlers, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);
  return scalar;
}

/* ---------------------------------------------------------------------- */

int ComputeRattlersAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = ubuf(ncontacts[i]).d;
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeRattlersAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    ncontacts[j] += (int) ubuf(buf[m++]).i;
  }
}

/* ---------------------------------------------------------------------- */

int ComputeRattlersAtom::pack_forward_comm(int n, int *list, double *buf,
                                          int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = rattler[j];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeRattlersAtom::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    rattler[i] = buf[m++];
  }
}
