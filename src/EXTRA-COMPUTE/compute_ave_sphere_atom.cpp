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

/* ----------------------------------------------------------------------
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "compute_ave_sphere_atom.h"

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

ComputeAveSphereAtom::ComputeAveSphereAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), result(nullptr)
{
  if (narg < 3 || narg > 5) error->all(FLERR, "Illegal compute ave/sphere/atom command");

  // process optional args

  cutoff = 0.0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "cutoff") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute ave/sphere/atom command");
      cutoff = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (cutoff <= 0.0) error->all(FLERR, "Illegal compute ave/sphere/atom command");
      iarg += 2;
    } else
      error->all(FLERR, "Illegal compute ave/sphere/atom command");
  }

  peratom_flag = 1;
  size_peratom_cols = 2;
  comm_forward = 3;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeAveSphereAtom::~ComputeAveSphereAtom()
{
  if (copymode) return;

  memory->destroy(result);
}

/* ---------------------------------------------------------------------- */

void ComputeAveSphereAtom::init()
{
  if (!force->pair && cutoff == 0.0)
    error->all(FLERR,
               "Compute ave/sphere/atom requires a cutoff be specified "
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
                 "Compute ave/sphere/atom cutoff exceeds ghost atom range - "
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

void ComputeAveSphereAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeAveSphereAtom::compute_peratom()
{
  int i, j, ii, jj, inum, jnum;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int count;
  double p[3], vcom[3], vnet[3];

  invoked_peratom = update->ntimestep;

  // grow result array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(result);
    nmax = atom->nmax;
    memory->create(result, nmax, 2, "ave/sphere/atom:result");
    array_atom = result;
  }

  // need velocities of ghost atoms

  comm->forward_comm(this);

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // compute properties for each atom in group
  // use full neighbor list to count atoms less than cutoff

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  double massone_i, massone_j, totalmass;

  double adof = domain->dimension;
  double mvv2e = force->mvv2e;
  double mv2d = force->mv2d;
  double boltz = force->boltz;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    if (mask[i] & groupbit) {
      if (rmass)
        massone_i = rmass[i];
      else
        massone_i = mass[type[i]];

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      // i atom contribution

      count = 1;
      totalmass = massone_i;
      p[0] = v[i][0] * massone_i;
      p[1] = v[i][1] * massone_i;
      p[2] = v[i][2] * massone_i;

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        if (rmass)
          massone_j = rmass[j];
        else
          massone_j = mass[type[j]];

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx * delx + dely * dely + delz * delz;
        if (rsq < cutsq) {
          count++;
          totalmass += massone_j;
          p[0] += v[j][0] * massone_j;
          p[1] += v[j][1] * massone_j;
          p[2] += v[j][2] * massone_j;
        }
      }

      vcom[0] = p[0] / totalmass;
      vcom[1] = p[1] / totalmass;
      vcom[2] = p[2] / totalmass;

      // i atom contribution

      vnet[0] = v[i][0] - vcom[0];
      vnet[1] = v[i][1] - vcom[1];
      vnet[2] = v[i][2] - vcom[2];
      double ke_sum = massone_i * (vnet[0] * vnet[0] + vnet[1] * vnet[1] + vnet[2] * vnet[2]);

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        if (rmass)
          massone_j = rmass[j];
        else
          massone_j = mass[type[j]];

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx * delx + dely * dely + delz * delz;
        if (rsq < cutsq) {
          vnet[0] = v[j][0] - vcom[0];
          vnet[1] = v[j][1] - vcom[1];
          vnet[2] = v[j][2] - vcom[2];
          ke_sum += massone_j * (vnet[0] * vnet[0] + vnet[1] * vnet[1] + vnet[2] * vnet[2]);
        }
      }
      double density = mv2d * totalmass / volume;
      double temp = mvv2e * ke_sum / (adof * count * boltz);
      result[i][0] = density;
      result[i][1] = temp;
    }
  }
}

/* ---------------------------------------------------------------------- */

int ComputeAveSphereAtom::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                            int * /*pbc*/)
{
  double **v = atom->v;

  int i, m = 0;
  for (i = 0; i < n; ++i) {
    buf[m++] = v[list[i]][0];
    buf[m++] = v[list[i]][1];
    buf[m++] = v[list[i]][2];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeAveSphereAtom::unpack_forward_comm(int n, int first, double *buf)
{
  double **v = atom->v;

  int i, last, m = 0;
  last = first + n;
  for (i = first; i < last; ++i) {
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeAveSphereAtom::memory_usage()
{
  double bytes = (double) 2 * nmax * sizeof(double);
  return bytes;
}
