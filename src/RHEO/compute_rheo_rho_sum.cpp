/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_rheo_rho_sum.h"

#include "atom.h"
#include "comm.h"
#include "compute_rheo_kernel.h"
#include "error.h"
#include "fix_rheo.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeRHEORhoSum::ComputeRHEORhoSum(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), fix_rheo(nullptr), compute_kernel(nullptr)
{
  if (narg != 3) error->all(FLERR,"Illegal compute RHEO/rho command");

  comm_forward = 1;
  comm_reverse = 1;
}

/* ---------------------------------------------------------------------- */

ComputeRHEORhoSum::~ComputeRHEORhoSum() {}

/* ---------------------------------------------------------------------- */

void ComputeRHEORhoSum::init()
{
  compute_kernel = fix_rheo->compute_kernel;
  cut = fix_rheo->cut;
  cutsq = cut * cut;

  // need an occasional half neighbor list
  neighbor->add_request(this, NeighConst::REQ_DEFAULT);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEORhoSum::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */


void ComputeRHEORhoSum::compute_peratom()
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double rsq, w;

  int nlocal = atom->nlocal;

  double **x = atom->x;
  double *rho = atom->rho;
  int *type = atom->type;
  double *mass = atom->mass;
  int newton = force->newton;

  double jmass;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  int nall = nlocal + atom->nghost;

  // initialize arrays, local with quintic self-contribution, ghosts are zeroed
  for (i = 0; i < nlocal; i++) {
    w = compute_kernel->calc_w_quintic(i, i, 0.0, 0.0, 0.0, 0.0);
    rho[i] = w * mass[type[i]];
  }

  for (i = nlocal; i < nall; i++) rho[i] = 0.0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq < cutsq) {
        w = compute_kernel->calc_w(i, j, delx, dely, delz, sqrt(rsq));
        rho[i] += w * mass[type[i]];
        if (newton || j < nlocal) {
          rho[j] += w * mass[type[j]];
        }
      }
    }
  }

  if (newton) comm->reverse_comm(this);
  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

int ComputeRHEORhoSum::pack_forward_comm(int n, int *list, double *buf,
                                        int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, k, m;
  double *rho = atom->rho;
  m = 0;

  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = rho[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */
void ComputeRHEORhoSum::unpack_forward_comm(int n, int first, double *buf)
{
  int i, k, m, last;
  double *rho = atom->rho;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    rho[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int ComputeRHEORhoSum::pack_reverse_comm(int n, int first, double *buf)
{
  int i, k, m, last;
  double *rho = atom->rho;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = rho[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEORhoSum::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, k, j, m;
  double *rho = atom->rho;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    rho[j] += buf[m++];
  }
}
