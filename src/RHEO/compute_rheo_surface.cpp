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
   Contributing authors:
   Joel Clemmer (SNL), Thomas O'Connor (CMU), Eric Palermo (CMU)
----------------------------------------------------------------------- */

#include "compute_rheo_surface.h"

#include "atom.h"
#include "comm.h"
#include "compute_rheo_kernel.h"
#include "compute_rheo_solids.h"
#include "domain.h"
#include "error.h"
#include "fix_rheo.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathExtra;

#define EPSILON 1e-10;

/* ---------------------------------------------------------------------- */

ComputeRHEOSurface::ComputeRHEOSurface(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), fix_rheo(nullptr), list(nullptr), compute_kernel(nullptr), compute_solids(nullptr),
  B(nullptr), gradC(nullptr), nsurface(nullptr), divr(nullptr), rsurface(nullptr)
{
  if (narg != 3) error->all(FLERR,"Illegal fix RHEO/SURFACE command");

  int dim = domain->dimension;
  comm_forward = 2;
  comm_reverse = dim * dim + 1;
}

/* ---------------------------------------------------------------------- */

ComputeRHEOSurface::~ComputeRHEOSurface()
{
  // Remove custom property if it exists
  int tmp1, tmp2, index;
  index = atom->find_custom("rheo_divr", tmp1, tmp2);
  if (index != -1) atom->remove_custom(index, 1, 0);

  index = atom->find_custom("rheo_rsurface", tmp1, tmp2);
  if (index != -1) atom->remove_custom(index, 1, 0);

  index = atom->find_custom("rheo_nsurface", tmp1, tmp2);
  if (index != -1) atom->remove_custom(index, 1, 3);

  memory->destroy(B);
  memory->destroy(gradC);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOSurface::init()
{
  compute_kernel = fix_rheo->compute_kernel;
  compute_solids = fix_rheo->compute_solids;
  cut = fix_rheo->cut;
  rho0 = fix_rheo->rho0;
  threshold_style = fix_rheo->surface_style;
  threshold_divr = fix_rheo->divrsurface;
  threshold_z = fix_rheo->zminsurface;

  cutsq = cut * cut;

  // Create rsurface, divr, nsurface arrays if they don't already exist
  // Create a custom atom property so it works with compute property/atom
  // Do not create grow callback as there's no reason to copy/exchange data
  // Manually grow if nmax_old exceeded
  // For B and gradC, create a local array since they are unlikely to be printed

  int tmp1, tmp2;
  int index = atom->find_custom("rheo_divr", tmp1, tmp2);
  if (index == -1)  index = atom->add_custom("rheo_divr", 1, 0);
  divr = atom->dvector[index];

  index = atom->find_custom("rheo_rsurface", tmp1, tmp2);
  if (index == -1)  index = atom->add_custom("rheo_rsurface", 1, 0);
  rsurface = atom->dvector[index];

  index = atom->find_custom("rheo_nsurface", tmp1, tmp2);
  if (index == -1)  index = atom->add_custom("rheo_nsurface", 1, 3);
  nsurface = atom->darray[index];

  nmax_old = atom->nmax;
  memory->create(B, nmax_old, dim * dim, "rheo/surface:B");
  memory->create(gradC, nmax_old, dim * dim, "rheo/surface:gradC");

  // need an occasional half neighbor list
  neighbor->add_request(this, NeighConst::REQ_HALF);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOSurface::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOSurface::compute_peratom()
{
  int i, j, ii, jj, inum, jnum, a, b, itype, jtype, fluidi, fluidj;
  double xtmp, ytmp, ztmp, rsq, Voli, Volj, rhoi, rhoj;
  double *dWij, *dWji;
  double dx[3];
  int *ilist, *jlist, *numneigh, **firstneigh;

  int nlocal = atom->nlocal;

  double **x = atom->x;
  int *status = atom->status;
  int newton = force->newton;
  int dim = domain->dimension;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rho = atom->rho;
  int *coordination = compute_kernel->coordination;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  int nmax = atom->nmax;
  if (nmax_old <= nmax) {
    memory->grow(divr, nmax, "atom:rheo_divr");
    memory->grow(rsurface, nmax, "atom:rheo_rsurface");
    memory->grow(nsurface, nmax, 3, "atom:rheo_nsurface");

    memory->grow(B, nmax, dim * dim, "rheo/surface:B");
    memory->grow(gradC, nmax, dim * dim, "rheo/surface:gradC");

    nmax_old = atom->nmax;
  }

  int nall = nlocal + atom->nghost;
  for (i = 0; i < nall; i++) {
    for (a = 0; a < dim; a++) {
      for (b = 0; b < dim; b++) {
        B[i][a * dim + b] = 0.0;
        gradC[i][a * dim + b] = 0.0;
      }
      nsurface[i][a] = 0.0;
    }
    divr[i] = 0.0;

    // Remove surface settings
    status[i] &= FixRHEO::surfacemask;
  }

  // loop over neighbors to calculate the average orientation of neighbors
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    itype = type[i];
    fluidi = status[i] & FixRHEO::STATUS_FLUID;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      dx[0] = xtmp - x[j][0];
      dx[1] = ytmp - x[j][1];
      dx[2] = ztmp - x[j][2];

      rsq = lensq(dx);
      if (rsq < cutsq) {
        jtype = type[j];
        fluidj = status[j] & FixRHEO::STATUS_FLUID;

        rhoi = rho[i];
        rhoj = rho[j];

        // Add corrections for walls
        if (fluidi && (!fluidj)) {
          rhoj = compute_solids->correct_rho(j, i);
        } else if ((!fluidi) && fluidj) {
          rhoi = compute_solids->correct_rho(i, j);
        } else if ((!fluidi) && (!fluidj)) {
          rhoi = rho0;
          rhoj = rho0;
        }

        Voli = mass[itype] / rhoi;
        Volj = mass[jtype] / rhoj;

        wp = compute_kernel->calc_dw_quintic(i, j, dx[0], dx[1], dx[2], sqrt(rsq),dWij, dWji);

        for (a = 0; a < dim; a++){
          divr[i] -= dWij[a] * dx[a] * Volj;
          gradC[i][a] += dWij[a] * Volj;
        }

        if (j < nlocal || newton) {
          for (a = 0; a < dim; a++){
            divr[j] += dWji[a] * dx[a] * Voli;
            gradC[j][a] += dWji[a] * Voli;
          }
        }
      }
    }
  }

  // reverse gradC and divr, forward divr
  comm_stage = 0;
  comm_reverse = dim * dim + 1;
  comm_forward = 1;
  if (newton) comm->reverse_comm(this);
  comm->forward_comm(this);

  // calculate nsurface for local atoms
  // Note, this isn't forwarded to ghosts
  double maggC;
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      maggC = 0.0;
      for (a = 0;a < dim; a++)
        maggC += gradC[i][a] * gradC[i][a];
      maggC = sqrt(maggC) + EPSILON;
      maggC = 1.0 / maggC;
      for (a = 0; a < dim; a++)
        nsurface[i][a] = -gradC[i][a] * maggC;
    }
  }

  // Find the free-surface
  if (threshold_style == FixRHEO::DIVR) {
    for (i = 0; i < nall; i++) {
      if (mask[i] & groupbit) {
        status[i] |= FixRHEO::STATUS_BULK;
        rsurface[i] = cut;
        if (divr[i] < threshold_divr) {
          status[i] |= FixRHEO::STATUS_SURFACE;
          rsurface[i] = 0.0;
          if (coordination[i] < threshold_z)
            status[i] |= FixRHEO::STATUS_SPLASH;
        }
      }
    }
  } else {
    for (i = 0; i < nall; i++) {
      if (mask[i] & groupbit) {
        status[i] |= FixRHEO::STATUS_BULK;
        rsurface[i] = cut;
        if (coordination[i] < divR_limit) {
          status[i] |= FixRHEO::STATUS_SURFACE;
          rsurface[i] = 0.0;
          if (coordination[i] < threshold_z)
            status[i] |= FixRHEO::STATUS_SPLASH;
        }
      }
    }
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      dx[0] = xtmp - x[j][0];
      dx[1] = ytmp - x[j][1];
      dx[2] = ztmp - x[j][2];
      rsq = lensq(dx);
      if (rsq < cutsq) {
        if ((status[i] & FixRHEO::STATUS_BULK) && (status[j] & FixRHEO::STATUS_SURFACE)) {
          status[i] &= FixRHEO::surfacemask;
          status[i] |= FixRHEO::STATUS_LAYER;
        }

        if (status[j] & FixRHEO::STATUS_SURFACE) rsurface[i] = MIN(rsurface[i], sqrt(rsq));


        if (j < nlocal || newton) {
          if ((status[j] & FixRHEO::STATUS_BULK) && (status[i] & FixRHEO::STATUS_SURFACE)) {
            status[j] &= FixRHEO::surfacemask;
            status[j] |= FixRHEO::STATUS_LAYER;
          }

          if (status[i] & FixRHEO::STATUS_SURFACE) rsurface[j] = MIN(rsurface[j], sqrt(rsq));
        }
      }
    }
  }

  // forward/reverse status and rsurface
  comm_stage = 1;
  comm_reverse = 2;
  comm_forward = 2;
  if (newton) comm->reverse_comm(this);
  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

int ComputeRHEOSurface::pack_reverse_comm(int n, int first, double *buf)
{
  int i,a,b,k,m,last;
  int dim = domain->dimension;
  int *status = atom->status;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (comm_stage == 0) {
      buf[m++] = divr[i];
      for (a = 0; a < dim; a ++ )
        for (b = 0; b < dim; b ++)
          buf[m++] = gradC[i][a * dim + b];
    } else if (comm_stage == 1) {
      buf[m++] = (double) status[i];
      buf[m++] = rsurface[i];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOSurface::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,a,b,k,j,m;
  int dim = domain->dimension;
  int *status = atom->status;
  int temp;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (comm_stage == 0) {
      divr[j] += buf[m++];
      for (a = 0; a < dim; a ++ )
        for (b = 0; b < dim; b ++)
          gradC[j][a * dim + b] += buf[m++];
    } else if (comm_stage == 1) {

      temp = (int) buf[m++];
      if ((status[j] & FixRHEO::STATUS_BULK) && (temp & FixRHEO::STATUS_LAYER))
        status[j] = temp;

      rsurface[j] = MIN(rsurface[j], buf[m++]);
    }
  }
}


/* ---------------------------------------------------------------------- */

int ComputeRHEOSurface::pack_forward_comm(int n, int *list, double *buf,
                                        int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,a,b,k,m;
  int *status = atom->status;
  m = 0;

  for (i = 0; i < n; i++) {
    j = list[i];
    if (comm_stage == 0) {
      buf[m++] = divr[j];
    } else if (comm_stage == 1) {
      buf[m++] = (double) status[j];
      buf[m++] = rsurface[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOSurface::unpack_forward_comm(int n, int first, double *buf)
{
  int i, k, a, b, m, last;
  int *status = atom->status;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (comm_stage == 0) {
      divr[i] = buf[m++];
    } else if (comm_stage == 1) {
      status[i] = (int) buf[m++];
      rsurface[i] = buf[m++];
    }
  }
}
