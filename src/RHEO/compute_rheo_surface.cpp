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
#include "compute_rheo_interface.h"
#include "compute_rheo_kernel.h"
#include "domain.h"
#include "error.h"
#include "fix_rheo.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"

using namespace LAMMPS_NS;
using namespace RHEO_NS;
using namespace FixConst;
using namespace MathExtra;

static constexpr double EPSILON = 1e-10;

/* ---------------------------------------------------------------------- */

ComputeRHEOSurface::ComputeRHEOSurface(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), nsurface(nullptr), rsurface(nullptr), divr(nullptr), fix_rheo(nullptr),
    rho0(nullptr), B(nullptr), gradC(nullptr), list(nullptr), compute_kernel(nullptr),
    compute_interface(nullptr)
{
  if (narg != 3) error->all(FLERR, "Illegal compute RHEO/SURFACE command");

  int dim = domain->dimension;
  comm_forward = 2;
  comm_reverse = dim * dim + 1;

  nmax_store = 0;
  grow_arrays(atom->nmax);
}

/* ---------------------------------------------------------------------- */

ComputeRHEOSurface::~ComputeRHEOSurface()
{
  memory->destroy(divr);
  memory->destroy(rsurface);
  memory->destroy(nsurface);
  memory->destroy(B);
  memory->destroy(gradC);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOSurface::init()
{
  compute_kernel = fix_rheo->compute_kernel;
  compute_interface = fix_rheo->compute_interface;
  cut = fix_rheo->cut;
  rho0 = fix_rheo->rho0;
  threshold_style = fix_rheo->surface_style;
  threshold_divr = fix_rheo->divr_surface;
  threshold_z = fix_rheo->zmin_surface;
  threshold_splash = fix_rheo->zmin_splash;
  interface_flag = fix_rheo->interface_flag;

  cutsq = cut * cut;

  // need an occasional half neighbor list
  neighbor->add_request(this, NeighConst::REQ_DEFAULT);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOSurface::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOSurface::compute_peratom()
{
  int i, j, ii, jj, inum, jnum, a, itype, jtype, fluidi, fluidj;
  double xtmp, ytmp, ztmp, rsq, Voli, Volj, rhoi, rhoj;
  double dWij[3], dWji[3], dx[3];
  int *ilist, *jlist, *numneigh, **firstneigh;

  int nlocal = atom->nlocal;

  double **x = atom->x;
  int *status = atom->rheo_status;
  int newton = force->newton;
  int dim = domain->dimension;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  double *rho = atom->rho;
  int *coordination = compute_kernel->coordination;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // Grow and zero arrays
  if (nmax_store < atom->nmax) grow_arrays(atom->nmax);

  size_t nbytes = nmax_store * sizeof(double);
  memset(&divr[0], 0, nbytes);
  memset(&rsurface[0], 0, nbytes);
  memset(&nsurface[0][0], 0, nbytes * dim);
  memset(&gradC[0][0], 0, nbytes * dim * dim);
  memset(&B[0][0], 0, nbytes * dim * dim);

  // loop over neighbors to calculate the average orientation of neighbors
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    itype = type[i];
    fluidi = !(status[i] & PHASECHECK);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      dx[0] = xtmp - x[j][0];
      dx[1] = ytmp - x[j][1];
      dx[2] = ztmp - x[j][2];

      rsq = lensq3(dx);
      if (rsq < cutsq) {
        jtype = type[j];
        fluidj = !(status[j] & PHASECHECK);

        rhoi = rho[i];
        rhoj = rho[j];

        // Add corrections for walls
        if (interface_flag) {
          if (fluidi && (!fluidj)) {
            rhoj = compute_interface->correct_rho(j);
          } else if ((!fluidi) && fluidj) {
            rhoi = compute_interface->correct_rho(i);
          } else if ((!fluidi) && (!fluidj)) {
            rhoi = rho0[itype];
            rhoj = rho0[jtype];
          }
        }

        if (rmass) {
          Voli = rmass[i] / rhoi;
          Volj = rmass[j] / rhoj;
        } else {
          Voli = mass[itype] / rhoi;
          Volj = mass[jtype] / rhoj;
        }
        compute_kernel->calc_dw_quintic(dx[0], dx[1], dx[2], sqrt(rsq), dWij, dWji);

        for (a = 0; a < dim; a++) {
          divr[i] -= dWij[a] * dx[a] * Volj;
          gradC[i][a] += dWij[a] * Volj;
        }

        if ((j < nlocal) || newton) {
          for (a = 0; a < dim; a++) {
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
      for (a = 0; a < dim; a++) maggC += gradC[i][a] * gradC[i][a];
      maggC = sqrt(maggC) + EPSILON;
      maggC = 1.0 / maggC;
      for (a = 0; a < dim; a++) nsurface[i][a] = -gradC[i][a] * maggC;
    }
  }

  // Remove surface settings and assign new values
  int nall = nlocal + atom->nghost;
  int test;

  for (i = 0; i < nall; i++) {
    status[i] &= SURFACEMASK;
    if (mask[i] & groupbit) {
      if (threshold_style == DIVR)
        test = divr[i] < threshold_divr;
      else
        test = coordination[i] < threshold_z;

      // Treat nonfluid particles as bulk
      if (status[i] & PHASECHECK) test = 0;

      if (test) {
        if (coordination[i] < threshold_splash)
          status[i] |= STATUS_SPLASH;
        else
          status[i] |= STATUS_SURFACE;
        rsurface[i] = 0.0;
      } else {
        status[i] |= STATUS_BULK;
        rsurface[i] = cut;
      }
    }
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fluidi = !(status[i] & PHASECHECK);

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      fluidj = !(status[j] & PHASECHECK);

      dx[0] = xtmp - x[j][0];
      dx[1] = ytmp - x[j][1];
      dx[2] = ztmp - x[j][2];
      rsq = lensq3(dx);
      if (rsq < cutsq) {
        if (fluidi) {
          if ((status[i] & STATUS_BULK) && (status[j] & STATUS_SURFACE)) {
            status[i] &= SURFACEMASK;
            status[i] |= STATUS_LAYER;
          }

          if (status[j] & STATUS_SURFACE) rsurface[i] = MIN(rsurface[i], sqrt(rsq));
        }

        if (fluidj && (j < nlocal || newton)) {
          if ((status[j] & STATUS_BULK) && (status[j] & PHASECHECK) &&
              (status[i] & STATUS_SURFACE)) {
            status[j] &= SURFACEMASK;
            status[j] |= STATUS_LAYER;
          }

          if (status[i] & STATUS_SURFACE) rsurface[j] = MIN(rsurface[j], sqrt(rsq));
        }
      }
    }
  }

  // clear normal vectors for non-surface particles

  for (i = 0; i < nall; i++) {
    if (mask[i] & groupbit) {
      if (!(status[i] & STATUS_SURFACE))
        for (a = 0; a < dim; a++) nsurface[i][a] = 0.0;
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
  int dim = domain->dimension;
  int *status = atom->rheo_status;
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) {
    if (comm_stage == 0) {
      buf[m++] = divr[i];
      for (int a = 0; a < dim; a++)
        for (int b = 0; b < dim; b++) buf[m++] = gradC[i][a * dim + b];
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
  int dim = domain->dimension;
  int *status = atom->rheo_status;
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    if (comm_stage == 0) {
      divr[j] += buf[m++];
      for (int a = 0; a < dim; a++)
        for (int b = 0; b < dim; b++) gradC[j][a * dim + b] += buf[m++];
    } else if (comm_stage == 1) {
      auto tmp1 = (int) buf[m++];
      if ((status[j] & STATUS_BULK) && (tmp1 & STATUS_LAYER)) {
        status[j] &= SURFACEMASK;
        status[j] |= STATUS_LAYER;
      }
      auto tmp2 = buf[m++];
      rsurface[j] = MIN(rsurface[j], tmp2);
    }
  }
}

/* ---------------------------------------------------------------------- */

int ComputeRHEOSurface::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                          int * /*pbc*/)
{
  int *status = atom->rheo_status;
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
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
  int *status = atom->rheo_status;
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) {
    if (comm_stage == 0) {
      divr[i] = buf[m++];
    } else if (comm_stage == 1) {
      status[i] = (int) buf[m++];
      rsurface[i] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOSurface::grow_arrays(int nmax)
{
  int dim = domain->dimension;

  memory->grow(divr, nmax, "rheo/surface:divr");
  memory->grow(rsurface, nmax, "rheo/surface:rsurface");
  memory->grow(nsurface, nmax, dim, "rheo/surface:nsurface");
  memory->grow(B, nmax, dim * dim, "rheo/surface:B");
  memory->grow(gradC, nmax, dim * dim, "rheo/surface:gradC");

  nmax_store = nmax;
}
