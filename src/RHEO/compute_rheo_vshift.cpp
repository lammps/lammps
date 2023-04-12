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

#include "compute_rheo_vshift.h"

#include "atom.h"
#include "comm.h"
#include "compute_rheo_interface.h"
#include "compute_rheo_kernel.h"
#include "domain.h"
#include "error.h"
#include "fix_rheo.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeRHEOVShift::ComputeRHEOVShift(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), vshift(nullptr), fix_rheo(nullptr), fix_rheo(nullptr),
  compute_kernel(nullptr), compute_interface(nullptr)
{
  if (narg != 3) error->all(FLERR,"Illegal compute RHEO/VShift command");

  comm_reverse = 3;
  surface_flag = 0;

  // Create vshift array if it doesn't already exist
  // Create a custom atom property so it works with compute property/atom
  // Do not create grow callback as there's no reason to copy/exchange data
  // Manually grow if nmax_old exceeded

  int tmp1, tmp2;
  int index = atom->find_custom("rheo_vshift", tmp1, tmp2);
  if (index == -1) {
    index = atom->add_custom("rheo_vshift", 1, 3);
    nmax_old = atom->nmax;
  }
  vshift = atom->dvector[index];
}

/* ---------------------------------------------------------------------- */

ComputeRHEOVShift::~ComputeRHEOVShift()
{
  // Remove custom property if it exists
  int tmp1, tmp2, index;
  index = atom->find_custom("rheo_vshift", tmp1, tmp2);
  if (index != -1) atom->remove_custom(index, 1, 3);

}

/* ---------------------------------------------------------------------- */

void ComputeRHEOVShift::init()
{
  neighbor->add_request(this, NeighConst::REQ_DEFAULT);

  surface_flag = 0;
  if (fix_rheo->surface_flag)
    surface_flag = 1;

  compute_kernel = fix_rheo->compute_kernel;
  compute_interface = fix_rheo->compute_interface;

  cut = fix_rheo->cut;
  cutsq = cut * cut;
  cutthird = cut / 3.0;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOVShift::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOVShift::compute_peratom()
{
  int i, j, a, b, ii, jj, jnum, itype, jtype;
  int fluidi, fluidj;
  double xtmp, ytmp, ztmp, rsq, r, rinv;
  double w, wp, dr, w0, w4, vmag, prefactor;
  double imass, jmass, voli, volj, rhoi, rhoj;
  double dx[3], vi[3], vj[3] = {0};
  int dim = domain->dimension;

  int *jlist;
  int inum, *ilist, *numneigh, **firstneigh;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;
  int *status = atom->status;
  int *surface = atom->surface;
  double *rho = atom->rho;
  double *mass = atom->mass;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  if (nmax_old < atom->nmax)
    memory->grow(vshift, atom->nmax, 3, "atom:rheo_vshift");
  nmax_old = atom->nmax;

  for (i = 0; i < nall; i++)
    for (a = 0; a < dim; a++)
      vshift[i][a] = 0.0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    imass = mass[itype];
    fluidi = status[i] & FixRHEO::STATUS_FLUID;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      fluidj = status[j] & FixRHEO::STATUS_FLUID;
      if ((!fluidi) && (!fluidj)) continue;
      if (!(status[i] & FixRHEO::STATUS_SHIFT) && !(status[j] & FixRHEO::STATUS_SHIFT)) continue;

      dx[0] = xtmp - x[j][0];
      dx[1] = ytmp - x[j][1];
      dx[2] = ztmp - x[j][2];
      rsq = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      if (rsq < cutsq) {
        jtype = type[j];
        jmass = mass[jtype];

        r = sqrt(rsq);
        rinv = 1 / r;

        for (a = 0; a < dim; a ++) {
          vi[a] = v[i][a];
          vj[a] = v[j][a];
        }

        rhoi = rho[i];
        rhoj = rho[j];

        // Add corrections for walls
        if (fluidi && (!fluidj)) {
          compute_interface->correct_v(v[i], v[j], vi, i, j);
          rhoj = compute_interface->correct_rho(j,i);
        } else if ((!fluidi) && fluidj) {
          compute_interface->correct_v(v[j], v[i], vj, j, i);
          rhoi = compute_interface->correct_rho(i,j);
        } else if ((!fluidi) && (!fluidj)) {
          rhoi = 1.0;
          rhoj = 1.0;
        }

        voli = imass / rhoi;
        volj = jmass / rhoj;

        wp = compute_kernel->calc_dw(i, j, dx[0], dx[1], dx[2], r);
        w = compute_kernel->calc_w(i, j, dx[0], dx[1], dx[2], r);
        w0 = compute_kernel->calc_w(i, j, 0, 0, 0, cutthird); // dx, dy, dz irrelevant
        w4 = w * w * w * w / (w0 * w0 * w0 * w0);
        dr = -2 * cutthird * (1 + 0.2 * w4) * wp * rinv;

        if (mask[i] & groupbit) {
          vmag = sqrt(vi[0] * vi[0] + vi[1] * vi[1] + vi[2] * vi[2]);
          prefactor = vmag * volj * dr;
          vshift[i][0] += prefactor * dx[0];
          vshift[i][1] += prefactor * dx[1];
          vshift[i][2] += prefactor * dx[2];
        }

        if (newton_pair || j < nlocal) {
          if (mask[j] & groupbit) {
            vmag = sqrt(vj[0] * vj[0] + vj[1] * vj[1] + vj[2] * vj[2]);
            prefactor = vmag * voli * dr;
            vshift[j][0] -= prefactor * dx[0];
            vshift[j][1] -= prefactor * dx[1];
            vshift[j][2] -= prefactor * dx[2];
          }
        }
      }
    }
  }

  if (newton_pair) comm->reverse_comm_compute(this);
}


/* ---------------------------------------------------------------------- */

void ComputeRHEOVShift::correct_surfaces()
{
  if (!surface_flag) return;

  int *status = atom->status;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int i, a, b;
  int dim = domain->dimension;

  int tmp1, tmp2;
  define after surf
  int index_nsurf = atom->find_custom("rheo_nsurf", tmp1, tmp2);
  if (index_nsurf == -1) error->all(FLERR, "Cannot find rheo nsurf");
  double **nsurf = atom->darray[index_nsurf];

  double nx,ny,nz,vx,vy,vz;
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if ((surface[i] & FixRHEO::STATUS_SURFACE) || (surface[i] & FixRHEO::STATUS_LAYER)) {
        nx = nsurf[i][0];
        ny = nsurf[i][1];
        vx = vshift[i][0];
        vy = vshift[i][1];
        vz = vshift[i][2];
        vshift[i][0] = (1 - nx * nx) * vx - nx * ny * vy;
        vshift[i][1] = (1 - ny * ny) * vy - nx * ny * vx;
        if (dim > 2) {
          nz = nsurf[i][2];
          vshift[i][0] -= nx * nz * vz;
          vshift[i][1] -= ny * nz * vz;
          vshift[i][2] = (1 - nz * nz) * vz - nz * ny * vy - nx * nz * vx;
        } else {
          vshift[i][2] = 0.0;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int ComputeRHEOVShift::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = vshift[i][0];
    buf[m++] = vshift[i][1];
    buf[m++] = vshift[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOVShift::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    vshift[j][0] += buf[m++];
    vshift[j][1] += buf[m++];
    vshift[j][2] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeRHEOVShift::memory_usage()
{
  double bytes = 3 * nmax_old * sizeof(double);
  return bytes;
}
