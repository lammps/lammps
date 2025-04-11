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
#include "compute_rheo_surface.h"
#include "domain.h"
#include "error.h"
#include "fix_rheo.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"

#include "update.h"

using namespace LAMMPS_NS;
using namespace RHEO_NS;
using namespace MathExtra;

/* ---------------------------------------------------------------------- */

ComputeRHEOVShift::ComputeRHEOVShift(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), vshift(nullptr), fix_rheo(nullptr), rho0(nullptr), wsame(nullptr),
    ct(nullptr), cgradt(nullptr), shift_type(nullptr), list(nullptr), compute_interface(nullptr),
    compute_kernel(nullptr), compute_surface(nullptr)
{
  if (narg != 3) error->all(FLERR, "Illegal compute RHEO/VShift command");

  comm_forward = 0;
  comm_reverse = 3;
  surface_flag = 0;

  nmax_store = 0;
}

/* ---------------------------------------------------------------------- */

ComputeRHEOVShift::~ComputeRHEOVShift()
{
  memory->destroy(vshift);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOVShift::init()
{
  neighbor->add_request(this, NeighConst::REQ_DEFAULT);

  surface_flag = fix_rheo->surface_flag;
  interface_flag = fix_rheo->interface_flag;

  compute_kernel = fix_rheo->compute_kernel;
  compute_interface = fix_rheo->compute_interface;
  compute_surface = fix_rheo->compute_surface;

  rho0 = fix_rheo->rho0;
  shift_type = fix_rheo->shift_type;
  cut = fix_rheo->cut;
  cutsq = cut * cut;
  cutthird = cut / 3.0;

  cross_type_flag = fix_rheo->shift_cross_type_flag;
  if (cross_type_flag) {
    scale = fix_rheo->shift_scale;
    wmin = fix_rheo->shift_wmin;
    cmin = fix_rheo->shift_cmin;
    comm_forward = 1;
    comm_reverse = 4;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOVShift::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOVShift::compute_peratom()
{
  int i, j, a, ii, jj, jnum, itype, jtype;
  int fluidi, fluidj;
  double xtmp, ytmp, ztmp, rsq, r, rinv;
  double w, wp, dr, w0, w4, vmag, prefactor;
  double imass, jmass, voli, volj, rhoi, rhoj;
  double dx[3], vi[3], vj[3];
  int dim = domain->dimension;

  int *jlist;
  int inum, *ilist, *numneigh, **firstneigh;

  int *type = atom->type;
  int *status = atom->rheo_status;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  double *rho = atom->rho;
  double *mass = atom->mass;
  double *rmass = atom->rmass;

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  if (nmax_store < atom->nmax) {
    memory->grow(vshift, atom->nmax, 3, "rheo:vshift");

    if (cross_type_flag) {
      memory->grow(ct, atom->nmax, "rheo:ct");
      memory->grow(cgradt, atom->nmax, 3, "rheo:cgradt");
      memory->grow(wsame, atom->nmax, "rheo:wsame");
    }

    nmax_store = atom->nmax;
  }

  for (i = 0; i < nall; i++)
    for (a = 0; a < dim; a++) vshift[i][a] = 0.0;

  for (a = 0; a < 3; a++) {
    vi[a] = 0.0;
    vj[a] = 0.0;
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    if (rmass)
      imass = rmass[i];
    else
      imass = mass[itype];
    fluidi = !(status[i] & PHASECHECK);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      fluidj = !(status[j] & PHASECHECK);
      if ((!fluidi) && (!fluidj)) continue;

      // Will skip shifting in FixRHEO initial integrate, but also skip here to save time
      if ((status[i] & STATUS_NO_SHIFT) && (status[j] & STATUS_NO_SHIFT)) continue;

      dx[0] = xtmp - x[j][0];
      dx[1] = ytmp - x[j][1];
      dx[2] = ztmp - x[j][2];
      rsq = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      if (rsq < cutsq) {
        jtype = type[j];
        if (rmass)
          jmass = rmass[j];
        else
          jmass = mass[jtype];

        r = sqrt(rsq);
        rinv = 1 / r;

        for (a = 0; a < dim; a++) {
          vi[a] = v[i][a];
          vj[a] = v[j][a];
        }

        rhoi = rho[i];
        rhoj = rho[j];

        // Add corrections for walls
        if (interface_flag) {
          if (fluidi && (!fluidj)) {
            compute_interface->correct_v(vj, vi, j, i);
            rhoj = compute_interface->correct_rho(j);
          } else if ((!fluidi) && fluidj) {
            compute_interface->correct_v(vi, vj, i, j);
            rhoi = compute_interface->correct_rho(i);
          } else if ((!fluidi) && (!fluidj)) {
            rhoi = rho0[itype];
            rhoj = rho0[jtype];
          }
        }

        voli = imass / rhoi;
        volj = jmass / rhoj;

        wp = compute_kernel->calc_dw(i, j, dx[0], dx[1], dx[2], r);
        w = compute_kernel->calc_w(i, j, dx[0], dx[1], dx[2], r);
        w0 = compute_kernel->calc_w(i, j, 0, 0, 0, cutthird);    // dx, dy, dz irrelevant
        w4 = w * w * w * w / (w0 * w0 * w0 * w0);
        dr = -2 * cutthird * (1 + 0.2 * w4) * wp * rinv;

        if ((mask[i] & groupbit) && fluidi) {
          vmag = sqrt(vi[0] * vi[0] + vi[1] * vi[1] + vi[2] * vi[2]);
          prefactor = vmag * volj * dr;

          vshift[i][0] += prefactor * dx[0];
          vshift[i][1] += prefactor * dx[1];
          vshift[i][2] += prefactor * dx[2];
        }

        if (newton_pair || j < nlocal) {
          if ((mask[j] & groupbit) && fluidj) {
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

  comm_stage = 0;
  if (newton_pair) comm->reverse_comm(this, 3);

  // Zero any excluded types

  for (i = 0; i < nlocal; i++)
    if (!shift_type[type[i]])
      for (a = 0; a < dim; a++)
        vshift[i][a] = 0.0;

  if (cross_type_flag) correct_type_interface();
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOVShift::correct_surfaces()
{
  if (!surface_flag) return;

  int *status = atom->rheo_status;
  int *mask = atom->mask;
  double **nsurface = compute_surface->nsurface;

  int nlocal = atom->nlocal;
  int dim = domain->dimension;

  double nx, ny, nz, vx, vy, vz, dot;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      if (status[i] & PHASECHECK) continue;

      if (status[i] & STATUS_SURFACE) {
        nx = nsurface[i][0];
        ny = nsurface[i][1];
        vx = vshift[i][0];
        vy = vshift[i][1];

        dot = nx * vx + ny * vy;
        if (dim == 3) {
          nz = nsurface[i][2];
          vz = vshift[i][2];
          dot += nz * vz;
        }

        // Allowing shifting into the bulk
        if (dot < 0.0) continue;

        vshift[i][0] = (1 - nx * nx) * vx - nx * ny * vy;
        vshift[i][1] = (1 - ny * ny) * vy - nx * ny * vx;
        if (dim == 3) {
          vshift[i][0] -= nx * nz * vz;
          vshift[i][1] -= ny * nz * vz;
          vshift[i][2] = (1 - nz * nz) * vz - nz * ny * vy - nx * nz * vx;
        } else {
          vshift[i][2] = 0.0;
        }
      } else if (status[i] & STATUS_SPLASH) {
        vshift[i][0] = 0.0;
        vshift[i][1] = 0.0;
        vshift[i][2] = 0.0;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOVShift::correct_type_interface()
{
  int i, j, a, ii, jj, jnum, itype, jtype;
  int fluidi, fluidj;
  double xtmp, ytmp, ztmp, rsq, r, w;
  double imass, jmass, voli, volj, rhoi, rhoj;
  double dx[3];
  int dim = domain->dimension;

  int *jlist;
  int inum, *ilist, *numneigh, **firstneigh;

  int *type = atom->type;
  int *status = atom->rheo_status;
  double **x = atom->x;
  double *rho = atom->rho;
  double *mass = atom->mass;
  double *rmass = atom->rmass;

  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  size_t nbytes = nmax_store * sizeof(double);
  memset(&ct[0], 0, nbytes);
  memset(&wsame[0], 0, nbytes);
  memset(&cgradt[0][0], 0, 3 * nbytes);
  double ctmp, *dWij, *dWji;

  // Calculate color gradient

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    fluidi = !(status[i] & PHASECHECK);
    jlist = firstneigh[i];
    jnum = numneigh[i];
    if (rmass)
      imass = rmass[i];
    else
      imass = mass[itype];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      dx[0] = xtmp - x[j][0];
      dx[1] = ytmp - x[j][1];
      dx[2] = ztmp - x[j][2];

      rsq = lensq3(dx);

      if (rsq > cutsq) continue;

      fluidj = !(status[j] & PHASECHECK);
      jtype = type[j];
      if (rmass)
        jmass = rmass[j];
      else
        jmass = mass[jtype];
      r = sqrt(rsq);

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

      voli = imass / rhoi;
      volj = jmass / rhoj;

      w = compute_kernel->calc_w(i, j, dx[0], dx[1], dx[2], r);

      if (itype != jtype) ctmp = 1;
      else ctmp = 0;

      ct[i] += ctmp * volj * w;
      if (newton_pair || j < nlocal)
        ct[j] += ctmp * voli * w;
    }
  }

  comm_stage = 1;
  if (newton_pair) comm->reverse_comm(this, 1);

  // Calculate color gradient
  // Note: in future might want to generalize this so color function can be used
  //   by other calculations (e.g. surface tension)
  //   maybe can create custom "calc_grad" method that takes an arbitrary field
  //   in ComputeRHEOGrad?

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    fluidi = !(status[i] & PHASECHECK);
    jlist = firstneigh[i];
    jnum = numneigh[i];
    imass = mass[itype];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      dx[0] = xtmp - x[j][0];
      dx[1] = ytmp - x[j][1];
      dx[2] = ztmp - x[j][2];
      rsq = lensq3(dx);

      if (rsq > cutsq) continue;

      fluidj = !(status[j] & PHASECHECK);
      jtype = type[j];
      if (rmass)
        jmass = rmass[j];
      else
        jmass = mass[jtype];
      r = sqrt(rsq);

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

      voli = imass / rhoi;
      volj = jmass / rhoj;

      w = compute_kernel->calc_w(i, j, dx[0], dx[1], dx[2], r);
      dWij = compute_kernel->dWij;
      dWji = compute_kernel->dWji;

      if (itype != jtype) ctmp = 1;
      else ctmp = 0;

      for (a = 0; a < dim; a++) {
        cgradt[i][a] -= ctmp * volj * dWij[a];
        if (newton_pair || j < nlocal)
          cgradt[j][a] -= ctmp * voli * dWji[a];
      }

      if (itype == jtype) {
        wsame[i] += w * r;
        if (newton_pair || j < nlocal)
          wsame[j] += w * r;
      }
    }
  }

  comm_stage = 2;
  if (newton_pair) comm->reverse_comm(this, 4);
  comm->forward_comm(this, 1);

  // Correct shifting at fluid-fluid interface
  // remove normal shifting component for interfacial particles
  // Based on Yang, Rakhsha, Hu, & Negrut 2022

  double ntmp[3], minv, dot;

  for (i = 0; i < nlocal; i++) {

    // If isolated, just don't shift
    if (wsame[i] < wmin) {
      for (a = 0; a < dim; a++)
        vshift[i][a] = 0.0;
      continue;
    }

    if (ct[i] < cmin) continue;

    minv = 0;
    for (a = 0; a < dim; a++)
      minv += cgradt[i][a] * cgradt[i][a];

    if (minv != 0)
      minv = 1 / sqrt(minv);

    for (a = 0; a < dim; a++)
      ntmp[a] = cgradt[i][a] * minv;

    dot = 0.0;
    for (a = 0; a < dim; a++)
      dot += ntmp[a] * vshift[i][a];

    // To allowing shifting into the same phase bulk
    // if (dot > 0.0) continue;

    for (a = 0; a < dim; a++)
      vshift[i][a] -= (1.0 - scale) * ntmp[a] * dot;
  }
}

/* ---------------------------------------------------------------------- */

int ComputeRHEOVShift::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;
  m = 0;

  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = wsame[j];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOVShift::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    wsame[i] = buf[m++];
}

/* ---------------------------------------------------------------------- */

int ComputeRHEOVShift::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, a, last;

  m = 0;
  last = first + n;
  if (comm_stage == 0) {
    for (i = first; i < last; i++) {
      buf[m++] = vshift[i][0];
      buf[m++] = vshift[i][1];
      buf[m++] = vshift[i][2];
    }
  } else if (comm_stage == 1) {
    for (i = first; i < last; i++)
      buf[m++] = ct[i];
  } else {
    for (i = first; i < last; i++) {
      for (a = 0; a < 3; a++)
        buf[m++] = cgradt[i][a];
      buf[m++] = wsame[i];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOVShift::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, a, m;

  m = 0;
  if (comm_stage == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      vshift[j][0] += buf[m++];
      vshift[j][1] += buf[m++];
      vshift[j][2] += buf[m++];
    }
  } else if (comm_stage == 1) {
    for (i = 0; i < n; i++) {
      j = list[i];
      ct[j] += buf[m++];
    }
  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      for (a = 0; a < 3; a++)
        cgradt[j][a] += buf[m++];
      wsame[j] += buf[m++];
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeRHEOVShift::memory_usage()
{
  double bytes = 3 * nmax_store * sizeof(double);

  if (cross_type_flag)
    bytes += 5 * nmax_store * sizeof(double);

  return bytes;
}
