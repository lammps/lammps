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

#include "compute_rheo_grad.h"

#include "atom.h"
#include "comm.h"
#include "compute_rheo_kernel.h"
#include "compute_rheo_solids.h"
#include "domain.h"
#include "error.h"
#include "fix_rheo.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
enum{COMMGRAD, COMMFIELD};

/* ---------------------------------------------------------------------- */

ComputeRHEOGrad::ComputeRHEOGrad(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), fix_rheo(nullptr), compute_interface(nullptr), compute_kernel(nullptr),
  gradv(nullptr), gradr(nullptr), gradt(nullptr), gradn(nullptr)
{
  if (narg < 4) error->all(FLERR,"Illegal compute rheo/grad command");

  velocity_flag = temperature_flag = rho_flag = eta_flag = 0;
  for (int iarg = 3; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"velocity") == 0) velocity_flag = 1;
    else if (strcmp(arg[iarg],"rho") == 0) rho_flag = 1;
    else if (strcmp(arg[iarg],"temperature") == 0) temperature_flag = 1;
    else if (strcmp(arg[iarg],"viscosity") == 0) eta_flag = 1;
    else error->all(FLERR, "Illegal compute rheo/grad command, {}", arg[iarg]);
  }

  ncomm_grad = 0;
  ncomm_field = 0;
  comm_reverse = 0;

  int dim = domain->dimension;
  if (velocity_flag) {
    ncomm_grad += dim * dim;
    ncomm_field += dim;
    comm_reverse += dim * dim;
  }

  if (rho_flag) {
    ncomm_grad += dim;
    ncomm_field += 1;
    comm_reverse += dim;
  }

  if (temperature_flag) {
    ncomm_grad += dim;
    ncomm_field += 1;
    comm_reverse += dim;
  }

  if (eta_flag) {
    ncomm_grad += dim;
    comm_reverse += dim;
  }

  comm_forward = ncomm_grad;
}

/* ---------------------------------------------------------------------- */

ComputeRHEOGrad::~ComputeRHEOGrad()
{
  int dim = domain->dimension;
  int tmp1, tmp2, index;

  index = atom->find_custom("rheo_grad_v", tmp1, tmp2);
  if (index != 1) atom->remove_custom(index, 1, dim * dim);
  index = atom->find_custom("rheo_grad_rho", tmp1, tmp2);
  if (index != 1)  atom->remove_custom(index, 1, dim);
  index = atom->find_custom("rheo_grad_t", tmp1, tmp2);
  if (index != 1)  atom->remove_custom(index, 1, dim);
  index = atom->find_custom("rheo_grad_eta", tmp1, tmp2);
  if (index != 1)  atom->remove_custom(index, 1, dim);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOGrad::init()
{
  neighbor->add_request(this, NeighConst::REQ_DEFAULT);

  cut = fix_rheo->cut;
  cutsq = cut * cut;
  rho0 = fix_rheo->rho0;
  compute_kernel = fix_rheo->compute_kernel;
  compute_interface = fix_rheo->compute_interface;

  // Create coordination array if it doesn't already exist
  // Create a custom atom property so it works with compute property/atom
  // Do not create grow callback as there's no reason to copy/exchange data
  // Manually grow if nmax_old exceeded

  int index;
  int dim = domain->dimension;
  if (velocity_flag) {
    index = atom->add_custom("rheo_grad_v", 1, dim * dim);
    gradv = atom->darray[index];
  }

  if (rho_flag) {
    index = atom->add_custom("rheo_grad_rho", 1, dim);
    gradr = atom->darray[index];
  }

  if (temperature_flag) {
    index= atom->add_custom("rheo_grad_temp", 1, dim);
    gradt = atom->darray[index];
  }

  if (eta_flag) {
    index = atom->add_custom("rheo_grad_eta", 1, dim);
    gradn = atom->darray[index];
  }

  nmax_old = 0;
  grow_arrays(atom->nmax);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOGrad::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOGrad::compute_peratom()
{
  int i, j, k, ii, jj, jnum, itype, jtype, a, b;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq, imass, jmass;
  double rhoi, rhoj, Voli, Volj, drho, dT, deta;
  double vij[3];
  double wp, *dWij, *dWji;

  int inum, *ilist, *numneigh, **firstneigh;
  int *jlist;
  int nlocal = atom->nlocal;

  double **x = atom->x;
  double **v = atom->v;
  double *rho = atom->rho;
  double *temperature = atom->temperature;
  int *status = atom->status;
  int *type = atom->type;
  double *mass = atom->mass;
  int newton = force->newton;
  int dim = domain->dimension;

  int tmp1, tmp2;
  int index_visc = atom->find_custom("rheo_viscosity", tmp1, tmp2);
  if (index_visc == -1) error->all(FLERR, "Cannot find rheo viscosity");
  double *viscosity = atom->dvector[index_visc];

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // initialize arrays
  if (nmax > nmax_old) grow_arrays(nmax);

  for (i = 0; i < nmax; i++) {
    if (velocity_flag) {
      for (k = 0; k < dim * dim; k++)
        gradv[i][k] = 0.0;
    }
    if (rho_flag) {
      for (k = 0; k < dim; k++)
        gradr[i][k] = 0.0;
    }
    if (temperature_flag) {
      for (k = 0; k < dim; k++)
        gradt[i][k] = 0.0;
    }
    if (eta_flag) {
      for (k = 0; k < dim; k++)
        gradn[i][k] = 0.0;
    }
  }

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

      jtype = type[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutsq) {
        rhoi = rho[i];
        rhoj = rho[j];

        // Add corrections for walls
        if ((status[i] & FixRHEO::STATUS_FLUID) && !(status[j] & FixRHEO::STATUS_FLUID)) {
          compute_interface->correct_v(v[i], v[j], vi, i, j);
          rhoj = compute_interface->correct_rho(j,i);
        } else if (!(status[i] & FixRHEO::STATUS_FLUID) && (status[j] & FixRHEO::STATUS_FLUID)) {
          compute_interface->correct_v(v[j], v[i], vj, j, i);
          rhoi = compute_interface->correct_rho(i,j);
        } else if (!(status[i] & FixRHEO::STATUS_FLUID) && !(status[j] & FixRHEO::STATUS_FLUID)) {
          rhoi = rho0;
          rhoj = rho0;
        }

        Voli = mass[itype] / rhoi;
        Volj = mass[jtype] / rhoj;

        vij[0] = v[i][0] - v[j][0];
        vij[1] = v[i][1] - v[j][1];
        vij[2] = v[i][2] - v[j][2];

        if (rho_flag) drho = rhoi - rhoj;
        if (temperature_flag) dT = temperature[i] - temperature[j];
        if (eta_flag) deta = viscosity[i] - viscosity[j];

        wp = compute_kernel->calc_dw(i, j, delx, dely, delz, sqrt(rsq));
        dWij = compute_kernel->dWij;
        dWji = compute_kernel->dWji;

        for (a = 0; a < dim; a++) {
          for (b = 0; b < dim; b++) {
            if (velocity_flag) // uxx uxy uxz uyx uyy uyz uzx uzy uzz
              gradv[i][a * dim + b] -= vij[a] * Volj * dWij[b];
          }

          if (rho_flag) // P,x  P,y  P,z
            gradr[i][a] -= drho * Volj * dWij[a];

          if (temperature_flag) // T,x  T,y  T,z
            gradt[i][a] -= dT * Volj * dWij[a];

          if (eta_flag) // n,x  n,y  n,z
            gradn[i][a] -= deta * Volj * dWij[a];
        }

        if (newton || j < nlocal) {
          for (a = 0; a < dim; a++) {
            for (b = 0; b < dim; b++) {
              if (velocity_flag) // uxx uxy uxz uyx uyy uyz uzx uzy uzz
                gradv[j][a * dim + b] += vij[a] * Voli * dWji[b];
            }

            if (rho_flag) // P,x  P,y  P,z
              gradr[j][a] += drho * Voli * dWji[a];

            if (temperature_flag) // T,x  T,y  T,z
              gradt[j][a] += dT * Voli * dWji[a];

            if (eta_flag) // n,x  n,y  n,z
              gradn[j][a] += deta * Voli * dWji[a];
          }
        }
      }
    }
  }

  if (newton) comm->reverse_comm(this);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOGrad::forward_gradients()
{
  comm_stage = COMMGRAD;
  comm_forward = ncomm_grad;
  comm->forward_comm(this);
}


/* ---------------------------------------------------------------------- */

void ComputeRHEOGrad::forward_fields()
{
  comm_stage = COMMFIELD;
  comm_forward = ncomm_field;
  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

int ComputeRHEOGrad::pack_forward_comm(int n, int *list, double *buf,
                                        int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,k,m;
  double *rho = atom->rho;
  double *temperature = atom->temperature;
  double *eta = atom->viscosity;
  double **v = atom->v;
  int dim = domain->dimension;

  m = 0;

  for (i = 0; i < n; i++) {
    j = list[i];
    if (comm_stage == COMMGRAD) {

      if (velocity_flag)
        for (k = 0; k < dim * dim; k++)
          buf[m++] = gradv[j][k];

      if (rho_flag)
        for (k = 0; k < dim; k++)
          buf[m++] = gradr[j][k];

      if (temperature_flag)
        for (k = 0; k < dim; k++)
          buf[m++] = gradt[j][k];

      if (eta_flag)
        for (k = 0; k < dim; k++)
          buf[m++] = gradn[j][k];

    } else if (comm_stage == COMMFIELD) {

      if (velocity_flag)
        for (k = 0; k < dim; k++)
          buf[m++] = v[j][k];

      if (rho_flag)
        buf[m++] = rho[j];

      if (temperature_flag)
        buf[m++] = temperature[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOGrad::unpack_forward_comm(int n, int first, double *buf)
{
  int i, k, m, last;
  double * rho = atom->rho;
  double * temperature = atom->temperature;
  double ** v = atom->v;
  int dim = domain->dimension;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (comm_stage == COMMGRAD) {
      if (velocity_flag)
        for (k = 0; k < dim * dim; k++)
          gradv[i][k] = buf[m++];

      if (rho_flag)
        for (k = 0; k < dim; k++)
          gradr[i][k] = buf[m++];

      if (temperature_flag)
        for (k = 0; k < dim; k++)
          gradt[i][k] = buf[m++];

      if (eta_flag)
        for (k = 0; k < dim; k++)
          gradn[i][k] = buf[m++];

    } else if (comm_stage == COMMFIELD) {
      if (velocity_flag)
        for (k = 0; k < dim; k++)
          v[i][k] = buf[m++];

      if (rho_flag)
        rho[i] = buf[m++];

      if (temperature_flag)
        temperature[i] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int ComputeRHEOGrad::pack_reverse_comm(int n, int first, double *buf)
{
  int i,k,m,last;
  int dim = domain->dimension;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (velocity_flag)
      for (k = 0; k < dim * dim; k++)
        buf[m++] = gradv[i][k];

    if (rho_flag)
      for (k = 0; k < dim; k++)
        buf[m++] = gradr[i][k];

    if (temperature_flag)
      for (k = 0; k < dim; k++)
        buf[m++] = gradt[i][k];

    if (eta_flag)
      for (k = 0; k < dim; k++)
        buf[m++] = gradn[i][k];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOGrad::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,k,j,m;
  int dim = domain->dimension;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (velocity_flag)
      for (k = 0; k < dim * dim; k++)
        gradv[j][k] += buf[m++];

    if (rho_flag)
      for (k = 0; k < dim; k++)
        gradr[j][k] += buf[m++];

    if (temperature_flag)
      for (k = 0; k < dim; k++)
        gradt[j][k] += buf[m++];

    if (eta_flag)
      for (k = 0; k < dim; k++)
        gradn[j][k] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOGrad::grow_arrays(int nmax)
{
  int dim = domain->dimension;
  if (velocity_flag)
    memory->grow(gradv, nmax, dim * dim, "atom:rheo_grad_v");

  if (rho_flag)
    memory->grow(gradr, nmax, dim, "atom:rheo_grad_rho");

  if (temperature_flag)
    memory->grow(gradt, nmax, dim, "atom:rheo_grad_temp");

  if (eta_flag)
    memory->grow(gradn, nmax, dim, "atom:rheo_grad_eta");
  nmax_old = nmax;
}

/* ---------------------------------------------------------------------- */

double ComputeRHEOKernel::memory_usage()
{
  double bytes = 0.0;
  if (velocity_flag)
    bytes = (size_t) nmax_old * dim * dim * sizeof(double);

  if (rho_flag)
    bytes = (size_t) nmax_old * dim * sizeof(double);

  if (temperature_flag)
    bytes = (size_t) nmax_old * dim * sizeof(double);

  if (eta_flag)
    bytes = (size_t) nmax_old * dim * sizeof(double);

  return bytes;
}
