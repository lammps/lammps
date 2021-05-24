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

#include <cmath>
#include <cstdlib>
#include "atom_vec_dielectric.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecDielectric::AtomVecDielectric(LAMMPS *lmp) : AtomVecFull(lmp),
  mu(NULL), area(NULL), ed(NULL), em(NULL), q_real(NULL), epsilon(NULL),
  curvature(NULL)
{
  comm_x_only = 0;        // 0 to exchange other data in addition to x in pack comm
  size_forward = 13;      // # items in func: pack_comm
  size_reverse = 3;
  size_border = 18;       // # items in func: pack_border
  size_velocity = 3;
  size_data_atom = 15;    // # items from read data file, func: data_atom.
  size_data_vel = 4;
  xcol_data = 5;          // column in data files where coordinates start
                          // the first four columns are tag, mol, type, q

  atom->mu_flag = 1;
}

/* ---------------------------------------------------------------------- */

AtomVecDielectric::~AtomVecDielectric()
{
  memory->destroy(area);
  memory->destroy(ed);
  memory->destroy(em);
  memory->destroy(epsilon);
  memory->destroy(q_real);
  memory->destroy(curvature);
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecDielectric::grow(int n)
{
  AtomVecFull::grow(n);

  mu = memory->grow(atom->mu,nmax,4,"atom:mu");

  // USER-DIELECTRIC specifics

  area = memory->grow(area,nmax,"atom:area");
  ed = memory->grow(ed,nmax,"atom:ed");
  em = memory->grow(em,nmax,"atom:em");
  epsilon = memory->grow(epsilon,nmax,"atom:epsilon");
  q_real = memory->grow(q_real,nmax,"atom:q_real");
  curvature = memory->grow(curvature,nmax,"atom:curvature");
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecDielectric::grow_reset()
{
  AtomVecFull::grow_reset();

  mu = atom->mu;
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecDielectric::copy(int i, int j, int delflag)
{
  AtomVecFull::copy(i, j, delflag);

  // USER-DIELECTRIC specifics

  mu[j][0] = mu[i][0];
  mu[j][1] = mu[i][1];
  mu[j][2] = mu[i][2];
  mu[j][3] = mu[i][3];

  area[j] = area[i];
  ed[j] = ed[i];
  em[j] = em[i];
  epsilon[j] = epsilon[i];
  q_real[j] = q_real[i];
  curvature[j] = curvature[i];
}

/* ---------------------------------------------------------------------- */

int AtomVecDielectric::pack_comm(int n, int *list, double *buf,
                             int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = q[j];
      buf[m++] = mu[j][0];
      buf[m++] = mu[j][1];
      buf[m++] = mu[j][2];
      buf[m++] = area[j];
      buf[m++] = ed[j];
      buf[m++] = em[j];
      buf[m++] = epsilon[j];
      buf[m++] = q_real[j];
      buf[m++] = curvature[j];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = q[j];
      buf[m++] = mu[j][0];
      buf[m++] = mu[j][1];
      buf[m++] = mu[j][2];
      buf[m++] = area[j];
      buf[m++] = ed[j];
      buf[m++] = em[j];
      buf[m++] = epsilon[j];
      buf[m++] = q_real[j];
      buf[m++] = curvature[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDielectric::pack_comm_vel(int n, int *list, double *buf,
                                 int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = q[j];
      buf[m++] = mu[j][0];
      buf[m++] = mu[j][1];
      buf[m++] = mu[j][2];
      buf[m++] = area[j];
      buf[m++] = ed[j];
      buf[m++] = em[j];
      buf[m++] = epsilon[j];
      buf[m++] = q_real[j];
      buf[m++] = curvature[j];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = q[j];
        buf[m++] = mu[j][0];
        buf[m++] = mu[j][1];
        buf[m++] = mu[j][2];
        buf[m++] = area[j];
        buf[m++] = ed[j];
        buf[m++] = em[j];
        buf[m++] = epsilon[j];
        buf[m++] = q_real[j];
        buf[m++] = curvature[j];
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = q[j];
        buf[m++] = mu[j][0];
        buf[m++] = mu[j][1];
        buf[m++] = mu[j][2];
        buf[m++] = area[j];
        buf[m++] = ed[j];
        buf[m++] = em[j];
        buf[m++] = epsilon[j];
        buf[m++] = q_real[j];
        buf[m++] = curvature[j];
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
      }
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDielectric::pack_comm_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = q[j];
    buf[m++] = mu[j][0];
    buf[m++] = mu[j][1];
    buf[m++] = mu[j][2];
    buf[m++] = area[j];
    buf[m++] = ed[j];
    buf[m++] = em[j];
    buf[m++] = epsilon[j];
    buf[m++] = q_real[j];
    buf[m++] = curvature[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecDielectric::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    q[i] = buf[m++];
    mu[i][0] = buf[m++];
    mu[i][1] = buf[m++];
    mu[i][2] = buf[m++];
    area[i] = buf[m++];
    ed[i] = buf[m++];
    em[i] = buf[m++];
    epsilon[i] = buf[m++];
    q_real[i] = buf[m++];
    curvature[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecDielectric::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    q[i] = buf[m++];
    mu[i][0] = buf[m++];
    mu[i][1] = buf[m++];
    mu[i][2] = buf[m++];
    area[i] = buf[m++];
    ed[i] = buf[m++];
    em[i] = buf[m++];
    epsilon[i] = buf[m++];
    q_real[i] = buf[m++];
    curvature[i] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecDielectric::unpack_comm_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    q[i] = buf[m++];
    mu[i][0] = buf[m++];
    mu[i][1] = buf[m++];
    mu[i][2] = buf[m++];
    area[i] = buf[m++];
    ed[i] = buf[m++];
    em[i] = buf[m++];
    epsilon[i] = buf[m++];
    q_real[i] = buf[m++];
    curvature[i] = buf[m++];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDielectric::pack_reverse(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecDielectric::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecDielectric::pack_border(int n, int *list, double *buf,
                               int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = ubuf(molecule[j]).d;
      buf[m++] = q[j];
      buf[m++] = mu[j][0];
      buf[m++] = mu[j][1];
      buf[m++] = mu[j][2];
      buf[m++] = mu[j][3];
      buf[m++] = area[j];
      buf[m++] = ed[j];
      buf[m++] = em[j];
      buf[m++] = epsilon[j];
      buf[m++] = q_real[j];
      buf[m++] = curvature[j];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = ubuf(molecule[j]).d;
      buf[m++] = q[j];
      buf[m++] = mu[j][0];
      buf[m++] = mu[j][1];
      buf[m++] = mu[j][2];
      buf[m++] = mu[j][3];
      buf[m++] = area[j];
      buf[m++] = ed[j];
      buf[m++] = em[j];
      buf[m++] = epsilon[j];
      buf[m++] = q_real[j];
      buf[m++] = curvature[j];
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDielectric::pack_border_vel(int n, int *list, double *buf,
                                   int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = ubuf(molecule[j]).d;
      buf[m++] = q[j];
      buf[m++] = mu[j][0];
      buf[m++] = mu[j][1];
      buf[m++] = mu[j][2];
      buf[m++] = mu[j][3];
      buf[m++] = area[j];
      buf[m++] = ed[j];
      buf[m++] = em[j];
      buf[m++] = epsilon[j];
      buf[m++] = q_real[j];
      buf[m++] = curvature[j];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = ubuf(molecule[j]).d;
        buf[m++] = q[j];
        buf[m++] = mu[j][0];
        buf[m++] = mu[j][1];
        buf[m++] = mu[j][2];
        buf[m++] = mu[j][3];
        buf[m++] = area[j];
        buf[m++] = ed[j];
        buf[m++] = em[j];
        buf[m++] = epsilon[j];
        buf[m++] = q_real[j];
        buf[m++] = curvature[j];
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = ubuf(molecule[j]).d;
        buf[m++] = q[j];
        buf[m++] = mu[j][0];
        buf[m++] = mu[j][1];
        buf[m++] = mu[j][2];
        buf[m++] = mu[j][3];
        buf[m++] = area[j];
        buf[m++] = ed[j];
        buf[m++] = em[j];
        buf[m++] = epsilon[j];
        buf[m++] = curvature[j];
        buf[m++] = q_real[j];
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDielectric::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = ubuf(molecule[j]).d;
    buf[m++] = q[j];
    buf[m++] = mu[j][0];
    buf[m++] = mu[j][1];
    buf[m++] = mu[j][2];
    buf[m++] = mu[j][3];
    buf[m++] = area[j];
    buf[m++] = ed[j];
    buf[m++] = em[j];
    buf[m++] = epsilon[j];
    buf[m++] = q_real[j];
    buf[m++] = curvature[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecDielectric::unpack_border(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    molecule[i] = (tagint) ubuf(buf[m++]).i;
    q[i] = buf[m++];
    mu[i][0] = buf[m++];
    mu[i][1] = buf[m++];
    mu[i][2] = buf[m++];
    mu[i][3] = buf[m++];
    area[i] = buf[m++];
    ed[i] = buf[m++];
    em[i] = buf[m++];
    epsilon[i] = buf[m++];
    q_real[i] = buf[m++];
    curvature[i] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecDielectric::unpack_border_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    molecule[i] = (tagint) ubuf(buf[m++]).i;
    q[i] = buf[m++];
    mu[i][0] = buf[m++];
    mu[i][1] = buf[m++];
    mu[i][2] = buf[m++];
    mu[i][3] = buf[m++];
    area[i] = buf[m++];
    ed[i] = buf[m++];
    em[i] = buf[m++];
    epsilon[i] = buf[m++];
    q_real[i] = buf[m++];
    curvature[i] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

int AtomVecDielectric::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    molecule[i] = (tagint) ubuf(buf[m++]).i;
    q[i] = buf[m++];
    mu[i][0] = buf[m++];
    mu[i][1] = buf[m++];
    mu[i][2] = buf[m++];
    mu[i][3] = buf[m++];
    area[i] = buf[m++];
    ed[i] = buf[m++];
    em[i] = buf[m++];
    epsilon[i] = buf[m++];
    q_real[i] = buf[m++];
    curvature[i] = buf[m++];
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack all atom quantities for shipping to another proc
   xyz must be 1st 3 values, so that comm::exchange can test on them
------------------------------------------------------------------------- */

int AtomVecDielectric::pack_exchange(int i, double *buf)
{
  int k;

  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;

  buf[m++] = ubuf(molecule[i]).d;

  buf[m++] = q[i];
  buf[m++] = mu[i][0];
  buf[m++] = mu[i][1];
  buf[m++] = mu[i][2];
  buf[m++] = mu[i][3];
  buf[m++] = area[i];
  buf[m++] = ed[i];
  buf[m++] = em[i];
  buf[m++] = epsilon[i];
  buf[m++] = q_real[i];
  buf[m++] = curvature[i];

  buf[m++] = ubuf(num_bond[i]).d;
  for (k = 0; k < num_bond[i]; k++) {
    buf[m++] = ubuf(bond_type[i][k]).d;
    buf[m++] = ubuf(bond_atom[i][k]).d;
  }

  buf[m++] = ubuf(num_angle[i]).d;
  for (k = 0; k < num_angle[i]; k++) {
    buf[m++] = ubuf(angle_type[i][k]).d;
    buf[m++] = ubuf(angle_atom1[i][k]).d;
    buf[m++] = ubuf(angle_atom2[i][k]).d;
    buf[m++] = ubuf(angle_atom3[i][k]).d;
  }

  buf[m++] = ubuf(num_dihedral[i]).d;
  for (k = 0; k < num_dihedral[i]; k++) {
    buf[m++] = ubuf(dihedral_type[i][k]).d;
    buf[m++] = ubuf(dihedral_atom1[i][k]).d;
    buf[m++] = ubuf(dihedral_atom2[i][k]).d;
    buf[m++] = ubuf(dihedral_atom3[i][k]).d;
    buf[m++] = ubuf(dihedral_atom4[i][k]).d;
  }

  buf[m++] = ubuf(num_improper[i]).d;
  for (k = 0; k < num_improper[i]; k++) {
    buf[m++] = ubuf(improper_type[i][k]).d;
    buf[m++] = ubuf(improper_atom1[i][k]).d;
    buf[m++] = ubuf(improper_atom2[i][k]).d;
    buf[m++] = ubuf(improper_atom3[i][k]).d;
    buf[m++] = ubuf(improper_atom4[i][k]).d;
  }

  buf[m++] = ubuf(nspecial[i][0]).d;
  buf[m++] = ubuf(nspecial[i][1]).d;
  buf[m++] = ubuf(nspecial[i][2]).d;
  for (k = 0; k < nspecial[i][2]; k++) buf[m++] = ubuf(special[i][k]).d;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDielectric::unpack_exchange(double *buf)
{
  int k;

  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;

  molecule[nlocal] = (tagint) ubuf(buf[m++]).i;

  q[nlocal] = buf[m++];
  mu[nlocal][0] = buf[m++];
  mu[nlocal][1] = buf[m++];
  mu[nlocal][2] = buf[m++];
  mu[nlocal][3] = buf[m++];
  area[nlocal] = buf[m++];
  ed[nlocal] = buf[m++];
  em[nlocal] = buf[m++];
  epsilon[nlocal] = buf[m++];
  q_real[nlocal] = buf[m++];
  curvature[nlocal] = buf[m++];

  num_bond[nlocal] = (int) ubuf(buf[m++]).i;
  for (k = 0; k < num_bond[nlocal]; k++) {
    bond_type[nlocal][k] = (int) ubuf(buf[m++]).i;
    bond_atom[nlocal][k] = (tagint) ubuf(buf[m++]).i;
  }

  num_angle[nlocal] = (int) ubuf(buf[m++]).i;
  for (k = 0; k < num_angle[nlocal]; k++) {
    angle_type[nlocal][k] = (int) ubuf(buf[m++]).i;
    angle_atom1[nlocal][k] = (tagint) ubuf(buf[m++]).i;
    angle_atom2[nlocal][k] = (tagint) ubuf(buf[m++]).i;
    angle_atom3[nlocal][k] = (tagint) ubuf(buf[m++]).i;
  }

  num_dihedral[nlocal] = (int) ubuf(buf[m++]).i;
  for (k = 0; k < num_dihedral[nlocal]; k++) {
    dihedral_type[nlocal][k] = (int) ubuf(buf[m++]).i;
    dihedral_atom1[nlocal][k] = (tagint) ubuf(buf[m++]).i;
    dihedral_atom2[nlocal][k] = (tagint) ubuf(buf[m++]).i;
    dihedral_atom3[nlocal][k] = (tagint) ubuf(buf[m++]).i;
    dihedral_atom4[nlocal][k] = (tagint) ubuf(buf[m++]).i;
  }

  num_improper[nlocal] = (int) ubuf(buf[m++]).i;
  for (k = 0; k < num_improper[nlocal]; k++) {
    improper_type[nlocal][k] = (int) ubuf(buf[m++]).i;
    improper_atom1[nlocal][k] = (tagint) ubuf(buf[m++]).i;
    improper_atom2[nlocal][k] = (tagint) ubuf(buf[m++]).i;
    improper_atom3[nlocal][k] = (tagint) ubuf(buf[m++]).i;
    improper_atom4[nlocal][k] = (tagint) ubuf(buf[m++]).i;
  }

  nspecial[nlocal][0] = (int) ubuf(buf[m++]).i;
  nspecial[nlocal][1] = (int) ubuf(buf[m++]).i;
  nspecial[nlocal][2] = (int) ubuf(buf[m++]).i;
  for (k = 0; k < nspecial[nlocal][2]; k++)
    special[nlocal][k] = (tagint) ubuf(buf[m++]).i;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->
        unpack_exchange(nlocal,&buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecDielectric::size_restart()
{
  int i;

  int nlocal = atom->nlocal;
  int n = 0;
  // n is calculated from pack_restart: similar to full (17)
  //   17 = 3 for x[], 4 for (tag, type, mask, image), 3 for v[], 1 for molecule ID, 1 for q,
  //        4 for num bonds/angles/dihedrals/impropers and 1 for the first element (buf[0])
  // plus 10 = 4 for mu (norm) and 6 (area, ed, em, epsilon, q_real, curvature) 
  // totaling 27
  for (i = 0; i < nlocal; i++)
    n += 27 + 2*num_bond[i] + 4*num_angle[i] +
      5*num_dihedral[i] + 5*num_improper[i];

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
        n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVecDielectric::pack_restart(int i, double *buf)
{
  int k;

  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = ubuf(molecule[i]).d;

  buf[m++] = q[i];
  buf[m++] = mu[i][0];
  buf[m++] = mu[i][1];
  buf[m++] = mu[i][2];
  buf[m++] = mu[i][3];
  buf[m++] = area[i];
  buf[m++] = ed[i];
  buf[m++] = em[i];
  buf[m++] = epsilon[i];
  buf[m++] = q_real[i];
  buf[m++] = curvature[i];

  buf[m++] = ubuf(num_bond[i]).d;
  for (k = 0; k < num_bond[i]; k++) {
    buf[m++] = ubuf(MAX(bond_type[i][k],-bond_type[i][k])).d;
    buf[m++] = ubuf(bond_atom[i][k]).d;
  }

  buf[m++] = ubuf(num_angle[i]).d;
  for (k = 0; k < num_angle[i]; k++) {
    buf[m++] = ubuf(MAX(angle_type[i][k],-angle_type[i][k])).d;
    buf[m++] = ubuf(angle_atom1[i][k]).d;
    buf[m++] = ubuf(angle_atom2[i][k]).d;
    buf[m++] = ubuf(angle_atom3[i][k]).d;
  }

  buf[m++] = ubuf(num_dihedral[i]).d;
  for (k = 0; k < num_dihedral[i]; k++) {
    buf[m++] = ubuf(MAX(dihedral_type[i][k],-dihedral_type[i][k])).d;
    buf[m++] = ubuf(dihedral_atom1[i][k]).d;
    buf[m++] = ubuf(dihedral_atom2[i][k]).d;
    buf[m++] = ubuf(dihedral_atom3[i][k]).d;
    buf[m++] = ubuf(dihedral_atom4[i][k]).d;
  }

  buf[m++] = ubuf(num_improper[i]).d;
  for (k = 0; k < num_improper[i]; k++) {
    buf[m++] = ubuf(MAX(improper_type[i][k],-improper_type[i][k])).d;
    buf[m++] = ubuf(improper_atom1[i][k]).d;
    buf[m++] = ubuf(improper_atom2[i][k]).d;
    buf[m++] = ubuf(improper_atom3[i][k]).d;
    buf[m++] = ubuf(improper_atom4[i][k]).d;
  }

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecDielectric::unpack_restart(double *buf)
{
  int k;

  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra,nmax,atom->nextra_store,"atom:extra");
  }

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  molecule[nlocal] = (tagint) ubuf(buf[m++]).i;

  q[nlocal] = buf[m++];
  mu[nlocal][0] = buf[m++];
  mu[nlocal][1] = buf[m++];
  mu[nlocal][2] = buf[m++];
  mu[nlocal][3] = buf[m++];
  area[nlocal] = buf[m++];
  ed[nlocal] = buf[m++];
  em[nlocal] = buf[m++];
  epsilon[nlocal] = buf[m++];
  q_real[nlocal] = buf[m++];
  curvature[nlocal] = buf[m++];

  num_bond[nlocal] = (int) ubuf(buf[m++]).i;
  for (k = 0; k < num_bond[nlocal]; k++) {
    bond_type[nlocal][k] = (int) ubuf(buf[m++]).i;
    bond_atom[nlocal][k] = (tagint) ubuf(buf[m++]).i;
  }

  num_angle[nlocal] = (int) ubuf(buf[m++]).i;
  for (k = 0; k < num_angle[nlocal]; k++) {
    angle_type[nlocal][k] = (int) ubuf(buf[m++]).i;
    angle_atom1[nlocal][k] = (tagint) ubuf(buf[m++]).i;
    angle_atom2[nlocal][k] = (tagint) ubuf(buf[m++]).i;
    angle_atom3[nlocal][k] = (tagint) ubuf(buf[m++]).i;
  }

  num_dihedral[nlocal] = (int) ubuf(buf[m++]).i;
  for (k = 0; k < num_dihedral[nlocal]; k++) {
    dihedral_type[nlocal][k] = (int) ubuf(buf[m++]).i;
    dihedral_atom1[nlocal][k] = (tagint) ubuf(buf[m++]).i;
    dihedral_atom2[nlocal][k] = (tagint) ubuf(buf[m++]).i;
    dihedral_atom3[nlocal][k] = (tagint) ubuf(buf[m++]).i;
    dihedral_atom4[nlocal][k] = (tagint) ubuf(buf[m++]).i;
  }

  num_improper[nlocal] = (int) ubuf(buf[m++]).i;
  for (k = 0; k < num_improper[nlocal]; k++) {
    improper_type[nlocal][k] = (int) ubuf(buf[m++]).i;
    improper_atom1[nlocal][k] = (tagint) ubuf(buf[m++]).i;
    improper_atom2[nlocal][k] = (tagint) ubuf(buf[m++]).i;
    improper_atom3[nlocal][k] = (tagint) ubuf(buf[m++]).i;
    improper_atom4[nlocal][k] = (tagint) ubuf(buf[m++]).i;
  }

  nspecial[nlocal][0] = nspecial[nlocal][1] = nspecial[nlocal][2] = 0;

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecDielectric::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  mask[nlocal] = 1;
  image[nlocal] = ((imageint) IMGMAX << IMG2BITS) |
    ((imageint) IMGMAX << IMGBITS) | IMGMAX;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  q[nlocal] = 0.0;
  mu[nlocal][0] = 0.0;
  mu[nlocal][1] = 0.0;
  mu[nlocal][2] = 1.0;
  mu[nlocal][3] = 1.0;
  area[nlocal] = 1.0;
  ed[nlocal] = 0.0;
  em[nlocal] = 1.0;
  epsilon[nlocal] = 1.0;
  q_real[nlocal] = 0.0;
  curvature[nlocal] = 0.0;

  molecule[nlocal] = 0;
  num_bond[nlocal] = 0;
  num_angle[nlocal] = 0;
  num_dihedral[nlocal] = 0;
  num_improper[nlocal] = 0;
  nspecial[nlocal][0] = nspecial[nlocal][1] = nspecial[nlocal][2] = 0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecDielectric::data_atom(double *coord, imageint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = ATOTAGINT(values[0]);
  molecule[nlocal] = ATOTAGINT(values[1]);
  type[nlocal] = atoi(values[2]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  q_real[nlocal] = atof(values[3]);

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  mu[nlocal][0] = atof(values[7]);
  mu[nlocal][1] = atof(values[8]);
  mu[nlocal][2] = atof(values[9]);
  mu[nlocal][3] = 1;

  area[nlocal] = atof(values[10]);
  ed[nlocal] = atof(values[11]);
  em[nlocal] = atof(values[12]);
  epsilon[nlocal] = atof(values[13]);
  curvature[nlocal] = atof(values[14]);
  q[nlocal] = q_real[nlocal] / epsilon[nlocal];

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  num_bond[nlocal] = 0;
  num_angle[nlocal] = 0;
  num_dihedral[nlocal] = 0;
  num_improper[nlocal] = 0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecDielectric::data_atom_hybrid(int nlocal, char **values)
{
  molecule[nlocal] = ATOTAGINT(values[0]);
  q_real[nlocal] = atof(values[1]);
  mu[nlocal][0] = atof(values[2]);
  mu[nlocal][1] = atof(values[3]);
  mu[nlocal][2] = atof(values[4]);
  mu[nlocal][3] = 1;
  area[nlocal] = atof(values[5]);
  ed[nlocal] = atof(values[6]);
  em[nlocal] = atof(values[7]);
  epsilon[nlocal] = atof(values[8]);
  curvature[nlocal] = atof(values[9]);
  q[nlocal] = q_real[nlocal] / epsilon[nlocal];

  num_bond[nlocal] = 0;
  num_angle[nlocal] = 0;
  num_dihedral[nlocal] = 0;
  num_improper[nlocal] = 0;

  return 10;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecDielectric::pack_data(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = ubuf(molecule[i]).d;
    buf[i][2] = ubuf(type[i]).d;
    buf[i][3] = q[i];
    buf[i][4] = x[i][0];
    buf[i][5] = x[i][1];
    buf[i][6] = x[i][2];
    buf[i][7] = mu[i][0];
    buf[i][8] = mu[i][1];
    buf[i][9] = mu[i][2];
    buf[i][10] = area[i];
    buf[i][11] = ed[i];
    buf[i][12] = em[i];
    buf[i][13] = epsilon[i];
    buf[i][14] = curvature[i];
    buf[i][15] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
    buf[i][16] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][17] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;
  }
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecDielectric::pack_data_hybrid(int i, double *buf)
{
  buf[0] = ubuf(molecule[i]).d;
  buf[1] = q[i];
  buf[2] = mu[i][0];
  buf[3] = mu[i][1];
  buf[4] = mu[i][2];
  buf[5] = area[i];
  buf[6] = ed[i];
  buf[7] = em[i];
  buf[8] = epsilon[i];
  buf[9] = q_real[i];
  buf[10] = curvature[i];
  return 11;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecDielectric::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT " " TAGINT_FORMAT
//            " %d %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e "
//            "%-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %d %d %d\n",
            " %d %g %g %g %g %g %g "
            "%g %g %g %g %g %g %d %d %d\n",
            (tagint) ubuf(buf[i][0]).i,(tagint) ubuf(buf[i][1]).i, (int) ubuf(buf[i][2]).i,
            buf[i][3],buf[i][4],buf[i][5], buf[i][6],buf[i][7],buf[i][8],
            buf[i][9],buf[i][10],buf[i][11], buf[i][12], buf[i][13], buf[i][14],
            (int) ubuf(buf[i][15]).i,(int) ubuf(buf[i][16]).i,
            (int) ubuf(buf[i][17]).i);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecDielectric::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp, TAGINT_FORMAT
          " %-1.16e %-1.16e %-1.16e %-1.16e %-1.16 %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e",
    buf[0],buf[1],buf[2],buf[3],buf[4],buf[5],buf[6],buf[7],buf[8], buf[9], buf[10]);
  return 11;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecDielectric::memory_usage()
{
  bigint bytes = AtomVecFull::memory_usage();

  if (atom->memcheck("mu")) bytes += memory->usage(mu,nmax,4);

  if (atom->memcheck("area")) bytes += memory->usage(area,nmax);
  if (atom->memcheck("ed")) bytes += memory->usage(ed,nmax);
  if (atom->memcheck("em")) bytes += memory->usage(em,nmax);
  if (atom->memcheck("epsilon")) bytes += memory->usage(epsilon,nmax);
  if (atom->memcheck("q_real")) bytes += memory->usage(q_real,nmax);
  if (atom->memcheck("curvature")) bytes += memory->usage(curvature,nmax);

  return bytes;
}
