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

/* ------------------------------------------------------------------------
   Contributing authors: Thomas Swinburne (CNRS & CINaM, Marseille, France)

   Please cite the related publication:
   T.D. Swinburne and M.-C. Marinica, Unsupervised calculation of free energy barriers in large crystalline systems, Physical Review Letters 2018
------------------------------------------------------------------------- */
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "atom_vec_pafi.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "memory.h"
#include "modify.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecPAFI::AtomVecPAFI(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;
  mass_type = 1;
  comm_x_only = 1;
  comm_f_only = 1;
  size_forward = 3;
  size_reverse = 3;
  size_border = 15;
  size_velocity = 3;
  size_data_atom = 5;//14;
  size_data_vel = 4;
  xcol_data = 3;
  atom->pafi_flag = 1;
}


/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecPAFI::grow(int n)
{
  if (n == 0) grow_nmax();
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  tag = memory->grow(atom->tag,nmax,"atom:tag");
  type = memory->grow(atom->type,nmax,"atom:type");
  mask = memory->grow(atom->mask,nmax,"atom:mask");
  image = memory->grow(atom->image,nmax,"atom:image");

  // allocating mech. quantities

  x = memory->grow(atom->x,nmax,3,"atom:x");
  v = memory->grow(atom->v,nmax,3,"atom:v");
  f = memory->grow(atom->f,nmax*comm->nthreads,3,"atom:f");

  // allocating path quantities
  path = memory->grow(atom->path,nmax,3,"atom:path");
  norm = memory->grow(atom->norm,nmax,3,"atom:norm");
  dnorm = memory->grow(atom->dnorm,nmax,3,"atom:dnorm");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecPAFI::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;
  path = atom->path; norm = atom->norm; dnorm = atom->dnorm;
}


/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecPAFI::copy(int i, int j, int delflag)
{
  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];

  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];

  path[j][0] = path[i][0];
  path[j][1] = path[i][1];
  path[j][2] = path[i][2];

  norm[j][0] = norm[i][0];
  norm[j][1] = norm[i][1];
  norm[j][2] = norm[i][2];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ---------------------------------------------------------------------- */

int AtomVecPAFI::pack_comm(int n, int *list, double *buf,
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
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecPAFI::pack_comm_vel(int n, int *list, double *buf,
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

void AtomVecPAFI::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecPAFI::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecPAFI::pack_reverse(int n, int first, double *buf)
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

void AtomVecPAFI::unpack_reverse(int n, int *list, double *buf)
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

int AtomVecPAFI::pack_border(int n, int *list, double *buf,
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
      buf[m++] = path[j][0];
      buf[m++] = path[j][1];
      buf[m++] = path[j][2];
      buf[m++] = norm[j][0];
      buf[m++] = norm[j][1];
      buf[m++] = norm[j][2];
      buf[m++] = dnorm[j][0];
      buf[m++] = dnorm[j][1];
      buf[m++] = dnorm[j][2];
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
      buf[m++] = path[j][0];
      buf[m++] = path[j][1];
      buf[m++] = path[j][2];
      buf[m++] = norm[j][0];
      buf[m++] = norm[j][1];
      buf[m++] = norm[j][2];
      buf[m++] = dnorm[j][0];
      buf[m++] = dnorm[j][1];
      buf[m++] = dnorm[j][2];
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecPAFI::pack_border_vel(int n, int *list, double *buf,
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
      buf[m++] = path[j][0];
      buf[m++] = path[j][1];
      buf[m++] = path[j][2];
      buf[m++] = norm[j][0];
      buf[m++] = norm[j][1];
      buf[m++] = norm[j][2];
      buf[m++] = dnorm[j][0];
      buf[m++] = dnorm[j][1];
      buf[m++] = dnorm[j][2];
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
        buf[m++] = path[j][0];
        buf[m++] = path[j][1];
        buf[m++] = path[j][2];
        buf[m++] = norm[j][0];
        buf[m++] = norm[j][1];
        buf[m++] = norm[j][2];
        buf[m++] = dnorm[j][0];
        buf[m++] = dnorm[j][1];
        buf[m++] = dnorm[j][2];
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
        buf[m++] = path[j][0];
        buf[m++] = path[j][1];
        buf[m++] = path[j][2];
        buf[m++] = norm[j][0];
        buf[m++] = norm[j][1];
        buf[m++] = norm[j][2];
        buf[m++] = dnorm[j][0];
        buf[m++] = dnorm[j][1];
        buf[m++] = dnorm[j][2];
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

int AtomVecPAFI::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = path[j][0];
    buf[m++] = path[j][1];
    buf[m++] = path[j][2];
    buf[m++] = norm[j][0];
    buf[m++] = norm[j][1];
    buf[m++] = norm[j][2];
    buf[m++] = dnorm[j][0];
    buf[m++] = dnorm[j][1];
    buf[m++] = dnorm[j][2];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecPAFI::unpack_border(int n, int first, double *buf)
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
    path[i][0] = buf[m++];
    path[i][1] = buf[m++];
    path[i][2] = buf[m++];
    norm[i][0] = buf[m++];
    norm[i][1] = buf[m++];
    norm[i][2] = buf[m++];
    dnorm[i][0] = buf[m++];
    dnorm[i][1] = buf[m++];
    dnorm[i][2] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);

}

/* ---------------------------------------------------------------------- */

void AtomVecPAFI::unpack_border_vel(int n, int first, double *buf)
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
    path[i][0] = buf[m++];
    path[i][1] = buf[m++];
    path[i][2] = buf[m++];
    norm[i][0] = buf[m++];
    norm[i][1] = buf[m++];
    norm[i][2] = buf[m++];
    dnorm[i][0] = buf[m++];
    dnorm[i][1] = buf[m++];
    dnorm[i][2] = buf[m++];
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

int AtomVecPAFI::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    path[i][0] = buf[m++];
    path[i][1] = buf[m++];
    path[i][2] = buf[m++];
    norm[i][0] = buf[m++];
    norm[i][1] = buf[m++];
    norm[i][2] = buf[m++];
    dnorm[i][0] = buf[m++];
    dnorm[i][1] = buf[m++];
    dnorm[i][2] = buf[m++];
  }

  return m;
}

/* ----------------------------------------------------------------------
   pack all atom quantities for shipping to another proc
   xyz must be 1st 3 values, so that comm::exchange can test on them
------------------------------------------------------------------------- */

int AtomVecPAFI::pack_exchange(int i, double *buf)
{
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

  buf[m++] = path[i][0];
  buf[m++] = path[i][1];
  buf[m++] = path[i][2];
  buf[m++] = norm[i][0];
  buf[m++] = norm[i][1];
  buf[m++] = norm[i][2];
  buf[m++] = dnorm[i][0];
  buf[m++] = dnorm[i][1];
  buf[m++] = dnorm[i][2];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecPAFI::unpack_exchange(double *buf)
{
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

  path[nlocal][0] = buf[m++];
  path[nlocal][1] = buf[m++];
  path[nlocal][2] = buf[m++];
  norm[nlocal][0] = buf[m++];
  norm[nlocal][1] = buf[m++];
  norm[nlocal][2] = buf[m++];
  dnorm[nlocal][0] = buf[m++];
  dnorm[nlocal][1] = buf[m++];
  dnorm[nlocal][2] = buf[m++];

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

int AtomVecPAFI::size_restart()
{
  int i;

  int nlocal = atom->nlocal;
  int n = 20 * nlocal;

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

int AtomVecPAFI::pack_restart(int i, double *buf)
{
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

  buf[m++] = path[i][0];
  buf[m++] = path[i][1];
  buf[m++] = path[i][2];
  buf[m++] = norm[i][0];
  buf[m++] = norm[i][1];
  buf[m++] = norm[i][2];
  buf[m++] = dnorm[i][0];
  buf[m++] = dnorm[i][1];
  buf[m++] = dnorm[i][2];

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecPAFI::unpack_restart(double *buf)
{
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

  path[nlocal][0] = buf[m++];
  path[nlocal][1] = buf[m++];
  path[nlocal][2] = buf[m++];
  norm[nlocal][0] = buf[m++];
  norm[nlocal][1] = buf[m++];
  norm[nlocal][2] = buf[m++];
  dnorm[nlocal][0] = buf[m++];
  dnorm[nlocal][1] = buf[m++];
  dnorm[nlocal][2] = buf[m++];

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

void AtomVecPAFI::create_atom(int itype, double *coord)
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

  path[nlocal][0] = 0.0;
  path[nlocal][1] = 0.0;
  path[nlocal][2] = 0.0;
  norm[nlocal][0] = 0.0;
  norm[nlocal][1] = 0.0;
  norm[nlocal][2] = 0.0;
  dnorm[nlocal][0] = 0.0;
  dnorm[nlocal][1] = 0.0;
  dnorm[nlocal][2] = 0.0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecPAFI::data_atom(double *coord, imageint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = ATOTAGINT(values[0]);
  type[nlocal] = atoi(values[1]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  path[nlocal][0] = 0.;
  path[nlocal][1] = 0.;
  path[nlocal][2] = 0.;

  norm[nlocal][0] = 0.;
  norm[nlocal][1] = 0.;
  norm[nlocal][2] = 0.;

  dnorm[nlocal][0] = 0.;
  dnorm[nlocal][1] = 0.;
  dnorm[nlocal][2] = 0.;

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecPAFI::data_atom_hybrid(int nlocal, char **values)
{

  path[nlocal][0] = atof(values[0]);
  path[nlocal][1] = atof(values[1]);
  path[nlocal][2] = atof(values[2]);

  norm[nlocal][0] = atof(values[3]);
  norm[nlocal][1] = atof(values[4]);
  norm[nlocal][2] = atof(values[5]);

  dnorm[nlocal][0] = atof(values[6]);
  dnorm[nlocal][1] = atof(values[7]);
  dnorm[nlocal][2] = atof(values[8]);

  return 9;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecPAFI::pack_data(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = ubuf(type[i]).d;
    buf[i][2] = x[i][0];
    buf[i][3] = x[i][1];
    buf[i][4] = x[i][2];

    buf[i][5] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
    buf[i][6] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][7] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;
  }
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecPAFI::pack_data_hybrid(int i, double *buf)
{
  buf[0] = path[i][0];
  buf[1] = path[i][1];
  buf[2] = path[i][2];

  buf[3] = norm[i][0];
  buf[4] = norm[i][1];
  buf[5] = norm[i][2];

  buf[6] = dnorm[i][0];
  buf[7] = dnorm[i][1];
  buf[8] = dnorm[i][2];
  return 9;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecPAFI::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT \
            " %d %-1.16e %-1.16e %-1.16e %d %d %d\n",
            (tagint) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,
            buf[i][2],buf[i][3],buf[i][4],
            (int) ubuf(buf[i][5]).i,(int) ubuf(buf[i][6]).i,
            (int) ubuf(buf[i][7]).i);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecPAFI::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp,"%-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e",
  buf[0],buf[1],buf[2],
  buf[3],buf[4],buf[5],
  buf[6],buf[7],buf[8]);
  return 9;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecPAFI::memory_usage()
{
  bigint bytes = 0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);

  if (atom->memcheck("path")) bytes += memory->usage(path,nmax,3);
  if (atom->memcheck("norm")) bytes += memory->usage(norm,nmax,3);
  if (atom->memcheck("dnorm")) bytes += memory->usage(dnorm,nmax,3);

  return bytes;
}
