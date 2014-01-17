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

/* ----------------------------------------------------------------------
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "atom_vec_ellipsoid.h"
#include "math_extra.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define DELTA 10000
#define DELTA_BONUS 10000

/* ---------------------------------------------------------------------- */

AtomVecEllipsoid::AtomVecEllipsoid(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;

  comm_x_only = comm_f_only = 0;
  size_forward = 7;
  size_reverse = 6;
  size_border = 14;
  size_velocity = 6;
  size_data_atom = 7;
  size_data_vel = 7;
  size_data_bonus = 8;
  xcol_data = 5;

  atom->ellipsoid_flag = 1;
  atom->rmass_flag = atom->angmom_flag = atom->torque_flag = 1;

  nlocal_bonus = nghost_bonus = nmax_bonus = 0;
  bonus = NULL;
}

/* ---------------------------------------------------------------------- */

AtomVecEllipsoid::~AtomVecEllipsoid()
{
  memory->sfree(bonus);
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecEllipsoid::grow(int n)
{
  if (n == 0) nmax += DELTA;
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  tag = memory->grow(atom->tag,nmax,"atom:tag");
  type = memory->grow(atom->type,nmax,"atom:type");
  mask = memory->grow(atom->mask,nmax,"atom:mask");
  image = memory->grow(atom->image,nmax,"atom:image");
  x = memory->grow(atom->x,nmax,3,"atom:x");
  v = memory->grow(atom->v,nmax,3,"atom:v");
  f = memory->grow(atom->f,nmax*comm->nthreads,3,"atom:f");

  rmass = memory->grow(atom->rmass,nmax,"atom:rmass");
  angmom = memory->grow(atom->angmom,nmax,3,"atom:angmom");
  torque = memory->grow(atom->torque,nmax*comm->nthreads,3,"atom:torque");
  ellipsoid = memory->grow(atom->ellipsoid,nmax,"atom:ellipsoid");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecEllipsoid::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;
  rmass = atom->rmass; angmom = atom->angmom; torque = atom->torque;
  ellipsoid = atom->ellipsoid;
}

/* ----------------------------------------------------------------------
   grow bonus data structure
------------------------------------------------------------------------- */

void AtomVecEllipsoid::grow_bonus()
{
  nmax_bonus += DELTA_BONUS;
  if (nmax_bonus < 0 || nmax_bonus > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  bonus = (Bonus *) memory->srealloc(bonus,nmax_bonus*sizeof(Bonus),
                                     "atom:bonus");
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecEllipsoid::copy(int i, int j, int delflag)
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

  rmass[j] = rmass[i];
  angmom[j][0] = angmom[i][0];
  angmom[j][1] = angmom[i][1];
  angmom[j][2] = angmom[i][2];

  // if deleting atom J via delflag and J has bonus data, then delete it

  if (delflag && ellipsoid[j] >= 0) {
    copy_bonus(nlocal_bonus-1,ellipsoid[j]);
    nlocal_bonus--;
  }

  // if atom I has bonus data, reset I's bonus.ilocal to loc J
  // do NOT do this if self-copy (I=J) since I's bonus data is already deleted

  if (ellipsoid[i] >= 0 && i != j) bonus[ellipsoid[i]].ilocal = j;
  ellipsoid[j] = ellipsoid[i];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ----------------------------------------------------------------------
   copy bonus data from I to J, effectively deleting the J entry
   also reset ellipsoid that points to I to now point to J
------------------------------------------------------------------------- */

void AtomVecEllipsoid::copy_bonus(int i, int j)
{
  ellipsoid[bonus[i].ilocal] = j;
  memcpy(&bonus[j],&bonus[i],sizeof(Bonus));
}

/* ----------------------------------------------------------------------
   clear ghost info in bonus data
   called before ghosts are recommunicated in comm and irregular
------------------------------------------------------------------------- */

void AtomVecEllipsoid::clear_bonus()
{
  nghost_bonus = 0;
}

/* ----------------------------------------------------------------------
   set shape values in bonus data for particle I
   oriented aligned with xyz axes
   this may create or delete entry in bonus data
------------------------------------------------------------------------- */

void AtomVecEllipsoid::set_shape(int i,
                                 double shapex, double shapey, double shapez)
{
  if (ellipsoid[i] < 0) {
    if (shapex == 0.0 && shapey == 0.0 && shapez == 0.0) return;
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    double *shape = bonus[nlocal_bonus].shape;
    double *quat = bonus[nlocal_bonus].quat;
    shape[0] = shapex;
    shape[1] = shapey;
    shape[2] = shapez;
    quat[0] = 1.0;
    quat[1] = 0.0;
    quat[2] = 0.0;
    quat[3] = 0.0;
    bonus[nlocal_bonus].ilocal = i;
    ellipsoid[i] = nlocal_bonus++;
  } else if (shapex == 0.0 && shapey == 0.0 && shapez == 0.0) {
    copy_bonus(nlocal_bonus-1,ellipsoid[i]);
    nlocal_bonus--;
    ellipsoid[i] = -1;
  } else {
    double *shape = bonus[ellipsoid[i]].shape;
    shape[0] = shapex;
    shape[1] = shapey;
    shape[2] = shapez;
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecEllipsoid::pack_comm(int n, int *list, double *buf,
                                int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;
  double *quat;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      if (ellipsoid[j] >= 0) {
        quat = bonus[ellipsoid[j]].quat;
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
      }
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
      if (ellipsoid[j] >= 0) {
        quat = bonus[ellipsoid[j]].quat;
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
      }
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecEllipsoid::pack_comm_vel(int n, int *list, double *buf,
                                    int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;
  double *quat;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      if (ellipsoid[j] >= 0) {
        quat = bonus[ellipsoid[j]].quat;
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
      }
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = angmom[j][0];
      buf[m++] = angmom[j][1];
      buf[m++] = angmom[j][2];
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
        if (ellipsoid[j] >= 0) {
          quat = bonus[ellipsoid[j]].quat;
          buf[m++] = quat[0];
          buf[m++] = quat[1];
          buf[m++] = quat[2];
          buf[m++] = quat[3];
        }
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
        buf[m++] = angmom[j][0];
        buf[m++] = angmom[j][1];
        buf[m++] = angmom[j][2];
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
        if (ellipsoid[j] >= 0) {
          quat = bonus[ellipsoid[j]].quat;
          buf[m++] = quat[0];
          buf[m++] = quat[1];
          buf[m++] = quat[2];
          buf[m++] = quat[3];
        }
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
        buf[m++] = angmom[j][0];
        buf[m++] = angmom[j][1];
        buf[m++] = angmom[j][2];
      }
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecEllipsoid::pack_comm_hybrid(int n, int *list, double *buf)
{
  int i,j,m;
  double *quat;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (ellipsoid[j] >= 0) {
      quat = bonus[ellipsoid[j]].quat;
      buf[m++] = quat[0];
      buf[m++] = quat[1];
      buf[m++] = quat[2];
      buf[m++] = quat[3];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoid::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;
  double *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    if (ellipsoid[i] >= 0) {
      quat = bonus[ellipsoid[i]].quat;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoid::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;
  double *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    if (ellipsoid[i] >= 0) {
      quat = bonus[ellipsoid[i]].quat;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
    }
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    angmom[i][0] = buf[m++];
    angmom[i][1] = buf[m++];
    angmom[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecEllipsoid::unpack_comm_hybrid(int n, int first, double *buf)
{
  int i,m,last;
  double *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (ellipsoid[i] >= 0) {
      quat = bonus[ellipsoid[i]].quat;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecEllipsoid::pack_reverse(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
    buf[m++] = torque[i][0];
    buf[m++] = torque[i][1];
    buf[m++] = torque[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecEllipsoid::pack_reverse_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = torque[i][0];
    buf[m++] = torque[i][1];
    buf[m++] = torque[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoid::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
    torque[j][0] += buf[m++];
    torque[j][1] += buf[m++];
    torque[j][2] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecEllipsoid::unpack_reverse_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    torque[j][0] += buf[m++];
    torque[j][1] += buf[m++];
    torque[j][2] += buf[m++];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecEllipsoid::pack_border(int n, int *list, double *buf,
                                  int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;
  double *shape,*quat;

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
      if (ellipsoid[j] < 0) buf[m++] = ubuf(0).d;
      else {
        buf[m++] = ubuf(1).d;
        shape = bonus[ellipsoid[j]].shape;
        quat = bonus[ellipsoid[j]].quat;
        buf[m++] = shape[0];
        buf[m++] = shape[1];
        buf[m++] = shape[2];
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
      }
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
      if (ellipsoid[j] < 0) buf[m++] = ubuf(0).d;
      else {
        buf[m++] = ubuf(1).d;
        shape = bonus[ellipsoid[j]].shape;
        quat = bonus[ellipsoid[j]].quat;
        buf[m++] = shape[0];
        buf[m++] = shape[1];
        buf[m++] = shape[2];
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecEllipsoid::pack_border_vel(int n, int *list, double *buf,
                                      int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;
  double *shape,*quat;

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
      if (ellipsoid[j] < 0) buf[m++] = ubuf(0).d;
      else {
        buf[m++] = ubuf(1).d;
        shape = bonus[ellipsoid[j]].shape;
        quat = bonus[ellipsoid[j]].quat;
        buf[m++] = shape[0];
        buf[m++] = shape[1];
        buf[m++] = shape[2];
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
      }
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = angmom[j][0];
      buf[m++] = angmom[j][1];
      buf[m++] = angmom[j][2];
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
        if (ellipsoid[j] < 0) buf[m++] = ubuf(0).d;
        else {
          buf[m++] = ubuf(1).d;
          shape = bonus[ellipsoid[j]].shape;
          quat  = bonus[ellipsoid[j]].quat;
          buf[m++] = shape[0];
          buf[m++] = shape[1];
          buf[m++] = shape[2];
          buf[m++] = quat[0];
          buf[m++] = quat[1];
          buf[m++] = quat[2];
          buf[m++] = quat[3];
        }
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
        buf[m++] = angmom[j][0];
        buf[m++] = angmom[j][1];
        buf[m++] = angmom[j][2];
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
        if (ellipsoid[j] < 0) buf[m++] = ubuf(0).d;
        else {
          buf[m++] = ubuf(1).d;
          shape = bonus[ellipsoid[j]].shape;
          quat = bonus[ellipsoid[j]].quat;
          buf[m++] = shape[0];
          buf[m++] = shape[1];
          buf[m++] = shape[2];
          buf[m++] = quat[0];
          buf[m++] = quat[1];
          buf[m++] = quat[2];
          buf[m++] = quat[3];
        }
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
        buf[m++] = angmom[j][0];
        buf[m++] = angmom[j][1];
        buf[m++] = angmom[j][2];
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecEllipsoid::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;
  double *shape,*quat;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (ellipsoid[j] < 0) buf[m++] = ubuf(0).d;
    else {
      buf[m++] = ubuf(1).d;
      shape = bonus[ellipsoid[j]].shape;
      quat = bonus[ellipsoid[j]].quat;
      buf[m++] = shape[0];
      buf[m++] = shape[1];
      buf[m++] = shape[2];
      buf[m++] = quat[0];
      buf[m++] = quat[1];
      buf[m++] = quat[2];
      buf[m++] = quat[3];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoid::unpack_border(int n, int first, double *buf)
{
  int i,j,m,last;
  double *shape,*quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (int) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    ellipsoid[i] = (int) ubuf(buf[m++]).i;
    if (ellipsoid[i] == 0) ellipsoid[i] = -1;
    else {
      j = nlocal_bonus + nghost_bonus;
      if (j == nmax_bonus) grow_bonus();
      shape = bonus[j].shape;
      quat = bonus[j].quat;
      shape[0] = buf[m++];
      shape[1] = buf[m++];
      shape[2] = buf[m++];
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
      bonus[j].ilocal = i;
      ellipsoid[i] = j;
      nghost_bonus++;
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoid::unpack_border_vel(int n, int first, double *buf)
{
  int i,j,m,last;
  double *shape,*quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (int) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    ellipsoid[i] = (int) ubuf(buf[m++]).i;
    if (ellipsoid[i] == 0) ellipsoid[i] = -1;
    else {
      j = nlocal_bonus + nghost_bonus;
      if (j == nmax_bonus) grow_bonus();
      shape = bonus[j].shape;
      quat = bonus[j].quat;
      shape[0] = buf[m++];
      shape[1] = buf[m++];
      shape[2] = buf[m++];
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
      bonus[j].ilocal = i;
      ellipsoid[i] = j;
      nghost_bonus++;
    }
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    angmom[i][0] = buf[m++];
    angmom[i][1] = buf[m++];
    angmom[i][2] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

int AtomVecEllipsoid::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,j,m,last;
  double *shape,*quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    ellipsoid[i] = (int) ubuf(buf[m++]).i;
    if (ellipsoid[i] == 0) ellipsoid[i] = -1;
    else {
      j = nlocal_bonus + nghost_bonus;
      if (j == nmax_bonus) grow_bonus();
      shape = bonus[j].shape;
      quat = bonus[j].quat;
      shape[0] = buf[m++];
      shape[1] = buf[m++];
      shape[2] = buf[m++];
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
      bonus[j].ilocal = i;
      ellipsoid[i] = j;
      nghost_bonus++;
    }
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecEllipsoid::pack_exchange(int i, double *buf)
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

  buf[m++] = rmass[i];
  buf[m++] = angmom[i][0];
  buf[m++] = angmom[i][1];
  buf[m++] = angmom[i][2];

  if (ellipsoid[i] < 0) buf[m++] = ubuf(0).d;
  else {
    buf[m++] = ubuf(1).d;
    int j = ellipsoid[i];
    double *shape = bonus[j].shape;
    double *quat = bonus[j].quat;
    buf[m++] = shape[0];
    buf[m++] = shape[1];
    buf[m++] = shape[2];
    buf[m++] = quat[0];
    buf[m++] = quat[1];
    buf[m++] = quat[2];
    buf[m++] = quat[3];
  }

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecEllipsoid::unpack_exchange(double *buf)
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
  tag[nlocal] = (int) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;

  rmass[nlocal] = buf[m++];
  angmom[nlocal][0] = buf[m++];
  angmom[nlocal][1] = buf[m++];
  angmom[nlocal][2] = buf[m++];

  ellipsoid[nlocal] = (int) ubuf(buf[m++]).i;
  if (ellipsoid[nlocal] == 0) ellipsoid[nlocal] = -1;
  else {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    double *shape = bonus[nlocal_bonus].shape;
    double *quat = bonus[nlocal_bonus].quat;
    shape[0] = buf[m++];
    shape[1] = buf[m++];
    shape[2] = buf[m++];
    quat[0] = buf[m++];
    quat[1] = buf[m++];
    quat[2] = buf[m++];
    quat[3] = buf[m++];
    bonus[nlocal_bonus].ilocal = nlocal;
    ellipsoid[nlocal] = nlocal_bonus++;
  }

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

int AtomVecEllipsoid::size_restart()
{
  int i;

  int n = 0;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++)
    if (ellipsoid[i] >= 0) n += 23;
    else n += 16;

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
        n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including bonus data
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVecEllipsoid::pack_restart(int i, double *buf)
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

  buf[m++] = rmass[i];
  buf[m++] = angmom[i][0];
  buf[m++] = angmom[i][1];
  buf[m++] = angmom[i][2];

  if (ellipsoid[i] < 0) buf[m++] = ubuf(0).d;
  else {
    buf[m++] = ubuf(1).d;
    int j = ellipsoid[i];
    buf[m++] = bonus[j].shape[0];
    buf[m++] = bonus[j].shape[1];
    buf[m++] = bonus[j].shape[2];
    buf[m++] = bonus[j].quat[0];
    buf[m++] = bonus[j].quat[1];
    buf[m++] = bonus[j].quat[2];
    buf[m++] = bonus[j].quat[3];
  }

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including bonus data
------------------------------------------------------------------------- */

int AtomVecEllipsoid::unpack_restart(double *buf)
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
  tag[nlocal] = (int) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  rmass[nlocal] = buf[m++];
  angmom[nlocal][0] = buf[m++];
  angmom[nlocal][1] = buf[m++];
  angmom[nlocal][2] = buf[m++];

  ellipsoid[nlocal] = (int) ubuf(buf[m++]).i;
  if (ellipsoid[nlocal] == 0) ellipsoid[nlocal] = -1;
  else {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    double *shape = bonus[nlocal_bonus].shape;
    double *quat = bonus[nlocal_bonus].quat;
    shape[0] = buf[m++];
    shape[1] = buf[m++];
    shape[2] = buf[m++];
    quat[0] = buf[m++];
    quat[1] = buf[m++];
    quat[2] = buf[m++];
    quat[3] = buf[m++];
    bonus[nlocal_bonus].ilocal = nlocal;
    ellipsoid[nlocal] = nlocal_bonus++;
  }

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

void AtomVecEllipsoid::create_atom(int itype, double *coord)
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

  rmass[nlocal] = 1.0;
  angmom[nlocal][0] = 0.0;
  angmom[nlocal][1] = 0.0;
  angmom[nlocal][2] = 0.0;
  ellipsoid[nlocal] = -1;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecEllipsoid::data_atom(double *coord, imageint imagetmp, 
                                 char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = atoi(values[0]);
  if (tag[nlocal] <= 0)
    error->one(FLERR,"Invalid atom ID in Atoms section of data file");

  type[nlocal] = atoi(values[1]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  ellipsoid[nlocal] = atoi(values[2]);
  if (ellipsoid[nlocal] == 0) ellipsoid[nlocal] = -1;
  else if (ellipsoid[nlocal] == 1) ellipsoid[nlocal] = 0;
  else error->one(FLERR,"Invalid atom type in Atoms section of data file");

  rmass[nlocal] = atof(values[3]);
  if (rmass[nlocal] <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  angmom[nlocal][0] = 0.0;
  angmom[nlocal][1] = 0.0;
  angmom[nlocal][2] = 0.0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecEllipsoid::data_atom_hybrid(int nlocal, char **values)
{
  ellipsoid[nlocal] = atoi(values[0]);
  if (ellipsoid[nlocal] == 0) ellipsoid[nlocal] = -1;
  else if (ellipsoid[nlocal] == 1) ellipsoid[nlocal] = 0;
  else error->one(FLERR,"Invalid atom type in Atoms section of data file");

  rmass[nlocal] = atof(values[1]);
  if (rmass[nlocal] <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  return 2;
}

/* ----------------------------------------------------------------------
   unpack one line from Ellipsoids section of data file
------------------------------------------------------------------------- */

void AtomVecEllipsoid::data_atom_bonus(int m, char **values)
{
  if (ellipsoid[m])
    error->one(FLERR,"Assigning ellipsoid parameters to non-ellipsoid atom");

  if (nlocal_bonus == nmax_bonus) grow_bonus();

  double *shape = bonus[nlocal_bonus].shape;
  shape[0] = 0.5 * atof(values[0]);
  shape[1] = 0.5 * atof(values[1]);
  shape[2] = 0.5 * atof(values[2]);
  if (shape[0] <= 0.0 || shape[1] <= 0.0 || shape[2] <= 0.0)
    error->one(FLERR,"Invalid shape in Ellipsoids section of data file");

  double *quat = bonus[nlocal_bonus].quat;
  quat[0] = atof(values[3]);
  quat[1] = atof(values[4]);
  quat[2] = atof(values[5]);
  quat[3] = atof(values[6]);
  MathExtra::qnormalize(quat);

  // reset ellipsoid mass
  // previously stored density in rmass

  rmass[m] *= 4.0*MY_PI/3.0 * shape[0]*shape[1]*shape[2];

  bonus[nlocal_bonus].ilocal = m;
  ellipsoid[m] = nlocal_bonus++;
}

/* ----------------------------------------------------------------------
   unpack one line from Velocities section of data file
------------------------------------------------------------------------- */

void AtomVecEllipsoid::data_vel(int m, char **values)
{
  v[m][0] = atof(values[0]);
  v[m][1] = atof(values[1]);
  v[m][2] = atof(values[2]);
  angmom[m][0] = atof(values[3]);
  angmom[m][1] = atof(values[4]);
  angmom[m][2] = atof(values[5]);
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Velocities section of data file
------------------------------------------------------------------------- */

int AtomVecEllipsoid::data_vel_hybrid(int m, char **values)
{
  angmom[m][0] = atof(values[0]);
  angmom[m][1] = atof(values[1]);
  angmom[m][2] = atof(values[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecEllipsoid::pack_data(double **buf)
{
  double *shape;

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = ubuf(type[i]).d;
    if (ellipsoid[i] < 0) buf[i][2] = ubuf(0).d;
    else buf[i][2] = ubuf(1).d;
    if (ellipsoid[i] < 0) buf[i][3] = rmass[i];
    else {
      shape = bonus[ellipsoid[i]].shape;
      buf[i][3] = rmass[i] / (4.0*MY_PI/3.0 * shape[0]*shape[1]*shape[2]);
    }
    buf[i][4] = x[i][0];
    buf[i][5] = x[i][1];
    buf[i][6] = x[i][2];
    buf[i][7] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
    buf[i][8] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][9] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;
  }
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecEllipsoid::pack_data_hybrid(int i, double *buf)
{
  if (ellipsoid[i] < 0) buf[0] = ubuf(0).d;
  else buf[0] = ubuf(1).d;
  if (ellipsoid[i] < 0) buf[1] = rmass[i];
  else {
    double *shape = bonus[ellipsoid[i]].shape;
    buf[1] = rmass[i] / (4.0*MY_PI/3.0 * shape[0]*shape[1]*shape[2]);
  }
  return 2;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecEllipsoid::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,"%d %d %d %-1.16e %-1.16e %-1.16e %-1.16e %d %d %d\n",
            (int) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,
            (int) ubuf(buf[i][2]).i,
            buf[i][3],buf[i][4],buf[i][5],buf[i][6],
            (int) ubuf(buf[i][7]).i,(int) ubuf(buf[i][8]).i,
            (int) ubuf(buf[i][9]).i);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecEllipsoid::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %d %-1.16e",(int) ubuf(buf[0]).i,buf[1]);
  return 2;
}

/* ----------------------------------------------------------------------
   pack velocity info for data file
------------------------------------------------------------------------- */

void AtomVecEllipsoid::pack_vel(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = v[i][0];
    buf[i][2] = v[i][1];
    buf[i][3] = v[i][2];
    buf[i][4] = angmom[i][0];
    buf[i][5] = angmom[i][1];
    buf[i][6] = angmom[i][2];
  }
}

/* ----------------------------------------------------------------------
   pack hybrid velocity info for data file
------------------------------------------------------------------------- */

int AtomVecEllipsoid::pack_vel_hybrid(int i, double *buf)
{
  buf[0] = angmom[i][0];
  buf[1] = angmom[i][1];
  buf[2] = angmom[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
   write velocity info to data file
------------------------------------------------------------------------- */

void AtomVecEllipsoid::write_vel(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,"%d %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e\n",
            (int) ubuf(buf[i][0]).i,buf[i][1],buf[i][2],buf[i][3],
            buf[i][4],buf[i][5],buf[i][6]);
}

/* ----------------------------------------------------------------------
   write hybrid velocity info to data file
------------------------------------------------------------------------- */

int AtomVecEllipsoid::write_vel_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e %-1.16e %-1.16e",buf[0],buf[1],buf[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecEllipsoid::memory_usage()
{
  bigint bytes = 0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);

  if (atom->memcheck("rmass")) bytes += memory->usage(rmass,nmax);
  if (atom->memcheck("angmom")) bytes += memory->usage(angmom,nmax,3);
  if (atom->memcheck("torque")) 
    bytes += memory->usage(torque,nmax*comm->nthreads,3);
  if (atom->memcheck("ellipsoid")) bytes += memory->usage(ellipsoid,nmax);

  bytes += nmax_bonus*sizeof(Bonus);

  return bytes;
}
