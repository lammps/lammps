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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "atom_vec_tri.h"
#include "math_extra.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "force.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 10000
#define DELTA_BONUS 10000
#define EPSILON 0.001

/* ---------------------------------------------------------------------- */

AtomVecTri::AtomVecTri(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;

  comm_x_only = comm_f_only = 0;
  size_forward = 7;
  size_reverse = 6;
  size_border = 24;
  size_velocity = 6;
  size_data_atom = 8;
  size_data_vel = 7;
  size_data_bonus = 10;
  xcol_data = 6;

  atom->tri_flag = 1;
  atom->molecule_flag = atom->rmass_flag = 1;
  atom->angmom_flag = atom->torque_flag = 1;

  nlocal_bonus = nghost_bonus = nmax_bonus = 0;
  bonus = NULL;
}

/* ---------------------------------------------------------------------- */

AtomVecTri::~AtomVecTri()
{
  memory->sfree(bonus);
}

/* ---------------------------------------------------------------------- */

void AtomVecTri::init()
{
  AtomVec::init();

  if (domain->dimension != 3)
    error->all(FLERR,"Atom_style tri can only be used in 3d simulations");
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecTri::grow(int n)
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

  molecule = memory->grow(atom->molecule,nmax,"atom:molecule");
  rmass = memory->grow(atom->rmass,nmax,"atom:rmass");
  angmom = memory->grow(atom->angmom,nmax,3,"atom:angmom");
  torque = memory->grow(atom->torque,nmax*comm->nthreads,3,"atom:torque");
  tri = memory->grow(atom->tri,nmax,"atom:tri");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecTri::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;
  molecule = atom->molecule; rmass = atom->rmass;
  angmom = atom->angmom; torque = atom->torque;
  tri = atom->tri;
}

/* ----------------------------------------------------------------------
   grow bonus data structure
------------------------------------------------------------------------- */

void AtomVecTri::grow_bonus()
{
  nmax_bonus += DELTA_BONUS;
  if (nmax_bonus < 0 || nmax_bonus > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  bonus = (Bonus *) memory->srealloc(bonus,nmax_bonus*sizeof(Bonus),
                                     "atom:bonus");
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
   if delflag and atom J has bonus data, then delete it
------------------------------------------------------------------------- */

void AtomVecTri::copy(int i, int j, int delflag)
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

  molecule[j] = molecule[i];
  rmass[j] = rmass[i];
  angmom[j][0] = angmom[i][0];
  angmom[j][1] = angmom[i][1];
  angmom[j][2] = angmom[i][2];

  // if deleting atom J via delflag and J has bonus data, then delete it

  if (delflag && tri[j] >= 0) {
    copy_bonus(nlocal_bonus-1,tri[j]);
    nlocal_bonus--;
  }

  // if atom I has bonus data, reset I's bonus.ilocal to loc J
  // do NOT do this if self-copy (I=J) since I's bonus data is already deleted

  if (tri[i] >= 0 && i != j) bonus[tri[i]].ilocal = j;
  tri[j] = tri[i];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ----------------------------------------------------------------------
   copy bonus data from I to J, effectively deleting the J entry
   also reset tri that points to I to now point to J
------------------------------------------------------------------------- */

void AtomVecTri::copy_bonus(int i, int j)
{
  tri[bonus[i].ilocal] = j;
  memcpy(&bonus[j],&bonus[i],sizeof(Bonus));
}

/* ----------------------------------------------------------------------
   clear ghost info in bonus data
   called before ghosts are recommunicated in comm and irregular
------------------------------------------------------------------------- */

void AtomVecTri::clear_bonus()
{
  nghost_bonus = 0;
}

/* ----------------------------------------------------------------------
   set equilateral tri of size in bonus data for particle I
   oriented symmetrically in xy plane
   this may create or delete entry in bonus data
------------------------------------------------------------------------- */

void AtomVecTri::set_equilateral(int i, double size)
{
  if (tri[i] < 0) {
    if (size == 0.0) return;
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    double *quat = bonus[nlocal_bonus].quat;
    double *c1 = bonus[nlocal_bonus].c1;
    double *c2 = bonus[nlocal_bonus].c2;
    double *c3 = bonus[nlocal_bonus].c3;
    double *inertia = bonus[nlocal_bonus].inertia;
    quat[0] = 1.0;
    quat[1] = 0.0;
    quat[2] = 0.0;
    quat[3] = 0.0;
    c1[0] = -size/2.0;
    c1[1] = -sqrt(3.0)/2.0 * size / 3.0;
    c1[2] = 0.0;
    c2[0] = size/2.0;
    c2[1] = -sqrt(3.0)/2.0 * size / 3.0;
    c2[2] = 0.0;
    c3[0] = 0.0;
    c3[1] = sqrt(3.0)/2.0 * size * 2.0/3.0;
    c3[2] = 0.0;
    inertia[0] = sqrt(3.0)/96.0 * size*size*size*size;
    inertia[1] = sqrt(3.0)/96.0 * size*size*size*size;
    inertia[2] = sqrt(3.0)/48.0 * size*size*size*size;
    bonus[nlocal_bonus].ilocal = i;
    tri[i] = nlocal_bonus++;
  } else if (size == 0.0) {
    copy_bonus(nlocal_bonus-1,tri[i]);
    nlocal_bonus--;
    tri[i] = -1;
  } else {
    double *c1 = bonus[tri[i]].c1;
    double *c2 = bonus[tri[i]].c2;
    double *c3 = bonus[tri[i]].c3;
    double *inertia = bonus[tri[i]].inertia;
    c1[0] = -size/2.0;
    c1[1] = -sqrt(3.0)/2.0 * size / 3.0;
    c1[2] = 0.0;
    c2[0] = size/2.0;
    c2[1] = -sqrt(3.0)/2.0 * size / 3.0;
    c2[2] = 0.0;
    c3[0] = 0.0;
    c3[1] = sqrt(3.0)/2.0 * size * 2.0/3.0;
    c3[2] = 0.0;
    inertia[0] = sqrt(3.0)/96.0 * size*size*size*size;
    inertia[1] = sqrt(3.0)/96.0 * size*size*size*size;
    inertia[2] = sqrt(3.0)/48.0 * size*size*size*size;
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecTri::pack_comm(int n, int *list, double *buf,
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
      if (tri[j] >= 0) {
        quat = bonus[tri[j]].quat;
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
      if (tri[j] >= 0) {
        quat = bonus[tri[j]].quat;
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

int AtomVecTri::pack_comm_vel(int n, int *list, double *buf,
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
      if (tri[j] >= 0) {
        quat = bonus[tri[j]].quat;
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
        if (tri[j] >= 0) {
          quat = bonus[tri[j]].quat;
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
        if (tri[j] >= 0) {
          quat = bonus[tri[j]].quat;
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

int AtomVecTri::pack_comm_hybrid(int n, int *list, double *buf)
{
  int i,j,m;
  double *quat;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (tri[j] >= 0) {
      quat = bonus[tri[j]].quat;
      buf[m++] = quat[0];
      buf[m++] = quat[1];
      buf[m++] = quat[2];
      buf[m++] = quat[3];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecTri::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;
  double *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    if (tri[i] >= 0) {
      quat = bonus[tri[i]].quat;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecTri::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;
  double *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    if (tri[i] >= 0) {
      quat = bonus[tri[i]].quat;
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

int AtomVecTri::unpack_comm_hybrid(int n, int first, double *buf)
{
  int i,m,last;
  double *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    if (tri[i] >= 0) {
      quat = bonus[tri[i]].quat;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
    }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecTri::pack_reverse(int n, int first, double *buf)
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

int AtomVecTri::pack_reverse_hybrid(int n, int first, double *buf)
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

void AtomVecTri::unpack_reverse(int n, int *list, double *buf)
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

int AtomVecTri::unpack_reverse_hybrid(int n, int *list, double *buf)
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

int AtomVecTri::pack_border(int n, int *list, double *buf,
                            int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;
  double *quat,*c1,*c2,*c3,*inertia;

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
      if (tri[j] < 0) buf[m++] = ubuf(0).d;
      else {
        buf[m++] = ubuf(1).d;
        quat = bonus[tri[j]].quat;
        c1 = bonus[tri[j]].c1;
        c2 = bonus[tri[j]].c2;
        c3 = bonus[tri[j]].c3;
        inertia = bonus[tri[j]].inertia;
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
        buf[m++] = c1[0];
        buf[m++] = c1[1];
        buf[m++] = c1[2];
        buf[m++] = c2[0];
        buf[m++] = c2[1];
        buf[m++] = c2[2];
        buf[m++] = c3[0];
        buf[m++] = c3[1];
        buf[m++] = c3[2];
        buf[m++] = inertia[0];
        buf[m++] = inertia[1];
        buf[m++] = inertia[2];
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
      buf[m++] = ubuf(molecule[j]).d;
      if (tri[j] < 0) buf[m++] = ubuf(0).d;
      else {
        buf[m++] = ubuf(1).d;
        quat = bonus[tri[j]].quat;
        c1 = bonus[tri[j]].c1;
        c2 = bonus[tri[j]].c2;
        c3 = bonus[tri[j]].c3;
        inertia = bonus[tri[j]].inertia;
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
        buf[m++] = c1[0];
        buf[m++] = c1[1];
        buf[m++] = c1[2];
        buf[m++] = c2[0];
        buf[m++] = c2[1];
        buf[m++] = c2[2];
        buf[m++] = c3[0];
        buf[m++] = c3[1];
        buf[m++] = c3[2];
        buf[m++] = inertia[0];
        buf[m++] = inertia[1];
        buf[m++] = inertia[2];
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecTri::pack_border_vel(int n, int *list, double *buf,
                                int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;
  double *quat,*c1,*c2,*c3,*inertia;

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
      if (tri[j] < 0) buf[m++] = ubuf(0).d;
      else {
        buf[m++] = ubuf(1).d;
        quat = bonus[tri[j]].quat;
        c1 = bonus[tri[j]].c1;
        c2 = bonus[tri[j]].c2;
        c3 = bonus[tri[j]].c3;
        inertia = bonus[tri[j]].inertia;
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
        buf[m++] = c1[0];
        buf[m++] = c1[1];
        buf[m++] = c1[2];
        buf[m++] = c2[0];
        buf[m++] = c2[1];
        buf[m++] = c2[2];
        buf[m++] = c3[0];
        buf[m++] = c3[1];
        buf[m++] = c3[2];
        buf[m++] = inertia[0];
        buf[m++] = inertia[1];
        buf[m++] = inertia[2];
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
        buf[m++] = ubuf(molecule[j]).d;
        if (tri[j] < 0) buf[m++] = ubuf(0).d;
        else {
          buf[m++] = ubuf(1).d;
          quat = bonus[tri[j]].quat;
          c1 = bonus[tri[j]].c1;
          c2 = bonus[tri[j]].c2;
          c3 = bonus[tri[j]].c3;
          inertia = bonus[tri[j]].inertia;
          buf[m++] = quat[0];
          buf[m++] = quat[1];
          buf[m++] = quat[2];
          buf[m++] = quat[3];
          buf[m++] = c1[0];
          buf[m++] = c1[1];
          buf[m++] = c1[2];
          buf[m++] = c2[0];
          buf[m++] = c2[1];
          buf[m++] = c2[2];
          buf[m++] = c3[0];
          buf[m++] = c3[1];
          buf[m++] = c3[2];
          buf[m++] = inertia[0];
          buf[m++] = inertia[1];
          buf[m++] = inertia[2];
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
        buf[m++] = ubuf(molecule[j]).d;
        if (tri[j] < 0) buf[m++] = ubuf(0).d;
        else {
          buf[m++] = ubuf(1).d;
          quat = bonus[tri[j]].quat;
          c1 = bonus[tri[j]].c1;
          c2 = bonus[tri[j]].c2;
          c3 = bonus[tri[j]].c3;
          inertia = bonus[tri[j]].inertia;
          buf[m++] = quat[0];
          buf[m++] = quat[1];
          buf[m++] = quat[2];
          buf[m++] = quat[3];
          buf[m++] = c1[0];
          buf[m++] = c1[1];
          buf[m++] = c1[2];
          buf[m++] = c2[0];
          buf[m++] = c2[1];
          buf[m++] = c2[2];
          buf[m++] = c3[0];
          buf[m++] = c3[1];
          buf[m++] = c3[2];
          buf[m++] = inertia[0];
          buf[m++] = inertia[1];
          buf[m++] = inertia[2];
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

int AtomVecTri::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;
  double *quat,*c1,*c2,*c3,*inertia;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = ubuf(molecule[j]).d;
    if (tri[j] < 0) buf[m++] = ubuf(0).d;
    else {
      buf[m++] = ubuf(1).d;
      quat = bonus[tri[j]].quat;
      c1 = bonus[tri[j]].c1;
      c2 = bonus[tri[j]].c2;
      c3 = bonus[tri[j]].c3;
      inertia = bonus[tri[j]].inertia;
      buf[m++] = quat[0];
      buf[m++] = quat[1];
      buf[m++] = quat[2];
      buf[m++] = quat[3];
      buf[m++] = c1[0];
      buf[m++] = c1[1];
      buf[m++] = c1[2];
      buf[m++] = c2[0];
      buf[m++] = c2[1];
      buf[m++] = c2[2];
      buf[m++] = c3[0];
      buf[m++] = c3[1];
      buf[m++] = c3[2];
      buf[m++] = inertia[0];
      buf[m++] = inertia[1];
      buf[m++] = inertia[2];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecTri::unpack_border(int n, int first, double *buf)
{
  int i,j,m,last;
  double *quat,*c1,*c2,*c3,*inertia;

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
    molecule[i] = (int) ubuf(buf[m++]).i;
    tri[i] = (int) ubuf(buf[m++]).i;
    if (tri[i] == 0) tri[i] = -1;
    else {
      j = nlocal_bonus + nghost_bonus;
      if (j == nmax_bonus) grow_bonus();
      quat = bonus[j].quat;
      c1 = bonus[j].c1;
      c2 = bonus[j].c2;
      c3 = bonus[j].c3;
      inertia = bonus[j].inertia;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
      c1[0] = buf[m++];
      c1[1] = buf[m++];
      c1[2] = buf[m++];
      c2[0] = buf[m++];
      c2[1] = buf[m++];
      c2[2] = buf[m++];
      c3[0] = buf[m++];
      c3[1] = buf[m++];
      c3[2] = buf[m++];
      inertia[0] = buf[m++];
      inertia[1] = buf[m++];
      inertia[2] = buf[m++];
      bonus[j].ilocal = i;
      tri[i] = j;
      nghost_bonus++;
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecTri::unpack_border_vel(int n, int first, double *buf)
{
  int i,j,m,last;
  double *quat,*c1,*c2,*c3,*inertia;

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
    molecule[i] = (int) ubuf(buf[m++]).i;
    tri[i] = (int) ubuf(buf[m++]).i;
    if (tri[i] == 0) tri[i] = -1;
    else {
      j = nlocal_bonus + nghost_bonus;
      if (j == nmax_bonus) grow_bonus();
      quat = bonus[j].quat;
      c1 = bonus[j].c1;
      c2 = bonus[j].c2;
      c3 = bonus[j].c3;
      inertia = bonus[j].inertia;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
      c1[0] = buf[m++];
      c1[1] = buf[m++];
      c1[2] = buf[m++];
      c2[0] = buf[m++];
      c2[1] = buf[m++];
      c2[2] = buf[m++];
      c3[0] = buf[m++];
      c3[1] = buf[m++];
      c3[2] = buf[m++];
      inertia[0] = buf[m++];
      inertia[1] = buf[m++];
      inertia[2] = buf[m++];
      bonus[j].ilocal = i;
      tri[i] = j;
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

int AtomVecTri::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,j,m,last;
  double *quat,*c1,*c2,*c3,*inertia;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    molecule[i] = (int) ubuf(buf[m++]).i;
    tri[i] = (int) ubuf(buf[m++]).i;
    if (tri[i] == 0) tri[i] = -1;
    else {
      j = nlocal_bonus + nghost_bonus;
      if (j == nmax_bonus) grow_bonus();
      quat = bonus[j].quat;
      c1 = bonus[j].c1;
      c2 = bonus[j].c2;
      c3 = bonus[j].c3;
      inertia = bonus[j].inertia;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
      c1[0] = buf[m++];
      c1[1] = buf[m++];
      c1[2] = buf[m++];
      c2[0] = buf[m++];
      c2[1] = buf[m++];
      c2[2] = buf[m++];
      c3[0] = buf[m++];
      c3[1] = buf[m++];
      c3[2] = buf[m++];
      inertia[0] = buf[m++];
      inertia[1] = buf[m++];
      inertia[2] = buf[m++];
      bonus[j].ilocal = i;
      tri[i] = j;
      nghost_bonus++;
    }
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecTri::pack_exchange(int i, double *buf)
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

  buf[m++] = ubuf(molecule[i]).d;
  buf[m++] = rmass[i];
  buf[m++] = angmom[i][0];
  buf[m++] = angmom[i][1];
  buf[m++] = angmom[i][2];

  if (tri[i] < 0) buf[m++] = ubuf(0).d;
  else {
    buf[m++] = ubuf(1).d;
    int j = tri[i];
    double *quat = bonus[j].quat;
    double *c1 = bonus[j].c1;
    double *c2 = bonus[j].c2;
    double *c3 = bonus[j].c3;
    double *inertia = bonus[j].inertia;
    buf[m++] = quat[0];
    buf[m++] = quat[1];
    buf[m++] = quat[2];
    buf[m++] = quat[3];
    buf[m++] = c1[0];
    buf[m++] = c1[1];
    buf[m++] = c1[2];
    buf[m++] = c2[0];
    buf[m++] = c2[1];
    buf[m++] = c2[2];
    buf[m++] = c3[0];
    buf[m++] = c3[1];
    buf[m++] = c3[2];
    buf[m++] = inertia[0];
    buf[m++] = inertia[1];
    buf[m++] = inertia[2];
  }

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecTri::unpack_exchange(double *buf)
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

  molecule[nlocal] = (int) ubuf(buf[m++]).i;
  rmass[nlocal] = buf[m++];
  angmom[nlocal][0] = buf[m++];
  angmom[nlocal][1] = buf[m++];
  angmom[nlocal][2] = buf[m++];

  tri[nlocal] = (int) ubuf(buf[m++]).i;
  if (tri[nlocal] == 0) tri[nlocal] = -1;
  else {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    double *quat = bonus[nlocal_bonus].quat;
    double *c1 = bonus[nlocal_bonus].c1;
    double *c2 = bonus[nlocal_bonus].c2;
    double *c3 = bonus[nlocal_bonus].c3;
    double *inertia = bonus[nlocal_bonus].inertia;
    quat[0] = buf[m++];
    quat[1] = buf[m++];
    quat[2] = buf[m++];
    quat[3] = buf[m++];
    c1[0] = buf[m++];
    c1[1] = buf[m++];
    c1[2] = buf[m++];
    c2[0] = buf[m++];
    c2[1] = buf[m++];
    c2[2] = buf[m++];
    c3[0] = buf[m++];
    c3[1] = buf[m++];
    c3[2] = buf[m++];
    inertia[0] = buf[m++];
    inertia[1] = buf[m++];
    inertia[2] = buf[m++];
    bonus[nlocal_bonus].ilocal = nlocal;
    tri[nlocal] = nlocal_bonus++;
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

int AtomVecTri::size_restart()
{
  int i;

  int n = 0;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++)
    if (tri[i] >= 0) n += 33;
    else n += 17;

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

int AtomVecTri::pack_restart(int i, double *buf)
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

  buf[m++] = ubuf(molecule[i]).d;
  buf[m++] = rmass[i];
  buf[m++] = angmom[i][0];
  buf[m++] = angmom[i][1];
  buf[m++] = angmom[i][2];

  if (tri[i] < 0) buf[m++] = ubuf(0).d;
  else {
    buf[m++] = ubuf(1).d;
    int j = tri[i];
    double *quat = bonus[j].quat;
    double *c1 = bonus[j].c1;
    double *c2 = bonus[j].c2;
    double *c3 = bonus[j].c3;
    double *inertia = bonus[j].inertia;
    buf[m++] = quat[0];
    buf[m++] = quat[1];
    buf[m++] = quat[2];
    buf[m++] = quat[3];
    buf[m++] = c1[0];
    buf[m++] = c1[1];
    buf[m++] = c1[2];
    buf[m++] = c2[0];
    buf[m++] = c2[1];
    buf[m++] = c2[2];
    buf[m++] = c3[0];
    buf[m++] = c3[1];
    buf[m++] = c3[2];
    buf[m++] = inertia[0];
    buf[m++] = inertia[1];
    buf[m++] = inertia[2];
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

int AtomVecTri::unpack_restart(double *buf)
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

  molecule[nlocal] = (int) ubuf(buf[m++]).i;
  rmass[nlocal] = buf[m++];
  angmom[nlocal][0] = buf[m++];
  angmom[nlocal][1] = buf[m++];
  angmom[nlocal][2] = buf[m++];

  tri[nlocal] = (int) ubuf(buf[m++]).i;
  if (tri[nlocal] == 0) tri[nlocal] = -1;
  else {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    double *quat = bonus[nlocal_bonus].quat;
    double *c1 = bonus[nlocal_bonus].c1;
    double *c2 = bonus[nlocal_bonus].c2;
    double *c3 = bonus[nlocal_bonus].c3;
    double *inertia = bonus[nlocal_bonus].inertia;
    quat[0] = buf[m++];
    quat[1] = buf[m++];
    quat[2] = buf[m++];
    quat[3] = buf[m++];
    c1[0] = buf[m++];
    c1[1] = buf[m++];
    c1[2] = buf[m++];
    c2[0] = buf[m++];
    c2[1] = buf[m++];
    c2[2] = buf[m++];
    c3[0] = buf[m++];
    c3[1] = buf[m++];
    c3[2] = buf[m++];
    inertia[0] = buf[m++];
    inertia[1] = buf[m++];
    inertia[2] = buf[m++];
    bonus[nlocal_bonus].ilocal = nlocal;
    tri[nlocal] = nlocal_bonus++;
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

void AtomVecTri::create_atom(int itype, double *coord)
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

  molecule[nlocal] = 0;
  rmass[nlocal] = 1.0;
  angmom[nlocal][0] = 0.0;
  angmom[nlocal][1] = 0.0;
  angmom[nlocal][2] = 0.0;
  tri[nlocal] = -1;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecTri::data_atom(double *coord, imageint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = ATOTAGINT(values[0]);
  molecule[nlocal] = atoi(values[1]);
  type[nlocal] = atoi(values[2]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  tri[nlocal] = atoi(values[3]);
  if (tri[nlocal] == 0) tri[nlocal] = -1;
  else if (tri[nlocal] == 1) tri[nlocal] = 0;
  else error->one(FLERR,"Invalid atom type in Atoms section of data file");

  rmass[nlocal] = atof(values[4]);
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
   unpack hybrid quantities from one tri in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecTri::data_atom_hybrid(int nlocal, char **values)
{
  molecule[nlocal] = atoi(values[0]);

  tri[nlocal] = atoi(values[1]);
  if (tri[nlocal] == 0) tri[nlocal] = -1;
  else if (tri[nlocal] == 1) tri[nlocal] = 0;
  else error->one(FLERR,"Invalid atom type in Atoms section of data file");

  rmass[nlocal] = atof(values[2]);
  if (rmass[nlocal] <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  return 3;
}

/* ----------------------------------------------------------------------
   unpack one line from Tris section of data file
------------------------------------------------------------------------- */

void AtomVecTri::data_atom_bonus(int m, char **values)
{
  if (tri[m]) error->one(FLERR,"Assigning tri parameters to non-tri atom");

  if (nlocal_bonus == nmax_bonus) grow_bonus();

  double c1[3],c2[3],c3[3];
  c1[0] = atof(values[0]);
  c1[1] = atof(values[1]);
  c1[2] = atof(values[2]);
  c2[0] = atof(values[3]);
  c2[1] = atof(values[4]);
  c2[2] = atof(values[5]);
  c3[0] = atof(values[6]);
  c3[1] = atof(values[7]);
  c3[2] = atof(values[8]);

  // check for duplicate points

  if (c1[0] == c2[0] && c1[1] == c2[1] && c1[2] == c2[2])
    error->one(FLERR,"Invalid shape in Triangles section of data file");
  if (c1[0] == c3[0] && c1[1] == c3[1] && c1[2] == c3[2])
    error->one(FLERR,"Invalid shape in Triangles section of data file");
  if (c2[0] == c3[0] && c2[1] == c3[1] && c2[2] == c3[2])
    error->one(FLERR,"Invalid shape in Triangles section of data file");

  // size = length of one edge

  double c2mc1[2],c3mc1[3];
  MathExtra::sub3(c2,c1,c2mc1);
  MathExtra::sub3(c3,c1,c3mc1);
  double size = MAX(MathExtra::len3(c2mc1),MathExtra::len3(c3mc1));

  // centroid = 1/3 of sum of vertices

  double centroid[3];
  centroid[0] = (c1[0]+c2[0]+c3[0]) / 3.0;
  centroid[1] = (c1[1]+c2[1]+c3[1]) / 3.0;
  centroid[2] = (c1[2]+c2[2]+c3[2]) / 3.0;

  double dx = centroid[0] - x[m][0];
  double dy = centroid[1] - x[m][1];
  double dz = centroid[2] - x[m][2];
  double delta = sqrt(dx*dx + dy*dy + dz*dz);

  if (delta/size > EPSILON)
    error->one(FLERR,"Inconsistent triangle in data file");

  x[m][0] = centroid[0];
  x[m][1] = centroid[1];
  x[m][2] = centroid[2];

  // reset tri mass
  // previously stored density in rmass
  // tri area = 0.5 len(U x V), where U,V are edge vectors from one vertex

  double norm[3];
  MathExtra::cross3(c2mc1,c3mc1,norm);
  double area = 0.5 * MathExtra::len3(norm);
  rmass[m] *= area;

  // inertia = inertia tensor of triangle as 6-vector in Voigt notation

  double inertia[6];
  MathExtra::inertia_triangle(c1,c2,c3,rmass[m],inertia);

  // diagonalize inertia tensor via Jacobi rotations
  // bonus[].inertia = 3 eigenvalues = principal moments of inertia
  // evectors and exzy_space = 3 evectors = principal axes of triangle

  double tensor[3][3],evectors[3][3];
  tensor[0][0] = inertia[0];
  tensor[1][1] = inertia[1];
  tensor[2][2] = inertia[2];
  tensor[1][2] = tensor[2][1] = inertia[3];
  tensor[0][2] = tensor[2][0] = inertia[4];
  tensor[0][1] = tensor[1][0] = inertia[5];

  int ierror = MathExtra::jacobi(tensor,bonus[nlocal_bonus].inertia,evectors);
  if (ierror) error->one(FLERR,"Insufficient Jacobi rotations for triangle");

  double ex_space[3],ey_space[3],ez_space[3];
  ex_space[0] = evectors[0][0];
  ex_space[1] = evectors[1][0];
  ex_space[2] = evectors[2][0];
  ey_space[0] = evectors[0][1];
  ey_space[1] = evectors[1][1];
  ey_space[2] = evectors[2][1];
  ez_space[0] = evectors[0][2];
  ez_space[1] = evectors[1][2];
  ez_space[2] = evectors[2][2];

  // enforce 3 orthogonal vectors as a right-handed coordinate system
  // flip 3rd vector if needed

  MathExtra::cross3(ex_space,ey_space,norm);
  if (MathExtra::dot3(norm,ez_space) < 0.0) MathExtra::negate3(ez_space);

  // create initial quaternion

  MathExtra::exyz_to_q(ex_space,ey_space,ez_space,bonus[nlocal_bonus].quat);

  // bonus c1,c2,c3 = displacement of c1,c2,c3 from centroid
  // in basis of principal axes

  double disp[3];
  MathExtra::sub3(c1,centroid,disp);
  MathExtra::transpose_matvec(ex_space,ey_space,ez_space,
                              disp,bonus[nlocal_bonus].c1);
  MathExtra::sub3(c2,centroid,disp);
  MathExtra::transpose_matvec(ex_space,ey_space,ez_space,
                              disp,bonus[nlocal_bonus].c2);
  MathExtra::sub3(c3,centroid,disp);
  MathExtra::transpose_matvec(ex_space,ey_space,ez_space,
                              disp,bonus[nlocal_bonus].c3);

  bonus[nlocal_bonus].ilocal = m;
  tri[m] = nlocal_bonus++;
}

/* ----------------------------------------------------------------------
   unpack one line from Velocities section of data file
------------------------------------------------------------------------- */

void AtomVecTri::data_vel(int m, char **values)
{
  v[m][0] = atof(values[0]);
  v[m][1] = atof(values[1]);
  v[m][2] = atof(values[2]);
  angmom[m][0] = atof(values[3]);
  angmom[m][1] = atof(values[4]);
  angmom[m][2] = atof(values[5]);
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one tri in Velocities section of data file
------------------------------------------------------------------------- */

int AtomVecTri::data_vel_hybrid(int m, char **values)
{
  angmom[m][0] = atof(values[0]);
  angmom[m][1] = atof(values[1]);
  angmom[m][2] = atof(values[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecTri::pack_data(double **buf)
{
  double c2mc1[2],c3mc1[3],norm[3];
  double area;

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = ubuf(molecule[i]).d;
    buf[i][2] = ubuf(type[i]).d;
    if (tri[i] < 0) buf[i][3] = ubuf(0).d;
    else buf[i][3] = ubuf(1).d;
    if (tri[i] < 0) buf[i][4] = rmass[i];
    else {
      MathExtra::sub3(bonus[tri[i]].c2,bonus[tri[i]].c1,c2mc1);
      MathExtra::sub3(bonus[tri[i]].c3,bonus[tri[i]].c1,c3mc1);
      MathExtra::cross3(c2mc1,c3mc1,norm);
      area = 0.5 * MathExtra::len3(norm);
      buf[i][4] = rmass[i]/area;
    }
    buf[i][5] = x[i][0];
    buf[i][6] = x[i][1];
    buf[i][7] = x[i][2];
    buf[i][8] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
    buf[i][9] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][10] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;
  }
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecTri::pack_data_hybrid(int i, double *buf)
{
  buf[0] = ubuf(molecule[i]).d;
  if (tri[i] < 0) buf[1] = ubuf(0).d;
  else buf[1] = ubuf(1).d;
  if (tri[i] < 0) buf[2] = rmass[i];
  else {
    double c2mc1[2],c3mc1[3],norm[3];
    MathExtra::sub3(bonus[tri[i]].c2,bonus[tri[i]].c1,c2mc1);
    MathExtra::sub3(bonus[tri[i]].c3,bonus[tri[i]].c1,c3mc1);
    MathExtra::cross3(c2mc1,c3mc1,norm);
    double area = 0.5 * MathExtra::len3(norm);
    buf[2] = rmass[i]/area;
  }
  return 3;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecTri::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT 
            " %d %d %d %-1.16e %-1.16e %-1.16e %-1.16e %d %d %d\n",
            (tagint) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,
            (int) ubuf(buf[i][2]).i,(int) ubuf(buf[i][3]).i,
            buf[i][4],buf[i][5],buf[i][6],buf[i][7],
            (int) ubuf(buf[i][8]).i,(int) ubuf(buf[i][9]).i,
            (int) ubuf(buf[i][10]).i);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecTri::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %d %d %-1.16e",(int) ubuf(buf[0]).i,(int) ubuf(buf[1]).i,buf[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   pack velocity info for data file
------------------------------------------------------------------------- */

void AtomVecTri::pack_vel(double **buf)
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

int AtomVecTri::pack_vel_hybrid(int i, double *buf)
{
  buf[0] = angmom[i][0];
  buf[1] = angmom[i][1];
  buf[2] = angmom[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
   write velocity info to data file
------------------------------------------------------------------------- */

void AtomVecTri::write_vel(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT 
            " %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e\n",
            (tagint) ubuf(buf[i][0]).i,buf[i][1],buf[i][2],buf[i][3],
            buf[i][4],buf[i][5],buf[i][6]);
}

/* ----------------------------------------------------------------------
   write hybrid velocity info to data file
------------------------------------------------------------------------- */

int AtomVecTri::write_vel_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e %-1.16e %-1.16e",buf[0],buf[1],buf[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecTri::memory_usage()
{
  bigint bytes = 0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);

  if (atom->memcheck("molecule")) bytes += memory->usage(molecule,nmax);
  if (atom->memcheck("rmass")) bytes += memory->usage(rmass,nmax);
  if (atom->memcheck("angmom")) bytes += memory->usage(angmom,nmax,3);
  if (atom->memcheck("torque")) bytes += 
                                  memory->usage(torque,nmax*comm->nthreads,3);
  if (atom->memcheck("tri")) bytes += memory->usage(tri,nmax);

  bytes += nmax_bonus*sizeof(Bonus);

  return bytes;
}
