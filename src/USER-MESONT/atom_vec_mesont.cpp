/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing author: Maxim Shugaev (UVA), mvs9t@virginia.edu
------------------------------------------------------------------------- */

#include "atom_vec_mesont.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecMesoNT::AtomVecMesoNT(LAMMPS *lmp) : AtomVec(lmp)
{
  mass_type = 1;
  atom->mesont_flag = 1;
  molecular = 0;

  comm_x_only = comm_f_only = 1;
  size_forward = 3;
  size_reverse = 3;
  size_border = 13;
  size_velocity = 3;
  size_data_atom = 12;
  size_data_vel = 4;
  xcol_data = 10;
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecMesoNT::grow(int n)
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
  x = memory->grow(atom->x,nmax,3,"atom:x");
  v = memory->grow(atom->v,nmax,3,"atom:v");
  f = memory->grow(atom->f,nmax*comm->nthreads,3,"atom:f");

  rmass = memory->grow(atom->rmass,nmax,"atom:rmass");
  radius = memory->grow(atom->radius,nmax,"atom:radius");
  length = memory->grow(atom->length,nmax,"atom:length");
  buckling = memory->grow(atom->buckling,nmax,"atom:buckling");
  molecule = memory->grow(atom->molecule,nmax,"atom:molecule");
  bond_nt = memory->grow(atom->bond_nt,nmax,2,"atom:bond_nt");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecMesoNT::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;

  rmass = atom->rmass;
  radius = atom->radius;
  length = atom->length;
  buckling = atom->buckling;
  molecule = atom->molecule;
  bond_nt = atom->bond_nt;
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecMesoNT::copy(int i, int j, int delflag)
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
  radius[j] = radius[i];
  length[j] = length[i];
  buckling[j] = buckling[i];
  molecule[j] = molecule[i];
  bond_nt[j][0] = bond_nt[i][0];
  bond_nt[j][1] = bond_nt[i][1];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ---------------------------------------------------------------------- */

int AtomVecMesoNT::pack_comm(int n, int *list, double *buf,
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

int AtomVecMesoNT::pack_comm_vel(int n, int *list, double *buf,
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

void AtomVecMesoNT::unpack_comm(int n, int first, double *buf)
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

void AtomVecMesoNT::unpack_comm_vel(int n, int first, double *buf)
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

int AtomVecMesoNT::pack_reverse(int n, int first, double *buf)
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

void AtomVecMesoNT::unpack_reverse(int n, int *list, double *buf)
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

int AtomVecMesoNT::pack_border(int n, int *list, double *buf,
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

      buf[m++] = rmass[j];
      buf[m++] = radius[j];
      buf[m++] = length[j];
      buf[m++] = ubuf(buckling[j]).d;
      buf[m++] = ubuf(molecule[j]).d;
      buf[m++] = ubuf(bond_nt[j][0]).d;
      buf[m++] = ubuf(bond_nt[j][1]).d;
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

      buf[m++] = rmass[j];
      buf[m++] = radius[j];
      buf[m++] = length[j];
      buf[m++] = ubuf(buckling[j]).d;
      buf[m++] = ubuf(molecule[j]).d;
      buf[m++] = ubuf(bond_nt[j][0]).d;
      buf[m++] = ubuf(bond_nt[j][1]).d;
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecMesoNT::pack_border_vel(int n, int *list, double *buf,
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
      buf[m++] = rmass[j];
      buf[m++] = radius[j];
      buf[m++] = length[j];
      buf[m++] = ubuf(buckling[j]).d;
      buf[m++] = ubuf(molecule[j]).d;
      buf[m++] = ubuf(bond_nt[j][0]).d;
      buf[m++] = ubuf(bond_nt[j][1]).d;
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
        buf[m++] = rmass[j];
        buf[m++] = radius[j];
        buf[m++] = length[j];
        buf[m++] = ubuf(buckling[j]).d;
        buf[m++] = ubuf(molecule[j]).d;
        buf[m++] = ubuf(bond_nt[j][0]).d;
        buf[m++] = ubuf(bond_nt[j][1]).d;
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
        buf[m++] = rmass[j];
        buf[m++] = radius[j];
        buf[m++] = length[j];
        buf[m++] = ubuf(buckling[j]).d;
        buf[m++] = ubuf(molecule[j]).d;
        buf[m++] = ubuf(bond_nt[j][0]).d;
        buf[m++] = ubuf(bond_nt[j][1]).d;
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
      m += modify->fix[atom->extra_border[iextra]]->
        pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecMesoNT::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = rmass[j];
    buf[m++] = radius[j];
    buf[m++] = length[j];
    buf[m++] = ubuf(buckling[j]).d;
    buf[m++] = ubuf(molecule[j]).d;
    buf[m++] = ubuf(bond_nt[j][0]).d;
    buf[m++] = ubuf(bond_nt[j][1]).d;
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecMesoNT::unpack_border(int n, int first, double *buf)
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

    rmass[i] = buf[m++];
    radius[i] = buf[m++];
    length[i] = buf[m++];
    buckling[i] = (int) ubuf(buf[m++]).i;
    molecule[i] = (tagint) ubuf(buf[m++]).i;
    bond_nt[i][0] = (tagint) ubuf(buf[m++]).i;
    bond_nt[i][1] = (tagint) ubuf(buf[m++]).i;
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecMesoNT::unpack_border_vel(int n, int first, double *buf)
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
    rmass[i] = buf[m++];
    radius[i] = buf[m++];
    length[i] = buf[m++];
    buckling[i] = (int) ubuf(buf[m++]).i;
    molecule[i] = (tagint) ubuf(buf[m++]).i;
    bond_nt[i][0] = (tagint) ubuf(buf[m++]).i;
    bond_nt[i][1] = (tagint) ubuf(buf[m++]).i;
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

int AtomVecMesoNT::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    rmass[i] = buf[m++];
    radius[i] = buf[m++];
    length[i] = buf[m++];
    buckling[i] = (int) ubuf(buf[m++]).i;
    molecule[i] = (tagint) ubuf(buf[m++]).i;
    bond_nt[i][0] = (tagint) ubuf(buf[m++]).i;
    bond_nt[i][1] = (tagint) ubuf(buf[m++]).i;
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecMesoNT::pack_exchange(int i, double *buf)
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
  buf[m++] = radius[i];
  buf[m++] = length[i];
  buf[m++] = ubuf(buckling[i]).d;
  buf[m++] = ubuf(molecule[i]).d;
  buf[m++] = ubuf(bond_nt[i][0]).d;
  buf[m++] = ubuf(bond_nt[i][1]).d;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecMesoNT::unpack_exchange(double *buf)
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

  rmass[nlocal] = buf[m++];
  radius[nlocal] = buf[m++];
  length[nlocal] = buf[m++];
  buckling[nlocal] = (int) ubuf(buf[m++]).i;
  molecule[nlocal] = (tagint) ubuf(buf[m++]).i;
  bond_nt[nlocal][0] = (tagint) ubuf(buf[m++]).i;
  bond_nt[nlocal][1] = (tagint) ubuf(buf[m++]).i;

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

int AtomVecMesoNT::size_restart()
{
  int i;

  int nlocal = atom->nlocal;
  int n = 18 * nlocal;

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

int AtomVecMesoNT::pack_restart(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  buf[m++] = rmass[i];
  buf[m++] = radius[i];
  buf[m++] = length[i];
  buf[m++] = ubuf(buckling[i]).d;
  buf[m++] = ubuf(molecule[i]).d;
  buf[m++] = ubuf(bond_nt[i][0]).d;
  buf[m++] = ubuf(bond_nt[i][1]).d;
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecMesoNT::unpack_restart(double *buf)
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
  rmass[nlocal] = buf[m++];
  radius[nlocal] = buf[m++];
  length[nlocal] = buf[m++];
  buckling[nlocal] = (int) ubuf(buf[m++]).i;
  molecule[nlocal] = (tagint) ubuf(buf[m++]).i;
  bond_nt[nlocal][0] = (tagint) ubuf(buf[m++]).i;
  bond_nt[nlocal][1] = (tagint) ubuf(buf[m++]).i;
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

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

void AtomVecMesoNT::create_atom(int itype, double *coord)
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
  rmass[nlocal] = 1.0;
  radius[nlocal] = 1.0;
  length[nlocal] = 1.0;
  buckling[nlocal] = 0;
  molecule[nlocal] = 0;
  bond_nt[nlocal][0] = -1;
  bond_nt[nlocal][1] = -1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecMesoNT::data_atom(double *coord, imageint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = utils::tnumeric(FLERR,values[0],true,lmp);
  molecule[nlocal] = utils::tnumeric(FLERR,values[1],true,lmp);
  type[nlocal] = utils::inumeric(FLERR,values[2],true,lmp);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");
  bond_nt[nlocal][0] = utils::inumeric(FLERR,values[3],true,lmp);
  bond_nt[nlocal][1] = utils::inumeric(FLERR,values[4],true,lmp);
  rmass[nlocal] = utils::numeric(FLERR,values[5],true,lmp);
  radius[nlocal] = utils::numeric(FLERR,values[6],true,lmp);
  length[nlocal] = utils::numeric(FLERR,values[7],true,lmp);
  buckling[nlocal] = utils::numeric(FLERR,values[8],true,lmp);

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

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

int AtomVecMesoNT::data_atom_hybrid(int nlocal, char **values)
{
  molecule[nlocal] = utils::tnumeric(FLERR,values[0],true,lmp);
  bond_nt[nlocal][0] = utils::inumeric(FLERR,values[1],true,lmp);
  bond_nt[nlocal][1] = utils::inumeric(FLERR,values[2],true,lmp);
  rmass[nlocal] = utils::numeric(FLERR,values[3],true,lmp);
  radius[nlocal] = utils::numeric(FLERR,values[4],true,lmp);
  length[nlocal] = utils::numeric(FLERR,values[5],true,lmp);
  buckling[nlocal] = utils::numeric(FLERR,values[6],true,lmp);

  return 7;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecMesoNT::pack_data(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    int m = 0;
    buf[i][m++] = ubuf(tag[i]).d;
    buf[i][m++] = ubuf(molecule[i]).d;
    buf[i][m++] = ubuf(type[i]).d;
    buf[i][m++] = ubuf(bond_nt[i][0]).d;
    buf[i][m++] = ubuf(bond_nt[i][1]).d;
    buf[i][m++] = rmass[i];
    buf[i][m++] = radius[i];
    buf[i][m++] = length[i];
    buf[i][m++] = ubuf(buckling[i]).d;
    buf[i][m++] = x[i][0];
    buf[i][m++] = x[i][1];
    buf[i][m++] = x[i][2];
    buf[i][m++] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
    buf[i][m++] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][m++] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;
  }
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecMesoNT::pack_data_hybrid(int i, double *buf)
{
  int m = 0;
  buf[m++] = ubuf(molecule[i]).d;
  buf[m++] = ubuf(bond_nt[i][0]).d;
  buf[m++] = ubuf(bond_nt[i][1]).d;
  buf[m++] = rmass[i];
  buf[m++] = radius[i];
  buf[m++] = length[i];
  buf[m++] = ubuf(buckling[i]).d;
  return m;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecMesoNT::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp, TAGINT_FORMAT " " TAGINT_FORMAT " %d "
      TAGINT_FORMAT " " TAGINT_FORMAT
      " %-1.16e %-1.16e %-1.16e %d %-1.16e %-1.16e %-1.16e %d %d %d\n",
      (tagint) ubuf(buf[i][0]).i, (tagint) ubuf(buf[i][1]).i,
      (int) ubuf(buf[i][2]).i, (tagint) ubuf(buf[i][3]).i,
      (tagint) ubuf(buf[i][4]).i,
      buf[i][5], buf[i][6], buf[i][7], (int) ubuf(buf[i][8]).i,
      buf[i][9], buf[i][10], buf[i][11],
      (int) ubuf(buf[i][12]).i, (int) ubuf(buf[i][13]).i,
      (int) ubuf(buf[i][14]).i);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecMesoNT::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT
    " %-1.16e %-1.16e %-1.16e %d",
    (tagint) ubuf(buf[0]).i, (tagint) ubuf(buf[1]).i, (tagint) ubuf(buf[2]).i,
    buf[3], buf[4], buf[5], (int)ubuf(buf[6]).i);
  return 7;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecMesoNT::memory_usage()
{
  bigint bytes = 0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);

  if (atom->memcheck("radius")) bytes += memory->usage(radius,nmax);
  if (atom->memcheck("rmass")) bytes += memory->usage(rmass,nmax);
  if (atom->memcheck("length")) bytes += memory->usage(length,nmax);
  if (atom->memcheck("buckling")) bytes += memory->usage(buckling,nmax);
  if (atom->memcheck("molecule")) bytes += memory->usage(molecule,nmax);
  if (atom->memcheck("bond_nt")) bytes += memory->usage(bond_nt,nmax,2);

  return bytes;
}
