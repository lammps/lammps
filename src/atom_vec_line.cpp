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
#include <cstring>
#include "atom_vec_line.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "force.h"
#include "fix.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define EPSILON 0.001

/* ---------------------------------------------------------------------- */

AtomVecLine::AtomVecLine(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;

  comm_x_only = comm_f_only = 0;
  size_forward = 4;
  size_reverse = 6;
  size_border = 12;
  size_velocity = 6;
  size_data_atom = 8;
  size_data_vel = 7;
  size_data_bonus = 5;
  xcol_data = 6;

  atom->line_flag = 1;
  atom->molecule_flag = atom->rmass_flag = 1;
  atom->radius_flag = atom->omega_flag = atom->torque_flag = 1;
  atom->sphere_flag = 1;

  nlocal_bonus = nghost_bonus = nmax_bonus = 0;
  bonus = NULL;
}

/* ---------------------------------------------------------------------- */

AtomVecLine::~AtomVecLine()
{
  memory->sfree(bonus);
}

/* ---------------------------------------------------------------------- */

void AtomVecLine::init()
{
  AtomVec::init();

  if (domain->dimension != 2)
    error->all(FLERR,"Atom_style line can only be used in 2d simulations");
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecLine::grow(int n)
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

  molecule = memory->grow(atom->molecule,nmax,"atom:molecule");
  rmass = memory->grow(atom->rmass,nmax,"atom:rmass");
  radius = memory->grow(atom->radius,nmax,"atom:radius");
  omega = memory->grow(atom->omega,nmax,3,"atom:omega");
  torque = memory->grow(atom->torque,nmax*comm->nthreads,3,"atom:torque");
  line = memory->grow(atom->line,nmax,"atom:line");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecLine::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;
  molecule = atom->molecule; rmass = atom->rmass;
  radius = atom->radius; omega = atom->omega; torque = atom->torque;
  line = atom->line;
}

/* ----------------------------------------------------------------------
   grow bonus data structure
------------------------------------------------------------------------- */

void AtomVecLine::grow_bonus()
{
  nmax_bonus = grow_nmax_bonus(nmax_bonus);
  if (nmax_bonus < 0)
    error->one(FLERR,"Per-processor system is too big");

  bonus = (Bonus *) memory->srealloc(bonus,nmax_bonus*sizeof(Bonus),
                                     "atom:bonus");
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecLine::copy(int i, int j, int delflag)
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
  radius[j] = radius[i];
  omega[j][0] = omega[i][0];
  omega[j][1] = omega[i][1];
  omega[j][2] = omega[i][2];

  // if deleting atom J via delflag and J has bonus data, then delete it

  if (delflag && line[j] >= 0) {
    copy_bonus(nlocal_bonus-1,line[j]);
    nlocal_bonus--;
  }

  // if atom I has bonus data, reset I's bonus.ilocal to loc J
  // do NOT do this if self-copy (I=J) since I's bonus data is already deleted

  if (line[i] >= 0 && i != j) bonus[line[i]].ilocal = j;
  line[j] = line[i];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ----------------------------------------------------------------------
   copy bonus data from I to J, effectively deleting the J entry
   also reset line that points to I to now point to J
------------------------------------------------------------------------- */

void AtomVecLine::copy_bonus(int i, int j)
{
  line[bonus[i].ilocal] = j;
  memcpy(&bonus[j],&bonus[i],sizeof(Bonus));
}

/* ----------------------------------------------------------------------
   clear ghost info in bonus data
   called before ghosts are recommunicated in comm and irregular
------------------------------------------------------------------------- */

void AtomVecLine::clear_bonus()
{
  nghost_bonus = 0;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->clear_bonus();
}

/* ----------------------------------------------------------------------
   set length value in bonus data for particle I
   oriented along x axis
   this may create or delete entry in bonus data
------------------------------------------------------------------------- */

void AtomVecLine::set_length(int i, double value)
{
  if (line[i] < 0) {
    if (value == 0.0) return;
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    bonus[nlocal_bonus].length = value;
    bonus[nlocal_bonus].theta = 0.0;
    bonus[nlocal_bonus].ilocal = i;
    line[i] = nlocal_bonus++;
  } else if (value == 0.0) {
    copy_bonus(nlocal_bonus-1,line[i]);
    nlocal_bonus--;
    line[i] = -1;
  } else bonus[line[i]].length = value;

  // also set radius = half of length
  // unless value = 0.0, then set diameter = 1.0

  radius[i] = 0.5 * value;
  if (value == 0.0) radius[i] = 0.5;
}

/* ---------------------------------------------------------------------- */

int AtomVecLine::pack_comm(int n, int *list, double *buf,
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
      if (line[j] >= 0) buf[m++] = bonus[line[j]].theta;
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
      if (line[j] >= 0) buf[m++] = bonus[line[j]].theta;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecLine::pack_comm_vel(int n, int *list, double *buf,
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
      if (line[j] >= 0) buf[m++] = bonus[line[j]].theta;
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = omega[j][0];
      buf[m++] = omega[j][1];
      buf[m++] = omega[j][2];
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
        if (line[j] >= 0) buf[m++] = bonus[line[j]].theta;
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
        buf[m++] = omega[j][0];
        buf[m++] = omega[j][1];
        buf[m++] = omega[j][2];
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
        if (line[j] >= 0) buf[m++] = bonus[line[j]].theta;
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
        buf[m++] = omega[j][0];
        buf[m++] = omega[j][1];
        buf[m++] = omega[j][2];
      }
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecLine::pack_comm_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (line[j] >= 0) buf[m++] = bonus[line[j]].theta;
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecLine::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    if (line[i] >= 0) bonus[line[i]].theta = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecLine::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    if (line[i] >= 0) bonus[line[i]].theta = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    omega[i][0] = buf[m++];
    omega[i][1] = buf[m++];
    omega[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecLine::unpack_comm_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    if (line[i] >= 0) bonus[line[i]].theta = buf[m++];
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecLine::pack_reverse(int n, int first, double *buf)
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

int AtomVecLine::pack_reverse_hybrid(int n, int first, double *buf)
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

void AtomVecLine::unpack_reverse(int n, int *list, double *buf)
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

int AtomVecLine::unpack_reverse_hybrid(int n, int *list, double *buf)
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

int AtomVecLine::pack_border(int n, int *list, double *buf,
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
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
      if (line[j] < 0) buf[m++] = ubuf(0).d;
      else {
        buf[m++] = ubuf(1).d;
        buf[m++] = bonus[line[j]].length;
        buf[m++] = bonus[line[j]].theta;
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
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
      if (line[j] < 0) buf[m++] = ubuf(0).d;
      else {
        buf[m++] = ubuf(1).d;
        buf[m++] = bonus[line[j]].length;
        buf[m++] = bonus[line[j]].theta;
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecLine::pack_border_vel(int n, int *list, double *buf,
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
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
      if (line[j] < 0) buf[m++] = ubuf(0).d;
      else {
        buf[m++] = ubuf(1).d;
        buf[m++] = bonus[line[j]].length;
        buf[m++] = bonus[line[j]].theta;
      }
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = omega[j][0];
      buf[m++] = omega[j][1];
      buf[m++] = omega[j][2];
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
        buf[m++] = radius[j];
        buf[m++] = rmass[j];
        if (line[j] < 0) buf[m++] = ubuf(0).d;
        else {
          buf[m++] = ubuf(1).d;
          buf[m++] = bonus[line[j]].length;
          buf[m++] = bonus[line[j]].theta;
        }
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
        buf[m++] = omega[j][0];
        buf[m++] = omega[j][1];
        buf[m++] = omega[j][2];
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
        buf[m++] = radius[j];
        buf[m++] = rmass[j];
        if (line[j] < 0) buf[m++] = ubuf(0).d;
        else {
          buf[m++] = ubuf(1).d;
          buf[m++] = bonus[line[j]].length;
          buf[m++] = bonus[line[j]].theta;
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
        buf[m++] = omega[j][0];
        buf[m++] = omega[j][1];
        buf[m++] = omega[j][2];
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecLine::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = ubuf(molecule[j]).d;
    buf[m++] = radius[j];
    buf[m++] = rmass[j];
    if (line[j] < 0) buf[m++] = ubuf(0).d;
    else {
      buf[m++] = ubuf(1).d;
      buf[m++] = bonus[line[j]].length;
      buf[m++] = bonus[line[j]].theta;
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecLine::unpack_border(int n, int first, double *buf)
{
  int i,j,m,last;

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
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
    line[i] = (int) ubuf(buf[m++]).i;
    if (line[i] == 0) line[i] = -1;
    else {
      j = nlocal_bonus + nghost_bonus;
      if (j == nmax_bonus) grow_bonus();
      bonus[j].length = buf[m++];
      bonus[j].theta = buf[m++];
      bonus[j].ilocal = i;
      line[i] = j;
      nghost_bonus++;
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecLine::unpack_border_vel(int n, int first, double *buf)
{
  int i,j,m,last;

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
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
    line[i] = (int) ubuf(buf[m++]).i;
    if (line[i] == 0) line[i] = -1;
    else {
      j = nlocal_bonus + nghost_bonus;
      if (j == nmax_bonus) grow_bonus();
      bonus[j].length = buf[m++];
      bonus[j].theta = buf[m++];
      bonus[j].ilocal = i;
      line[i] = j;
      nghost_bonus++;
    }
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    omega[i][0] = buf[m++];
    omega[i][1] = buf[m++];
    omega[i][2] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

int AtomVecLine::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,j,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    molecule[i] = (tagint) ubuf(buf[m++]).i;
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
    line[i] = (int) ubuf(buf[m++]).i;
    if (line[i] == 0) line[i] = -1;
    else {
      j = nlocal_bonus + nghost_bonus;
      if (j == nmax_bonus) grow_bonus();
      bonus[j].length = buf[m++];
      bonus[j].theta = buf[m++];
      bonus[j].ilocal = i;
      line[i] = j;
      nghost_bonus++;
    }
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecLine::pack_exchange(int i, double *buf)
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
  buf[m++] = radius[i];
  buf[m++] = omega[i][0];
  buf[m++] = omega[i][1];
  buf[m++] = omega[i][2];

  if (line[i] < 0) buf[m++] = ubuf(0).d;
  else {
    buf[m++] = ubuf(1).d;
    int j = line[i];
    buf[m++] = bonus[j].length;
    buf[m++] = bonus[j].theta;
  }

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecLine::unpack_exchange(double *buf)
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

  molecule[nlocal] = (tagint) ubuf(buf[m++]).i;
  rmass[nlocal] = buf[m++];
  radius[nlocal] = buf[m++];
  omega[nlocal][0] = buf[m++];
  omega[nlocal][1] = buf[m++];
  omega[nlocal][2] = buf[m++];

  line[nlocal] = (int) ubuf(buf[m++]).i;
  if (line[nlocal] == 0) line[nlocal] = -1;
  else {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    bonus[nlocal_bonus].length = buf[m++];
    bonus[nlocal_bonus].theta = buf[m++];
    bonus[nlocal_bonus].ilocal = nlocal;
    line[nlocal] = nlocal_bonus++;
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

int AtomVecLine::size_restart()
{
  int i;

  int n = 0;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++)
    if (line[i] >= 0) n += 20;
    else n += 18;

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

int AtomVecLine::pack_restart(int i, double *buf)
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
  buf[m++] = radius[i];
  buf[m++] = omega[i][0];
  buf[m++] = omega[i][1];
  buf[m++] = omega[i][2];

  if (line[i] < 0) buf[m++] = ubuf(0).d;
  else {
    buf[m++] = ubuf(1).d;
    int j = line[i];
    buf[m++] = bonus[j].length;
    buf[m++] = bonus[j].theta;
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

int AtomVecLine::unpack_restart(double *buf)
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

  molecule[nlocal] = (tagint) ubuf(buf[m++]).i;
  rmass[nlocal] = buf[m++];
  radius[nlocal] = buf[m++];
  omega[nlocal][0] = buf[m++];
  omega[nlocal][1] = buf[m++];
  omega[nlocal][2] = buf[m++];

  line[nlocal] = (int) ubuf(buf[m++]).i;
  if (line[nlocal] == 0) line[nlocal] = -1;
  else {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    bonus[nlocal_bonus].length = buf[m++];
    bonus[nlocal_bonus].theta = buf[m++];
    bonus[nlocal_bonus].ilocal = nlocal;
    line[nlocal] = nlocal_bonus++;
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

void AtomVecLine::create_atom(int itype, double *coord)
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
  radius[nlocal] = 0.5;
  rmass[nlocal] = 4.0*MY_PI/3.0 * radius[nlocal]*radius[nlocal]*radius[nlocal];
  omega[nlocal][0] = 0.0;
  omega[nlocal][1] = 0.0;
  omega[nlocal][2] = 0.0;
  line[nlocal] = -1;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecLine::data_atom(double *coord, imageint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = ATOTAGINT(values[0]);
  molecule[nlocal] = ATOTAGINT(values[1]);
  type[nlocal] = atoi(values[2]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  line[nlocal] = atoi(values[3]);
  if (line[nlocal] == 0) line[nlocal] = -1;
  else if (line[nlocal] == 1) line[nlocal] = 0;
  else error->one(FLERR,"Invalid atom type in Atoms section of data file");

  rmass[nlocal] = atof(values[4]);
  if (rmass[nlocal] <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  if (line[nlocal] < 0) {
    radius[nlocal] = 0.5;
    rmass[nlocal] *= 4.0*MY_PI/3.0 *
      radius[nlocal]*radius[nlocal]*radius[nlocal];
  } else radius[nlocal] = 0.0;

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  omega[nlocal][0] = 0.0;
  omega[nlocal][1] = 0.0;
  omega[nlocal][2] = 0.0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecLine::data_atom_hybrid(int nlocal, char **values)
{
  molecule[nlocal] = ATOTAGINT(values[0]);

  line[nlocal] = atoi(values[1]);
  if (line[nlocal] == 0) line[nlocal] = -1;
  else if (line[nlocal] == 1) line[nlocal] = 0;
  else error->one(FLERR,"Invalid atom type in Atoms section of data file");

  rmass[nlocal] = atof(values[2]);
  if (rmass[nlocal] <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  if (line[nlocal] < 0) {
    radius[nlocal] = 0.5;
    rmass[nlocal] *= 4.0*MY_PI/3.0 *
      radius[nlocal]*radius[nlocal]*radius[nlocal];
  } else radius[nlocal] = 0.0;

  return 3;
}

/* ----------------------------------------------------------------------
   unpack one line from Lines section of data file
------------------------------------------------------------------------- */

void AtomVecLine::data_atom_bonus(int m, char **values)
{
  if (line[m]) error->one(FLERR,"Assigning line parameters to non-line atom");

  if (nlocal_bonus == nmax_bonus) grow_bonus();

  double x1 = atof(values[0]);
  double y1 = atof(values[1]);
  double x2 = atof(values[2]);
  double y2 = atof(values[3]);
  double dx = x2 - x1;
  double dy = y2 - y1;
  double length = sqrt(dx*dx + dy*dy);

  bonus[nlocal_bonus].length = length;
  if (dy >= 0.0) bonus[nlocal_bonus].theta = acos(dx/length);
  else bonus[nlocal_bonus].theta = -acos(dx/length);

  double xc = 0.5*(x1+x2);
  double yc = 0.5*(y1+y2);
  dx = xc - x[m][0];
  dy = yc - x[m][1];
  double delta = sqrt(dx*dx + dy*dy);

  if (delta/length > EPSILON)
    error->one(FLERR,"Inconsistent line segment in data file");

  x[m][0] = xc;
  x[m][1] = yc;

  // reset line radius and mass
  // rmass currently holds density

  radius[m] = 0.5 * length;
  rmass[m] *= length;

  bonus[nlocal_bonus].ilocal = m;
  line[m] = nlocal_bonus++;
}

/* ----------------------------------------------------------------------
   unpack one line from Velocities section of data file
------------------------------------------------------------------------- */

void AtomVecLine::data_vel(int m, char **values)
{
  v[m][0] = atof(values[0]);
  v[m][1] = atof(values[1]);
  v[m][2] = atof(values[2]);
  omega[m][0] = atof(values[3]);
  omega[m][1] = atof(values[4]);
  omega[m][2] = atof(values[5]);
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Velocities section of data file
------------------------------------------------------------------------- */

int AtomVecLine::data_vel_hybrid(int m, char **values)
{
  omega[m][0] = atof(values[0]);
  omega[m][1] = atof(values[1]);
  omega[m][2] = atof(values[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecLine::pack_data(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = ubuf(molecule[i]).d;
    buf[i][2] = ubuf(type[i]).d;
    if (line[i] < 0) buf[i][3] = ubuf(0).d;
    else buf[i][3] = ubuf(1).d;
    if (line[i] < 0)
      buf[i][4] = rmass[i] / (4.0*MY_PI/3.0 * radius[i]*radius[i]*radius[i]);
    else buf[i][4] = rmass[i]/bonus[line[i]].length;
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

int AtomVecLine::pack_data_hybrid(int i, double *buf)
{
  buf[0] = ubuf(molecule[i]).d;
  if (line[i] < 0) buf[1] = ubuf(0).d;
  else buf[1] = ubuf(1).d;
  if (line[i] < 0)
    buf[2] = rmass[i] / (4.0*MY_PI/3.0 * radius[i]*radius[i]*radius[i]);
  else buf[2] = rmass[i]/bonus[line[i]].length;
  return 3;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecLine::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT " " TAGINT_FORMAT
            " %d %d %-1.16e %-1.16e %-1.16e %-1.16e %d %d %d\n",
            (tagint) ubuf(buf[i][0]).i,(tagint) ubuf(buf[i][1]).i,
            (int) ubuf(buf[i][2]).i,(int) ubuf(buf[i][3]).i,
            buf[i][4],buf[i][5],buf[i][6],buf[i][7],
            (int) ubuf(buf[i][8]).i,(int) ubuf(buf[i][9]).i,
            (int) ubuf(buf[i][10]).i);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecLine::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," " TAGINT_FORMAT " %d %-1.16e",
          (tagint) ubuf(buf[0]).i,(int) ubuf(buf[1]).i,buf[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   pack velocity info for data file
------------------------------------------------------------------------- */

void AtomVecLine::pack_vel(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = v[i][0];
    buf[i][2] = v[i][1];
    buf[i][3] = v[i][2];
    buf[i][4] = omega[i][0];
    buf[i][5] = omega[i][1];
    buf[i][6] = omega[i][2];
  }
}

/* ----------------------------------------------------------------------
   pack hybrid velocity info for data file
------------------------------------------------------------------------- */

int AtomVecLine::pack_vel_hybrid(int i, double *buf)
{
  buf[0] = omega[i][0];
  buf[1] = omega[i][1];
  buf[2] = omega[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
   write velocity info to data file
------------------------------------------------------------------------- */

void AtomVecLine::write_vel(FILE *fp, int n, double **buf)
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

int AtomVecLine::write_vel_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e %-1.16e %-1.16e",buf[0],buf[1],buf[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecLine::memory_usage()
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
  if (atom->memcheck("radius")) bytes += memory->usage(radius,nmax);
  if (atom->memcheck("omega")) bytes += memory->usage(omega,nmax,3);
  if (atom->memcheck("torque"))
    bytes += memory->usage(torque,nmax*comm->nthreads,3);
  if (atom->memcheck("line")) bytes += memory->usage(line,nmax);

  bytes += nmax_bonus*sizeof(Bonus);

  return bytes;
}

/* ----------------------------------------------------------------------
   check consistency of internal Bonus data structure
   n = # of atoms in regular structure to check against
------------------------------------------------------------------------- */

/*
void AtomVecLine::consistency_check(int n, char *str)
{
  int iflag = 0;
  int count = 0;
  for (int i = 0; i < n; i++) {

    if (line[i] >= 0) {
      count++;
      if (line[i] >= nlocal_bonus) iflag++;
      if (bonus[line[i]].ilocal != i) iflag++;
      //if (comm->me == 1 && update->ntimestep == 873)
      //        printf("CCHK %s: %d %d: %d %d: %d %d\n",
      //       str,i,n,line[i],nlocal_bonus,bonus[line[i]].ilocal,iflag);
    }
  }

  if (iflag) {
    printf("BAD vecline ptrs: %s: %d %d: %d\n",str,comm->me,
           update->ntimestep,iflag);
    MPI_Abort(world,1);
  }

  if (count != nlocal_bonus) {
    char msg[128];
    printf("BAD vecline count: %s: %d %d: %d %d\n",
           str,comm->me,update->ntimestep,count,nlocal_bonus);
    MPI_Abort(world,1);
  }
}
*/
