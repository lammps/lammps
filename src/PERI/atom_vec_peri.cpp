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
   Contributing author: Mike Parks (SNL)
------------------------------------------------------------------------- */

#include "float.h"
#include "stdlib.h"
#include "atom_vec_peri.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 10000

/* ---------------------------------------------------------------------- */

AtomVecPeri::AtomVecPeri(LAMMPS *lmp, int narg, char **arg) : 
  AtomVec(lmp, narg, arg) 
{
  molecular = 0;

  comm_x_only = 0;
  comm_f_only = 1;
  size_forward = 4;
  size_reverse = 3;
  size_border = 11;
  size_velocity = 3;
  size_data_atom = 7;
  size_data_vel = 4;
  xcol_data = 5;

  atom->peri_flag = 1;
  atom->vfrac_flag = atom->rmass_flag = 1;
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */
 
void AtomVecPeri::grow(int n)
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

  vfrac = memory->grow(atom->vfrac,nmax,"atom:vfrac");
  rmass = memory->grow(atom->rmass,nmax,"atom:rmass");
  s0 = memory->grow(atom->s0,nmax,"atom:s0");
  x0 = memory->grow(atom->x0,nmax,3,"atom:x0");
 
  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecPeri::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;
  vfrac = atom->vfrac; rmass = atom->rmass;
  s0 = atom->s0; x0 = atom->x0;
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecPeri::copy(int i, int j, int delflag)
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

  vfrac[j] = vfrac[i];
  rmass[j] = rmass[i];
  s0[j] = s0[i];
  x0[j][0] = x0[i][0];
  x0[j][1] = x0[i][1];
  x0[j][2] = x0[i][2];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++) 
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j);
}

/* ---------------------------------------------------------------------- */

int AtomVecPeri::pack_comm(int n, int *list, double *buf,
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
      buf[m++] = s0[j];
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
      buf[m++] = s0[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecPeri::pack_comm_vel(int n, int *list, double *buf,
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
      buf[m++] = s0[j];
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
	buf[m++] = s0[j];
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
	buf[m++] = s0[j];
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

int AtomVecPeri::pack_comm_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = s0[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecPeri::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    s0[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecPeri::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    s0[i] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecPeri::unpack_comm_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    s0[i] = buf[m++];
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecPeri::pack_reverse(int n, int first, double *buf)
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

void AtomVecPeri::unpack_reverse(int n, int *list, double *buf)
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

int AtomVecPeri::pack_border(int n, int *list, double *buf,
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
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = vfrac[j];
      buf[m++] = s0[j];
      buf[m++] = x0[j][0];
      buf[m++] = x0[j][1];
      buf[m++] = x0[j][2];
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
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = vfrac[j];
      buf[m++] = s0[j];
      buf[m++] = x0[j][0];
      buf[m++] = x0[j][1];
      buf[m++] = x0[j][2];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecPeri::pack_border_vel(int n, int *list, double *buf,
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
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = vfrac[j];
      buf[m++] = s0[j];
      buf[m++] = x0[j][0];
      buf[m++] = x0[j][1];
      buf[m++] = x0[j][2];
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
	buf[m++] = tag[j];
	buf[m++] = type[j];
	buf[m++] = mask[j];
	buf[m++] = vfrac[j];
	buf[m++] = s0[j];
	buf[m++] = x0[j][0];
	buf[m++] = x0[j][1];
	buf[m++] = x0[j][2];
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
	buf[m++] = tag[j];
	buf[m++] = type[j];
	buf[m++] = mask[j];
	buf[m++] = vfrac[j];
	buf[m++] = s0[j];
	buf[m++] = x0[j][0];
	buf[m++] = x0[j][1];
	buf[m++] = x0[j][2];
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

int AtomVecPeri::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = vfrac[j];
    buf[m++] = s0[j];
    buf[m++] = x0[j][0];
    buf[m++] = x0[j][1];
    buf[m++] = x0[j][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecPeri::unpack_border(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = static_cast<int> (buf[m++]);
    type[i] = static_cast<int> (buf[m++]);
    mask[i] = static_cast<int> (buf[m++]);
    vfrac[i] = buf[m++];
    s0[i] = buf[m++];
    x0[i][0] = buf[m++];
    x0[i][1] = buf[m++];
    x0[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecPeri::unpack_border_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = static_cast<int> (buf[m++]);
    type[i] = static_cast<int> (buf[m++]);
    mask[i] = static_cast<int> (buf[m++]);
    vfrac[i] = buf[m++];
    s0[i] = buf[m++];
    x0[i][0] = buf[m++];
    x0[i][1] = buf[m++];
    x0[i][2] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */
 
int AtomVecPeri::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    vfrac[i] = buf[m++];
    s0[i] = buf[m++];
    x0[i][0] = buf[m++];
    x0[i][1] = buf[m++];
    x0[i][2] = buf[m++];
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecPeri::pack_exchange(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = tag[i];
  buf[m++] = type[i];
  buf[m++] = mask[i];
  buf[m++] = image[i];

  buf[m++] = vfrac[i];
  buf[m++] = rmass[i];
  buf[m++] = s0[i];
  buf[m++] = x0[i][0];
  buf[m++] = x0[i][1];
  buf[m++] = x0[i][2];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++) 
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecPeri::unpack_exchange(double *buf)
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
  tag[nlocal] = static_cast<int> (buf[m++]);
  type[nlocal] = static_cast<int> (buf[m++]);
  mask[nlocal] = static_cast<int> (buf[m++]);
  image[nlocal] = static_cast<int> (buf[m++]);

  vfrac[nlocal] = buf[m++];
  rmass[nlocal] = buf[m++];
  s0[nlocal] = buf[m++];
  x0[nlocal][0] = buf[m++];
  x0[nlocal][1] = buf[m++];
  x0[nlocal][2] = buf[m++];

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
 
int AtomVecPeri::size_restart()
{
  int i;
 
  int nlocal = atom->nlocal;
  int n = 17 * nlocal;
 
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
 
int AtomVecPeri::pack_restart(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = tag[i];
  buf[m++] = type[i];
  buf[m++] = mask[i];
  buf[m++] = image[i];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
 
  buf[m++] = vfrac[i];
  buf[m++] = rmass[i];
  buf[m++] = s0[i];
  buf[m++] = x0[i][0];
  buf[m++] = x0[i][1];
  buf[m++] = x0[i][2];
 
  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);
 
  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */
 
int AtomVecPeri::unpack_restart(double *buf)
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
  tag[nlocal] = static_cast<int> (buf[m++]);
  type[nlocal] = static_cast<int> (buf[m++]);
  mask[nlocal] = static_cast<int> (buf[m++]);
  image[nlocal] = static_cast<int> (buf[m++]);
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
 
  vfrac[nlocal] = buf[m++];
  rmass[nlocal] = buf[m++];
  s0[nlocal] = buf[m++];
  x0[nlocal][0] = buf[m++];
  x0[nlocal][1] = buf[m++];
  x0[nlocal][2] = buf[m++];
 
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
 
void AtomVecPeri::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);
 
  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  mask[nlocal] = 1;
  image[nlocal] = (512 << 20) | (512 << 10) | 512;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
 
  vfrac[nlocal] = 1.0;
  rmass[nlocal] = 1.0;
  s0[nlocal] = DBL_MAX;
  x0[nlocal][0] = coord[0];
  x0[nlocal][1] = coord[1];
  x0[nlocal][2] = coord[2];
 
  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */
 
void AtomVecPeri::data_atom(double *coord, int imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = atoi(values[0]);
  if (tag[nlocal] <= 0)
    error->one(FLERR,"Invalid atom ID in Atoms section of data file");
 
  type[nlocal] = atoi(values[1]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  vfrac[nlocal] = atof(values[2]);
  rmass[nlocal] = atof(values[3]);
  if (rmass[nlocal] <= 0.0) error->one(FLERR,"Invalid mass value");

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;
 
  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  s0[nlocal] = DBL_MAX;
  x0[nlocal][0] = coord[0];
  x0[nlocal][1] = coord[1];
  x0[nlocal][2] = coord[2];
 
  atom->nlocal++;
}
 

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */
 
int AtomVecPeri::data_atom_hybrid(int nlocal, char **values)
{
  vfrac[nlocal] = atof(values[0]);
  rmass[nlocal] = atof(values[1]);
  if (rmass[nlocal] <= 0.0) error->one(FLERR,"Invalid mass value");

  s0[nlocal] = DBL_MAX;
  x0[nlocal][0] = x[nlocal][0];
  x0[nlocal][1] = x[nlocal][1];
  x0[nlocal][2] = x[nlocal][2];
 
  return 2;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */
 
bigint AtomVecPeri::memory_usage()
{
  bigint bytes = 0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);

  if (atom->memcheck("vfrac")) bytes += memory->usage(vfrac,nmax);
  if (atom->memcheck("rmass")) bytes += memory->usage(rmass,nmax);
  if (atom->memcheck("s0")) bytes += memory->usage(s0,nmax);
  if (atom->memcheck("x0")) bytes += memory->usage(x0,nmax,3);
 
  return bytes;
}
