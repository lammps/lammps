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
#include "atom_vec_body.h"
#include "style_body.h"
#include "body.h"
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

/* ---------------------------------------------------------------------- */

AtomVecBody::AtomVecBody(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;

  // size_forward and size_border set in settings(), via Body class

  comm_x_only = comm_f_only = 0;
  size_forward = 0;
  size_reverse = 6;
  size_border = 0;
  size_velocity = 6;
  size_data_atom = 7;
  size_data_vel = 7;
  xcol_data = 5;

  atom->body_flag = 1;
  atom->rmass_flag = 1;
  atom->angmom_flag = atom->torque_flag = 1;

  nlocal_bonus = nghost_bonus = nmax_bonus = 0;
  bonus = NULL;

  bptr = NULL;

  nargcopy = 0;
  argcopy = NULL;
  copyflag = 1;

  if (sizeof(double) == sizeof(int)) intdoubleratio = 1;
  else if (sizeof(double) == 2*sizeof(int)) intdoubleratio = 2;
  else error->all(FLERR,"Internal error in atom_style body");
}

/* ---------------------------------------------------------------------- */

AtomVecBody::~AtomVecBody()
{
  int nall = nlocal_bonus + nghost_bonus;
  for (int i = 0; i < nall; i++) {
    icp->put(bonus[i].iindex);
    dcp->put(bonus[i].dindex);
  }
  memory->sfree(bonus);

  delete bptr;

  for (int i = 0; i < nargcopy; i++) delete [] argcopy[i];
  delete [] argcopy;
}

/* ----------------------------------------------------------------------
   process additional args
   instantiate Body class
   set size_forward and size_border to max sizes
------------------------------------------------------------------------- */

void AtomVecBody::settings(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Invalid atom_style body command");

  if (0) bptr = NULL;

#define BODY_CLASS
#define BodyStyle(key,Class) \
  else if (strcmp(arg[0],#key) == 0) bptr = new Class(lmp,narg,arg);
#include "style_body.h"
#undef BodyStyle
#undef BODY_CLASS

  else error->all(FLERR,"Invalid body style");

  bptr->avec = this;
  icp = bptr->icp;
  dcp = bptr->dcp;

  // max size of forward/border comm
  // 7,16 are packed in pack_comm/pack_border
  // bptr values = max number of additional ivalues/dvalues from Body class

  size_forward = 7 + bptr->size_forward;
  size_border = 16 + bptr->size_border;

  // make copy of args if called externally, so can write to restart file
  // make no copy of args if called from read_restart()

  if (copyflag) {
    nargcopy = narg;
    argcopy = new char*[nargcopy];
    for (int i = 0; i < nargcopy; i++) {
      int n = strlen(arg[i]) + 1;
      argcopy[i] = new char[n];
      strcpy(argcopy[i],arg[i]);
    }
  }
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecBody::grow(int n)
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
  body = memory->grow(atom->body,nmax,"atom:body");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecBody::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;
  rmass = atom->rmass; angmom = atom->angmom; torque = atom->torque;
  body = atom->body;
}

/* ----------------------------------------------------------------------
   grow bonus data structure
------------------------------------------------------------------------- */

void AtomVecBody::grow_bonus()
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

void AtomVecBody::copy(int i, int j, int delflag)
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

  if (delflag && body[j] >= 0) {
    icp->put(bonus[body[j]].iindex);
    dcp->put(bonus[body[j]].dindex);
    copy_bonus(nlocal_bonus-1,body[j]);
    nlocal_bonus--;
  }

  // if atom I has bonus data, reset I's bonus.ilocal to loc J
  // do NOT do this if self-copy (I=J) since I's bonus data is already deleted

  if (body[i] >= 0 && i != j) bonus[body[i]].ilocal = j;
  body[j] = body[i];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ----------------------------------------------------------------------
   copy bonus data from I to J, effectively deleting the J entry
   also reset body that points to I to now point to J
------------------------------------------------------------------------- */

void AtomVecBody::copy_bonus(int i, int j)
{
  body[bonus[i].ilocal] = j;
  memcpy(&bonus[j],&bonus[i],sizeof(Bonus));
}

/* ----------------------------------------------------------------------
   clear ghost info in bonus data
   called before ghosts are recommunicated in comm and irregular
------------------------------------------------------------------------- */

void AtomVecBody::clear_bonus()
{
  int nall = nlocal_bonus + nghost_bonus;
  for (int i = nlocal_bonus; i < nall; i++) {
    icp->put(bonus[i].iindex);
    dcp->put(bonus[i].dindex);
  }
  nghost_bonus = 0;
}

/* ---------------------------------------------------------------------- */

int AtomVecBody::pack_comm(int n, int *list, double *buf,
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
      if (body[j] >= 0) {
        quat = bonus[body[j]].quat;
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
        m += bptr->pack_comm_body(&bonus[body[j]],&buf[m]);
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
      if (body[j] >= 0) {
        quat = bonus[body[j]].quat;
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
        m += bptr->pack_comm_body(&bonus[body[j]],&buf[m]);
      }
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBody::pack_comm_vel(int n, int *list, double *buf,
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
      if (body[j] >= 0) {
        quat = bonus[body[j]].quat;
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
        m += bptr->pack_comm_body(&bonus[body[j]],&buf[m]);
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
        if (body[j] >= 0) {
          quat = bonus[body[j]].quat;
          buf[m++] = quat[0];
          buf[m++] = quat[1];
          buf[m++] = quat[2];
          buf[m++] = quat[3];
          m += bptr->pack_comm_body(&bonus[body[j]],&buf[m]);
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
        if (body[j] >= 0) {
          quat = bonus[body[j]].quat;
          buf[m++] = quat[0];
          buf[m++] = quat[1];
          buf[m++] = quat[2];
          buf[m++] = quat[3];
          m += bptr->pack_comm_body(&bonus[body[j]],&buf[m]);
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

int AtomVecBody::pack_comm_hybrid(int n, int *list, double *buf)
{
  int i,j,m;
  double *quat;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (body[j] >= 0) {
      quat = bonus[body[j]].quat;
      buf[m++] = quat[0];
      buf[m++] = quat[1];
      buf[m++] = quat[2];
      buf[m++] = quat[3];
      m += bptr->pack_comm_body(&bonus[body[j]],&buf[m]);
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecBody::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;
  double *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    if (body[i] >= 0) {
      quat = bonus[body[i]].quat;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
      m += bptr->unpack_comm_body(&bonus[body[i]],&buf[m]);
    }
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecBody::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;
  double *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    if (body[i] >= 0) {
      quat = bonus[body[i]].quat;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
      m += bptr->unpack_comm_body(&bonus[body[i]],&buf[m]);
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

int AtomVecBody::unpack_comm_hybrid(int n, int first, double *buf)
{
  int i,m,last;
  double *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    if (body[i] >= 0) {
      quat = bonus[body[i]].quat;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
      m += bptr->unpack_comm_body(&bonus[body[i]],&buf[m]);
    }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBody::pack_reverse(int n, int first, double *buf)
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

int AtomVecBody::pack_reverse_hybrid(int n, int first, double *buf)
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

void AtomVecBody::unpack_reverse(int n, int *list, double *buf)
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

int AtomVecBody::unpack_reverse_hybrid(int n, int *list, double *buf)
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

int AtomVecBody::pack_border(int n, int *list, double *buf,
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
      if (body[j] < 0) buf[m++] = ubuf(0).d;
      else {
        buf[m++] = ubuf(1).d;
        quat = bonus[body[j]].quat;
        inertia = bonus[body[j]].inertia;
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
        buf[m++] = inertia[0];
        buf[m++] = inertia[1];
        buf[m++] = inertia[2];
        buf[m++] = ubuf(bonus[body[j]].ninteger).d;
        buf[m++] = ubuf(bonus[body[j]].ndouble).d;
        m += bptr->pack_border_body(&bonus[body[j]],&buf[m]);
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
      if (body[j] < 0) buf[m++] = ubuf(0).d;
      else {
        buf[m++] = ubuf(1).d;
        quat = bonus[body[j]].quat;
        inertia = bonus[body[j]].inertia;
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
        buf[m++] = inertia[0];
        buf[m++] = inertia[1];
        buf[m++] = inertia[2];
        buf[m++] = ubuf(bonus[body[j]].ninteger).d;
        buf[m++] = ubuf(bonus[body[j]].ndouble).d;
        m += bptr->pack_border_body(&bonus[body[j]],&buf[m]);
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBody::pack_border_vel(int n, int *list, double *buf,
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
      if (body[j] < 0) buf[m++] = ubuf(0).d;
      else {
        buf[m++] = ubuf(1).d;
        quat = bonus[body[j]].quat;
        inertia = bonus[body[j]].inertia;
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
        buf[m++] = inertia[0];
        buf[m++] = inertia[1];
        buf[m++] = inertia[2];
        buf[m++] = ubuf(bonus[body[j]].ninteger).d;
        buf[m++] = ubuf(bonus[body[j]].ndouble).d;
        m += bptr->pack_border_body(&bonus[body[j]],&buf[m]);
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
        if (body[j] < 0) buf[m++] = ubuf(0).d;
        else {
          buf[m++] = ubuf(1).d;
          quat = bonus[body[j]].quat;
          inertia = bonus[body[j]].inertia;
          buf[m++] = quat[0];
          buf[m++] = quat[1];
          buf[m++] = quat[2];
          buf[m++] = quat[3];
          buf[m++] = inertia[0];
          buf[m++] = inertia[1];
          buf[m++] = inertia[2];
          buf[m++] = ubuf(bonus[body[j]].ninteger).d;
          buf[m++] = ubuf(bonus[body[j]].ndouble).d;
          m += bptr->pack_border_body(&bonus[body[j]],&buf[m]);
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
        if (body[j] < 0) buf[m++] = ubuf(0).d;
        else {
          buf[m++] = ubuf(1).d;
          quat = bonus[body[j]].quat;
          inertia = bonus[body[j]].inertia;
          buf[m++] = quat[0];
          buf[m++] = quat[1];
          buf[m++] = quat[2];
          buf[m++] = quat[3];
          buf[m++] = inertia[0];
          buf[m++] = inertia[1];
          buf[m++] = inertia[2];
          buf[m++] = ubuf(bonus[body[j]].ninteger).d;
          buf[m++] = ubuf(bonus[body[j]].ndouble).d;
          m += bptr->pack_border_body(&bonus[body[j]],&buf[m]);
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

int AtomVecBody::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;
  double *quat,*c1,*c2,*c3,*inertia;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (body[j] < 0) buf[m++] = ubuf(0).d;
    else {
      buf[m++] = ubuf(1).d;
      quat = bonus[body[j]].quat;
      inertia = bonus[body[j]].inertia;
      buf[m++] = quat[0];
      buf[m++] = quat[1];
      buf[m++] = quat[2];
      buf[m++] = quat[3];
      buf[m++] = inertia[0];
      buf[m++] = inertia[1];
      buf[m++] = inertia[2];
      buf[m++] = ubuf(bonus[body[j]].ninteger).d;
      buf[m++] = ubuf(bonus[body[j]].ndouble).d;
      m += bptr->pack_border_body(&bonus[body[j]],&buf[m]);
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecBody::unpack_border(int n, int first, double *buf)
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
    body[i] = (int) ubuf(buf[m++]).i;
    if (body[i] == 0) body[i] = -1;
    else {
      j = nlocal_bonus + nghost_bonus;
      if (j == nmax_bonus) grow_bonus();
      quat = bonus[j].quat;
      inertia = bonus[j].inertia;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
      inertia[0] = buf[m++];
      inertia[1] = buf[m++];
      inertia[2] = buf[m++];
      bonus[j].ninteger = (int) ubuf(buf[m++]).i;
      bonus[j].ndouble = (int) ubuf(buf[m++]).i;
      bonus[j].ivalue = icp->get(bonus[j].ninteger,bonus[j].iindex);
      bonus[j].dvalue = dcp->get(bonus[j].ndouble,bonus[j].dindex);
      m += bptr->unpack_border_body(&bonus[j],&buf[m]);
      bonus[j].ilocal = i;
      body[i] = j;
      nghost_bonus++;
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecBody::unpack_border_vel(int n, int first, double *buf)
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
    body[i] = (int) ubuf(buf[m++]).i;
    if (body[i] == 0) body[i] = -1;
    else {
      j = nlocal_bonus + nghost_bonus;
      if (j == nmax_bonus) grow_bonus();
      quat = bonus[j].quat;
      inertia = bonus[j].inertia;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
      inertia[0] = buf[m++];
      inertia[1] = buf[m++];
      inertia[2] = buf[m++];
      bonus[j].ninteger = (int) ubuf(buf[m++]).i;
      bonus[j].ndouble = (int) ubuf(buf[m++]).i;
      bonus[j].ivalue = icp->get(bonus[j].ninteger,bonus[j].iindex);
      bonus[j].dvalue = dcp->get(bonus[j].ndouble,bonus[j].dindex);
      m += bptr->unpack_border_body(&bonus[j],&buf[m]);
      bonus[j].ilocal = i;
      body[i] = j;
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

int AtomVecBody::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,j,m,last;
  double *quat,*c1,*c2,*c3,*inertia;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    body[i] = (int) ubuf(buf[m++]).i;
    if (body[i] == 0) body[i] = -1;
    else {
      j = nlocal_bonus + nghost_bonus;
      if (j == nmax_bonus) grow_bonus();
      quat = bonus[j].quat;
      inertia = bonus[j].inertia;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
      inertia[0] = buf[m++];
      inertia[1] = buf[m++];
      inertia[2] = buf[m++];
      bonus[j].ninteger = (int) ubuf(buf[m++]).i;
      bonus[j].ndouble = (int) ubuf(buf[m++]).i;
      bonus[j].ivalue = icp->get(bonus[j].ninteger,bonus[j].iindex);
      bonus[j].dvalue = dcp->get(bonus[j].ndouble,bonus[j].dindex);
      m += bptr->unpack_border_body(&bonus[j],&buf[m]);
      bonus[j].ilocal = i;
      body[i] = j;
      nghost_bonus++;
    }
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecBody::pack_exchange(int i, double *buf)
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

  if (body[i] < 0) buf[m++] = ubuf(0).d;
  else {
    buf[m++] = ubuf(1).d;
    int j = body[i];
    double *quat = bonus[j].quat;
    double *inertia = bonus[j].inertia;
    buf[m++] = quat[0];
    buf[m++] = quat[1];
    buf[m++] = quat[2];
    buf[m++] = quat[3];
    buf[m++] = inertia[0];
    buf[m++] = inertia[1];
    buf[m++] = inertia[2];
    buf[m++] = ubuf(bonus[j].ninteger).d;
    buf[m++] = ubuf(bonus[j].ndouble).d;
    memcpy(&buf[m],bonus[j].ivalue,bonus[j].ninteger*sizeof(int));
    if (intdoubleratio == 1) m += bonus[j].ninteger;
    else m += (bonus[j].ninteger+1)/2;
    memcpy(&buf[m],bonus[j].dvalue,bonus[j].ndouble*sizeof(double));
    m += bonus[j].ndouble;
  }

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBody::unpack_exchange(double *buf)
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
  angmom[nlocal][0] = buf[m++];
  angmom[nlocal][1] = buf[m++];
  angmom[nlocal][2] = buf[m++];

  body[nlocal] = (int) ubuf(buf[m++]).i;
  if (body[nlocal] == 0) body[nlocal] = -1;
  else {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    double *quat = bonus[nlocal_bonus].quat;
    double *inertia = bonus[nlocal_bonus].inertia;
    quat[0] = buf[m++];
    quat[1] = buf[m++];
    quat[2] = buf[m++];
    quat[3] = buf[m++];
    inertia[0] = buf[m++];
    inertia[1] = buf[m++];
    inertia[2] = buf[m++];
    bonus[nlocal_bonus].ninteger = (int) ubuf(buf[m++]).i;
    bonus[nlocal_bonus].ndouble = (int) ubuf(buf[m++]).i;
    bonus[nlocal_bonus].ivalue = icp->get(bonus[nlocal_bonus].ninteger,
					  bonus[nlocal_bonus].iindex);
    bonus[nlocal_bonus].dvalue = dcp->get(bonus[nlocal_bonus].ndouble,
					  bonus[nlocal_bonus].dindex);
    memcpy(bonus[nlocal_bonus].ivalue,&buf[m],
           bonus[nlocal_bonus].ninteger*sizeof(int));
    if (intdoubleratio == 1) m += bonus[nlocal_bonus].ninteger;
    else m += (bonus[nlocal_bonus].ninteger+1)/2;
    memcpy(bonus[nlocal_bonus].dvalue,&buf[m],
           bonus[nlocal_bonus].ndouble*sizeof(double));
    m += bonus[nlocal_bonus].ndouble;

    bonus[nlocal_bonus].ilocal = nlocal;
    body[nlocal] = nlocal_bonus++;
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

int AtomVecBody::size_restart()
{
  int i;

  int n = 0;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++)
    if (body[i] >= 0) {
      n += 25;
      if (intdoubleratio == 1) n += bonus[body[i]].ninteger;
      else n += (bonus[body[i]].ninteger+1)/2;
      n += bonus[body[i]].ndouble;
    } else n += 16;

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

int AtomVecBody::pack_restart(int i, double *buf)
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

  if (body[i] < 0) buf[m++] = ubuf(0).d;
  else {
    buf[m++] = ubuf(1).d;
    int j = body[i];
    double *quat = bonus[j].quat;
    double *inertia = bonus[j].inertia;
    buf[m++] = quat[0];
    buf[m++] = quat[1];
    buf[m++] = quat[2];
    buf[m++] = quat[3];
    buf[m++] = inertia[0];
    buf[m++] = inertia[1];
    buf[m++] = inertia[2];
    buf[m++] = ubuf(bonus[j].ninteger).d;
    buf[m++] = ubuf(bonus[j].ndouble).d;
    memcpy(&buf[m],bonus[j].ivalue,bonus[j].ninteger*sizeof(int));
    if (intdoubleratio == 1) m += bonus[j].ninteger;
    else m += (bonus[j].ninteger+1)/2;
    memcpy(&buf[m],bonus[j].dvalue,bonus[j].ndouble*sizeof(double));
    m += bonus[j].ndouble;
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

int AtomVecBody::unpack_restart(double *buf)
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

  rmass[nlocal] = buf[m++];
  angmom[nlocal][0] = buf[m++];
  angmom[nlocal][1] = buf[m++];
  angmom[nlocal][2] = buf[m++];

  body[nlocal] = (int) ubuf(buf[m++]).i;
  if (body[nlocal] == 0) body[nlocal] = -1;
  else {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    double *quat = bonus[nlocal_bonus].quat;
    double *inertia = bonus[nlocal_bonus].inertia;
    quat[0] = buf[m++];
    quat[1] = buf[m++];
    quat[2] = buf[m++];
    quat[3] = buf[m++];
    inertia[0] = buf[m++];
    inertia[1] = buf[m++];
    inertia[2] = buf[m++];
    bonus[nlocal_bonus].ninteger = (int) ubuf(buf[m++]).i;
    bonus[nlocal_bonus].ndouble = (int) ubuf(buf[m++]).i;
    bonus[nlocal_bonus].ivalue = icp->get(bonus[nlocal_bonus].ninteger,
					  bonus[nlocal_bonus].iindex);
    bonus[nlocal_bonus].dvalue = dcp->get(bonus[nlocal_bonus].ndouble,
					  bonus[nlocal_bonus].dindex);
    memcpy(bonus[nlocal_bonus].ivalue,&buf[m],
           bonus[nlocal_bonus].ninteger*sizeof(int));
    if (intdoubleratio == 1) m += bonus[nlocal_bonus].ninteger;
    else m += (bonus[nlocal_bonus].ninteger+1)/2;
    memcpy(bonus[nlocal_bonus].dvalue,&buf[m],
           bonus[nlocal_bonus].ndouble*sizeof(double));
    m += bonus[nlocal_bonus].ndouble;
    bonus[nlocal_bonus].ilocal = nlocal;
    body[nlocal] = nlocal_bonus++;
  }

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecBody::write_restart_settings(FILE *fp)
{
  fwrite(&nargcopy,sizeof(int),1,fp);
  for (int i = 0; i < nargcopy; i++) {
    int n = strlen(argcopy[i]) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(argcopy[i],sizeof(char),n,fp);
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecBody::read_restart_settings(FILE *fp)
{
  int n;

  int me = comm->me;
  if (me == 0) fread(&nargcopy,sizeof(int),1,fp);
  MPI_Bcast(&nargcopy,1,MPI_INT,0,world);
  argcopy = new char*[nargcopy];
    
  for (int i = 0; i < nargcopy; i++) {
    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    argcopy[i] = new char[n];
    if (me == 0) fread(argcopy[i],sizeof(char),n,fp);
    MPI_Bcast(argcopy[i],n,MPI_CHAR,0,world);
  }

  copyflag = 0;
  settings(nargcopy,argcopy);
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecBody::create_atom(int itype, double *coord)
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
  body[nlocal] = -1;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecBody::data_atom(double *coord, imageint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = ATOTAGINT(values[0]);
  type[nlocal] = atoi(values[1]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  body[nlocal] = atoi(values[2]);
  if (body[nlocal] == 0) body[nlocal] = -1;
  else if (body[nlocal] == 1) body[nlocal] = 0;
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

int AtomVecBody::data_atom_hybrid(int nlocal, char **values)
{
  body[nlocal] = atoi(values[0]);
  if (body[nlocal] == 0) body[nlocal] = -1;
  else if (body[nlocal] == 1) body[nlocal] = 0;
  else error->one(FLERR,"Invalid atom type in Atoms section of data file");

  rmass[nlocal] = atof(values[1]);
  if (rmass[nlocal] <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  return 2;
}

/* ----------------------------------------------------------------------
   unpack one body from Bodies section of data file
------------------------------------------------------------------------- */

void AtomVecBody::data_body(int m, int ninteger, int ndouble, 
                             char **ivalues, char **dvalues)
{
  if (body[m]) error->one(FLERR,"Assigning body parameters to non-body atom");
  if (nlocal_bonus == nmax_bonus) grow_bonus();
  bptr->data_body(nlocal_bonus,ninteger,ndouble,ivalues,dvalues);
  bonus[nlocal_bonus].ilocal = m;
  body[m] = nlocal_bonus++;
}

/* ----------------------------------------------------------------------
   unpack one tri from Velocities section of data file
------------------------------------------------------------------------- */

void AtomVecBody::data_vel(int m, char **values)
{
  v[m][0] = atof(values[0]);
  v[m][1] = atof(values[1]);
  v[m][2] = atof(values[2]);
  angmom[m][0] = atof(values[3]);
  angmom[m][1] = atof(values[4]);
  angmom[m][2] = atof(values[5]);
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one body in Velocities section of data file
------------------------------------------------------------------------- */

int AtomVecBody::data_vel_hybrid(int m, char **values)
{
  angmom[m][0] = atof(values[0]);
  angmom[m][1] = atof(values[1]);
  angmom[m][2] = atof(values[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecBody::pack_data(double **buf)
{
  double c2mc1[2],c3mc1[3],norm[3];
  double area;

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = ubuf(type[i]).d;
    if (body[i] < 0) buf[i][2] = ubuf(0).d;
    else buf[i][2] = ubuf(1).d;
    buf[i][3] = rmass[i];
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

int AtomVecBody::pack_data_hybrid(int i, double *buf)
{
  if (body[i] < 0) buf[0] = ubuf(0).d;
  else buf[0] = ubuf(1).d;
  buf[1] = rmass[i];
  return 2;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecBody::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT " %d %d %g %g %g %g %d %d %d\n",
            (tagint) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,
            (int) ubuf(buf[i][2]).i,
            buf[i][3],buf[i][4],buf[i][5],buf[i][6],
            (int) ubuf(buf[i][7]).i,(int) ubuf(buf[i][8]).i,
            (int) ubuf(buf[i][9]).i);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecBody::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %d %g",(int) ubuf(buf[0]).i,buf[1]);
  return 2;
}

/* ----------------------------------------------------------------------
   pack velocity info for data file
------------------------------------------------------------------------- */

void AtomVecBody::pack_vel(double **buf)
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

int AtomVecBody::pack_vel_hybrid(int i, double *buf)
{
  buf[0] = angmom[i][0];
  buf[1] = angmom[i][1];
  buf[2] = angmom[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
   write velocity info to data file
------------------------------------------------------------------------- */

void AtomVecBody::write_vel(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT " %g %g %g %g %g %g\n",
            (tagint) ubuf(buf[i][0]).i,buf[i][1],buf[i][2],buf[i][3],
            buf[i][4],buf[i][5],buf[i][6]);
}

/* ----------------------------------------------------------------------
   write hybrid velocity info to data file
------------------------------------------------------------------------- */

int AtomVecBody::write_vel_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %g %g %g",buf[0],buf[1],buf[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecBody::memory_usage()
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
  if (atom->memcheck("torque")) bytes += 
                                  memory->usage(torque,nmax*comm->nthreads,3);
  if (atom->memcheck("body")) bytes += memory->usage(body,nmax);

  bytes += nmax_bonus*sizeof(Bonus);
  bytes += icp->size + dcp->size;

  int nall = nlocal_bonus + nghost_bonus;
  for (int i = 0; i < nall; i++) {
    bytes += bonus[i].ninteger * sizeof(int);
    bytes += bonus[i].ndouble * sizeof(double);
  }

  return bytes;
}

/* ----------------------------------------------------------------------
   debug method for sanity checking of own/bonus data pointers
------------------------------------------------------------------------- */

/*
void AtomVecBody::check(int flag)
{
  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->body[i] >= 0 && atom->body[i] >= nlocal_bonus) {
      printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
      errorx->one(FLERR,"BAD AAA");
    }
  }
  for (int i = atom->nlocal; i < atom->nlocal+atom->nghost; i++) {
    if (atom->body[i] >= 0 && 
        (atom->body[i] < nlocal_bonus || 
         atom->body[i] >= nlocal_bonus+nghost_bonus)) {
      printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
      errorx->one(FLERR,"BAD BBB");
    }
  }
  for (int i = 0; i < nlocal_bonus; i++) {
    if (bonus[i].ilocal < 0 || bonus[i].ilocal >= atom->nlocal) {
      printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
      errorx->one(FLERR,"BAD CCC");
    }
  }
  for (int i = 0; i < nlocal_bonus; i++) {
    if (atom->body[bonus[i].ilocal] != i) {
      printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
      errorx->one(FLERR,"BAD DDD");
    }
  }
  for (int i = nlocal_bonus; i < nlocal_bonus+nghost_bonus; i++) {
    if (bonus[i].ilocal < atom->nlocal || 
        bonus[i].ilocal >= atom->nlocal+atom->nghost) {
      printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
      errorx->one(FLERR,"BAD EEE");
    }
  }
  for (int i = nlocal_bonus; i < nlocal_bonus+nghost_bonus; i++) {
    if (atom->body[bonus[i].ilocal] != i) {
      printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
      errorx->one(FLERR,"BAD FFF");
    }
  }
}
*/
