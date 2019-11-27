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

#include "atom_vec_ellipsoid.h"
#include <cstring>
#include "math_extra.h"
#include "atom.h"
#include "modify.h"
#include "fix.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

AtomVecEllipsoid::AtomVecEllipsoid(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;
  bonus_flag = 1;

  size_forward_bonus = 4;
  size_border_bonus = 8;
  size_restart_bonus_one = 7;
  size_data_bonus = 8;

  atom->ellipsoid_flag = 1;
  atom->rmass_flag = atom->angmom_flag = atom->torque_flag = 1;

  nlocal_bonus = nghost_bonus = nmax_bonus = 0;
  bonus = NULL;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in the string does not matter
  //   except fields_data_atom and fields_data_vel which must match data file

  fields_grow = (char *) "rmass angmom torque ellipsoid";
  fields_copy = (char *) "rmass angmom";
  fields_comm = NULL;
  fields_comm_vel = (char *) "angmom";
  fields_reverse = (char *) "torque";
  fields_border = (char *) "rmass";
  fields_border_vel = (char *) "rmass angmom";
  fields_exchange = (char *) "rmass angmom";
  fields_restart = (char *) "rmass angmom";
  fields_create = (char *) "rmass angmom ellipsoid";
  fields_data_atom = (char *) "id type ellipsoid rmass x";
  fields_data_vel = (char *) "angmom";

  setup_fields();
}

/* ---------------------------------------------------------------------- */

AtomVecEllipsoid::~AtomVecEllipsoid()
{
  memory->sfree(bonus);
}

/* ----------------------------------------------------------------------
   grow bonus data structure
------------------------------------------------------------------------- */

void AtomVecEllipsoid::grow_bonus()
{
  nmax_bonus = grow_nmax_bonus(nmax_bonus);
  if (nmax_bonus < 0)
    error->one(FLERR,"Per-processor system is too big");

  bonus = (Bonus *) memory->srealloc(bonus,nmax_bonus*sizeof(Bonus),
                                     "atom:bonus");
}

/* ----------------------------------------------------------------------
   copy atom I bonus info to atom J
------------------------------------------------------------------------- */

void AtomVecEllipsoid::copy_bonus(int i, int j, int delflag)
{
  int *ellipsoid = atom->ellipsoid;

  // if deleting atom J via delflag and J has bonus data, then delete it

  if (delflag && ellipsoid[j] >= 0) {
    copy_bonus_all(nlocal_bonus-1,ellipsoid[j]);
    nlocal_bonus--;
  }

  // if atom I has bonus data, reset I's bonus.ilocal to loc J
  // do NOT do this if self-copy (I=J) since I's bonus data is already deleted

  if (ellipsoid[i] >= 0 && i != j) bonus[ellipsoid[i]].ilocal = j;
  ellipsoid[j] = ellipsoid[i];
}

/* ----------------------------------------------------------------------
   copy bonus data from I to J, effectively deleting the J entry
   also reset ellipsoid that points to I to now point to J
------------------------------------------------------------------------- */

void AtomVecEllipsoid::copy_bonus_all(int i, int j)
{
  atom->ellipsoid[bonus[i].ilocal] = j;
  memcpy(&bonus[j],&bonus[i],sizeof(Bonus));
}

/* ----------------------------------------------------------------------
   clear ghost info in bonus data
   called before ghosts are recommunicated in comm and irregular
------------------------------------------------------------------------- */

void AtomVecEllipsoid::clear_bonus()
{
  nghost_bonus = 0;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->clear_bonus();
}

/* ---------------------------------------------------------------------- */

int AtomVecEllipsoid::pack_comm_bonus(int n, int *list, double *buf)
{
  int i,j,m;
  double *quat;

  int *ellipsoid = atom->ellipsoid;

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

void AtomVecEllipsoid::unpack_comm_bonus(int n, int first, double *buf)
{
  int i,m,last;
  double *quat;

  int *ellipsoid = atom->ellipsoid;

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
}

/* ---------------------------------------------------------------------- */

int AtomVecEllipsoid::pack_border_bonus(int n, int *list, double *buf)
{
  int i,j,m;
  double dx,dy,dz;
  double *shape,*quat;

  int *ellipsoid = atom->ellipsoid;

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

int AtomVecEllipsoid::unpack_border_bonus(int n, int first, double *buf)
{
  int i,j,m,last;
  double *shape,*quat;

  int *ellipsoid = atom->ellipsoid;

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

int AtomVecEllipsoid::pack_exchange_bonus(int i, double *buf)
{
  int m = 0;

  int *ellipsoid = atom->ellipsoid;

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

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecEllipsoid::unpack_exchange_bonus(int ilocal, double *buf)
{
  int m = 0;

  int *ellipsoid = atom->ellipsoid;

  ellipsoid[ilocal] = (int) ubuf(buf[m++]).i;
  if (ellipsoid[ilocal] == 0) ellipsoid[ilocal] = -1;
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
    bonus[nlocal_bonus].ilocal = ilocal;
    ellipsoid[ilocal] = nlocal_bonus++;
  }

  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecEllipsoid::size_restart_bonus()
{
  int i;

  int *ellipsoid = atom->ellipsoid;

  int n = 0;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) {
    if (ellipsoid[i] >= 0) n += size_restart_bonus_one;
    n++;
  }

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including bonus data
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVecEllipsoid::pack_restart_bonus(int i, double *buf)
{
  int m = 0;

  int *ellipsoid = atom->ellipsoid;

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

  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including bonus data
------------------------------------------------------------------------- */

int AtomVecEllipsoid::unpack_restart_bonus(int ilocal, double *buf)
{
  int m = 0;

  int *ellipsoid = atom->ellipsoid;

  ellipsoid[ilocal] = (int) ubuf(buf[m++]).i;
  if (ellipsoid[ilocal] == 0) ellipsoid[ilocal] = -1;
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
    bonus[nlocal_bonus].ilocal = ilocal;
    ellipsoid[ilocal] = nlocal_bonus++;
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack one line from Ellipsoids section of data file
------------------------------------------------------------------------- */

void AtomVecEllipsoid::data_atom_bonus(int m, char **values)
{
  int *ellipsoid = atom->ellipsoid;

  if (ellipsoid[m])
    error->one(FLERR,"Assigning ellipsoid parameters to non-ellipsoid atom");

  if (nlocal_bonus == nmax_bonus) grow_bonus();

  double *shape = bonus[nlocal_bonus].shape;
  shape[0] = 0.5 * utils::numeric(FLERR,values[0],true,lmp);
  shape[1] = 0.5 * utils::numeric(FLERR,values[1],true,lmp);
  shape[2] = 0.5 * utils::numeric(FLERR,values[2],true,lmp);
  if (shape[0] <= 0.0 || shape[1] <= 0.0 || shape[2] <= 0.0)
    error->one(FLERR,"Invalid shape in Ellipsoids section of data file");

  double *quat = bonus[nlocal_bonus].quat;
  quat[0] = utils::numeric(FLERR,values[3],true,lmp);
  quat[1] = utils::numeric(FLERR,values[4],true,lmp);
  quat[2] = utils::numeric(FLERR,values[5],true,lmp);
  quat[3] = utils::numeric(FLERR,values[6],true,lmp);
  MathExtra::qnormalize(quat);

  // reset ellipsoid mass
  // previously stored density in rmass

  atom->rmass[m] *= 4.0*MY_PI/3.0 * shape[0]*shape[1]*shape[2];

  bonus[nlocal_bonus].ilocal = m;
  ellipsoid[m] = nlocal_bonus++;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated bonus memory
------------------------------------------------------------------------- */

bigint AtomVecEllipsoid::memory_usage_bonus()
{
  bigint bytes = 0;
  bytes += nmax_bonus*sizeof(Bonus);
  return bytes;
}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecEllipsoid::create_atom_post(int ilocal)
{
  atom->rmass[ilocal] = 1.0;
  atom->ellipsoid[ilocal] = -1;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecEllipsoid::data_atom_post(int ilocal)
{
  ellipsoid_flag = atom->ellipsoid[ilocal];
  if (ellipsoid_flag == 0) ellipsoid_flag = -1;
  else if (ellipsoid_flag == 1) ellipsoid_flag = 0;
  else error->one(FLERR,"Invalid ellipsoid flag in Atoms section of data file");
  atom->ellipsoid[ilocal] = ellipsoid_flag;

  if (atom->rmass[ilocal] <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");
}

/* ----------------------------------------------------------------------
   modify values for AtomVec::pack_data() to pack
------------------------------------------------------------------------- */

void AtomVecEllipsoid::pack_data_pre(int ilocal)
{ 
  double *shape;

  ellipsoid_flag = atom->ellipsoid[ilocal];
  rmass = atom->rmass[ilocal];

  if (ellipsoid_flag < 0) atom->ellipsoid[ilocal] = 0;
  else atom->ellipsoid[ilocal] = 1;

  if (ellipsoid_flag >= 0) {
    shape = bonus[ellipsoid_flag].shape;
    atom->rmass[ilocal] /= 4.0*MY_PI/3.0 * shape[0]*shape[1]*shape[2];
  }
}

/* ----------------------------------------------------------------------
   unmodify values packed by AtomVec::pack_data()
------------------------------------------------------------------------- */

void AtomVecEllipsoid::pack_data_post(int ilocal)
{ 
  atom->ellipsoid[ilocal] = ellipsoid_flag;
  atom->rmass[ilocal] = rmass;
}

/* ----------------------------------------------------------------------
   set shape values in bonus data for particle I
   oriented aligned with xyz axes
   this may create or delete entry in bonus data
------------------------------------------------------------------------- */

void AtomVecEllipsoid::
set_shape(int i, double shapex, double shapey, double shapez)
{
  int *ellipsoid = atom->ellipsoid;

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
    copy_bonus_all(nlocal_bonus-1,ellipsoid[i]);
    nlocal_bonus--;
    ellipsoid[i] = -1;
  } else {
    double *shape = bonus[ellipsoid[i]].shape;
    shape[0] = shapex;
    shape[1] = shapey;
    shape[2] = shapez;
  }
}
