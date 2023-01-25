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
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "atom_vec_ellipsoid.h"

#include "atom.h"
#include "error.h"
#include "fix.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"

#include <cstring>

using namespace LAMMPS_NS;
using MathConst::MY_PI;

/* ---------------------------------------------------------------------- */

AtomVecEllipsoid::AtomVecEllipsoid(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::ATOMIC;
  bonus_flag = 1;

  size_forward_bonus = 4;
  size_border_bonus = 8;
  size_restart_bonus_one = 8;
  size_data_bonus = 8;

  atom->ellipsoid_flag = 1;
  atom->rmass_flag = atom->angmom_flag = atom->torque_flag = 1;

  nlocal_bonus = nghost_bonus = nmax_bonus = 0;
  bonus = nullptr;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = {"rmass", "angmom", "torque", "ellipsoid"};
  fields_copy = {"rmass", "angmom"};
  fields_comm_vel = {"angmom"};
  fields_reverse = {"torque"};
  fields_border = {"rmass"};
  fields_border_vel = {"rmass", "angmom"};
  fields_exchange = {"rmass", "angmom"};
  fields_restart = {"rmass", "angmom"};
  fields_create = {"rmass", "angmom", "ellipsoid"};
  fields_data_atom = {"id", "type", "ellipsoid", "rmass", "x"};
  fields_data_vel = {"id", "v", "angmom"};

  setup_fields();
}

/* ---------------------------------------------------------------------- */

AtomVecEllipsoid::~AtomVecEllipsoid()
{
  memory->sfree(bonus);
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecEllipsoid::grow_pointers()
{
  ellipsoid = atom->ellipsoid;
  rmass = atom->rmass;
  angmom = atom->angmom;
}

/* ----------------------------------------------------------------------
   grow bonus data structure
------------------------------------------------------------------------- */

void AtomVecEllipsoid::grow_bonus()
{
  nmax_bonus = grow_nmax_bonus(nmax_bonus);
  if (nmax_bonus < 0) error->one(FLERR, "Per-processor system is too big");

  bonus = (Bonus *) memory->srealloc(bonus, nmax_bonus * sizeof(Bonus), "atom:bonus");
}

/* ----------------------------------------------------------------------
   copy atom I bonus info to atom J
------------------------------------------------------------------------- */

void AtomVecEllipsoid::copy_bonus(int i, int j, int delflag)
{
  // if deleting atom J via delflag and J has bonus data, then delete it

  if (delflag && ellipsoid[j] >= 0) {
    copy_bonus_all(nlocal_bonus - 1, ellipsoid[j]);
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
  ellipsoid[bonus[i].ilocal] = j;
  memcpy(&bonus[j], &bonus[i], sizeof(Bonus));
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
  int i, j, m;
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

void AtomVecEllipsoid::unpack_comm_bonus(int n, int first, double *buf)
{
  int i, m, last;
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
}

/* ---------------------------------------------------------------------- */

int AtomVecEllipsoid::pack_border_bonus(int n, int *list, double *buf)
{
  int i, j, m;
  double *shape, *quat;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (ellipsoid[j] < 0)
      buf[m++] = ubuf(0).d;
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
  int i, j, m, last;
  double *shape, *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    ellipsoid[i] = (int) ubuf(buf[m++]).i;
    if (ellipsoid[i] == 0)
      ellipsoid[i] = -1;
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

  if (ellipsoid[i] < 0)
    buf[m++] = ubuf(0).d;
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

  ellipsoid[ilocal] = (int) ubuf(buf[m++]).i;
  if (ellipsoid[ilocal] == 0)
    ellipsoid[ilocal] = -1;
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

  int n = 0;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) {
    if (ellipsoid[i] >= 0)
      n += size_restart_bonus_one;
    else
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

  if (ellipsoid[i] < 0)
    buf[m++] = ubuf(0).d;
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

  ellipsoid[ilocal] = (int) ubuf(buf[m++]).i;
  if (ellipsoid[ilocal] == 0)
    ellipsoid[ilocal] = -1;
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

void AtomVecEllipsoid::data_atom_bonus(int m, const std::vector<std::string> &values)
{
  if (ellipsoid[m]) error->one(FLERR, "Assigning ellipsoid parameters to non-ellipsoid atom");

  if (nlocal_bonus == nmax_bonus) grow_bonus();

  double *shape = bonus[nlocal_bonus].shape;
  int ivalue = 1;
  shape[0] = 0.5 * utils::numeric(FLERR, values[ivalue++], true, lmp);
  shape[1] = 0.5 * utils::numeric(FLERR, values[ivalue++], true, lmp);
  shape[2] = 0.5 * utils::numeric(FLERR, values[ivalue++], true, lmp);
  if (shape[0] <= 0.0 || shape[1] <= 0.0 || shape[2] <= 0.0)
    error->one(FLERR, "Invalid shape in Ellipsoids section of data file");

  double *quat = bonus[nlocal_bonus].quat;
  quat[0] = utils::numeric(FLERR, values[ivalue++], true, lmp);
  quat[1] = utils::numeric(FLERR, values[ivalue++], true, lmp);
  quat[2] = utils::numeric(FLERR, values[ivalue++], true, lmp);
  quat[3] = utils::numeric(FLERR, values[ivalue++], true, lmp);
  MathExtra::qnormalize(quat);

  // reset ellipsoid mass
  // previously stored density in rmass

  rmass[m] *= 4.0 * MY_PI / 3.0 * shape[0] * shape[1] * shape[2];

  bonus[nlocal_bonus].ilocal = m;
  ellipsoid[m] = nlocal_bonus++;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated bonus memory
------------------------------------------------------------------------- */

double AtomVecEllipsoid::memory_usage_bonus()
{
  double bytes = 0;
  bytes += nmax_bonus * sizeof(Bonus);
  return bytes;
}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecEllipsoid::create_atom_post(int ilocal)
{
  rmass[ilocal] = 1.0;
  ellipsoid[ilocal] = -1;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecEllipsoid::data_atom_post(int ilocal)
{
  ellipsoid_flag = ellipsoid[ilocal];
  if (ellipsoid_flag == 0)
    ellipsoid_flag = -1;
  else if (ellipsoid_flag == 1)
    ellipsoid_flag = 0;
  else
    error->one(FLERR, "Invalid ellipsoid flag in Atoms section of data file");
  ellipsoid[ilocal] = ellipsoid_flag;

  if (rmass[ilocal] <= 0.0) error->one(FLERR, "Invalid density in Atoms section of data file");

  angmom[ilocal][0] = 0.0;
  angmom[ilocal][1] = 0.0;
  angmom[ilocal][2] = 0.0;
}

/* ----------------------------------------------------------------------
   modify values for AtomVec::pack_data() to pack
------------------------------------------------------------------------- */

void AtomVecEllipsoid::pack_data_pre(int ilocal)
{
  double *shape;

  ellipsoid_flag = atom->ellipsoid[ilocal];
  rmass_one = atom->rmass[ilocal];

  if (ellipsoid_flag < 0)
    ellipsoid[ilocal] = 0;
  else
    ellipsoid[ilocal] = 1;

  if (ellipsoid_flag >= 0) {
    shape = bonus[ellipsoid_flag].shape;
    rmass[ilocal] /= 4.0 * MY_PI / 3.0 * shape[0] * shape[1] * shape[2];
  }
}

/* ----------------------------------------------------------------------
   unmodify values packed by AtomVec::pack_data()
------------------------------------------------------------------------- */

void AtomVecEllipsoid::pack_data_post(int ilocal)
{
  ellipsoid[ilocal] = ellipsoid_flag;
  rmass[ilocal] = rmass_one;
}

/* ----------------------------------------------------------------------
   pack bonus ellipsoid info for writing to data file
   if buf is nullptr, just return buffer size
------------------------------------------------------------------------- */

int AtomVecEllipsoid::pack_data_bonus(double *buf, int /*flag*/)
{
  int i, j;

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int m = 0;
  for (i = 0; i < nlocal; i++) {
    if (ellipsoid[i] < 0) continue;
    if (buf) {
      buf[m++] = ubuf(tag[i]).d;
      j = ellipsoid[i];
      buf[m++] = 2.0 * bonus[j].shape[0];
      buf[m++] = 2.0 * bonus[j].shape[1];
      buf[m++] = 2.0 * bonus[j].shape[2];
      buf[m++] = bonus[j].quat[0];
      buf[m++] = bonus[j].quat[1];
      buf[m++] = bonus[j].quat[2];
      buf[m++] = bonus[j].quat[3];
    } else
      m += size_data_bonus;
  }

  return m;
}

/* ----------------------------------------------------------------------
   write bonus ellipsoid info to data file
------------------------------------------------------------------------- */

void AtomVecEllipsoid::write_data_bonus(FILE *fp, int n, double *buf, int /*flag*/)
{
  int i = 0;
  while (i < n) {
    fmt::print(fp, "{} {} {} {} {} {} {} {}\n", ubuf(buf[i]).i, buf[i + 1], buf[i + 2], buf[i + 3],
               buf[i + 4], buf[i + 5], buf[i + 6], buf[i + 7]);
    i += size_data_bonus;
  }
}

/* ----------------------------------------------------------------------
   set shape values in bonus data for particle I
   oriented aligned with xyz axes
   this may create or delete entry in bonus data
------------------------------------------------------------------------- */

void AtomVecEllipsoid::set_shape(int i, double shapex, double shapey, double shapez)
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
    copy_bonus_all(nlocal_bonus - 1, ellipsoid[i]);
    nlocal_bonus--;
    ellipsoid[i] = -1;
  } else {
    double *shape = bonus[ellipsoid[i]].shape;
    shape[0] = shapex;
    shape[1] = shapey;
    shape[2] = shapez;
  }
}
