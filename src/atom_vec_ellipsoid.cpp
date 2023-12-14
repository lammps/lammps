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
#include <tr1/cmath>

using namespace LAMMPS_NS;
using MathConst::MY_PI;

#define EPSBLOCK2 1.0e-3

/* ---------------------------------------------------------------------- */

AtomVecEllipsoid::AtomVecEllipsoid(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::ATOMIC;
  bonus_flag = 1;

  size_forward_bonus = 4;
  size_border_bonus = 10;
  size_restart_bonus_one = 10;
  size_data_bonus = 10;

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
  double *block;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (ellipsoid[j] < 0)
      buf[m++] = ubuf(0).d;
    else {
      buf[m++] = ubuf(1).d;
      shape = bonus[ellipsoid[j]].shape;
      quat = bonus[ellipsoid[j]].quat;
      block = bonus[ellipsoid[j]].block;
      buf[m++] = shape[0];
      buf[m++] = shape[1];
      buf[m++] = shape[2];
      buf[m++] = quat[0];
      buf[m++] = quat[1];
      buf[m++] = quat[2];
      buf[m++] = quat[3];
      buf[m++] = block[0];
      buf[m++] = block[1];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecEllipsoid::unpack_border_bonus(int n, int first, double *buf)
{
  int i, j, m, last;
  double *shape, *quat;
  double *block;
  bool flag_super;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (ubuf(buf[m++]).i == 0)
      ellipsoid[i] = -1;
    else {
      j = nlocal_bonus + nghost_bonus;
      if (j == nmax_bonus) grow_bonus();
      shape = bonus[j].shape;
      quat = bonus[j].quat;
      block = bonus[j].block;
      shape[0] = buf[m++];
      shape[1] = buf[m++];
      shape[2] = buf[m++];
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
      block[0] = buf[m++];
      block[1] = buf[m++];
      block[2] = block[0] / block[1];
      flag_super = ((std::fabs(block[0] - 2) > EPSBLOCK2) && (std::fabs(block[1] - 2) > EPSBLOCK2));
      bonus[j].flag_super = flag_super;
      bonus[j].ilocal = i;
      bonus[j].radcirc = compute_radcirc(shape, block, flag_super);
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
    double *block = bonus[j].block;
    buf[m++] = shape[0];
    buf[m++] = shape[1];
    buf[m++] = shape[2];
    buf[m++] = quat[0];
    buf[m++] = quat[1];
    buf[m++] = quat[2];
    buf[m++] = quat[3];
    buf[m++] = block[0];
    buf[m++] = block[1];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecEllipsoid::unpack_exchange_bonus(int ilocal, double *buf)
{
  int m = 0;

  if (ubuf(buf[m++]).i == 0)
    ellipsoid[ilocal] = -1;
  else {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    double *shape = bonus[nlocal_bonus].shape;
    double *quat = bonus[nlocal_bonus].quat;
    double *block = bonus[nlocal_bonus].block;
    bool &flag_super = bonus[nlocal_bonus].flag_super;
    shape[0] = buf[m++];
    shape[1] = buf[m++];
    shape[2] = buf[m++];
    quat[0] = buf[m++];
    quat[1] = buf[m++];
    quat[2] = buf[m++];
    quat[3] = buf[m++];
    block[0] = buf[m++];
    block[1] = buf[m++];
    block[2] = block[0] / block[1];
    flag_super = ((std::fabs(block[0] - 2) > EPSBLOCK2) && (std::fabs(block[1] - 2) > EPSBLOCK2));
    bonus[nlocal_bonus].radcirc = compute_radcirc(shape, block, flag_super);
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
    buf[m++] = bonus[j].block[0];
    buf[m++] = bonus[j].block[1];
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
    double *block = bonus[nlocal_bonus].block;
    bool &flag_super = bonus[nlocal_bonus].flag_super;
    shape[0] = buf[m++];
    shape[1] = buf[m++];
    shape[2] = buf[m++];
    quat[0] = buf[m++];
    quat[1] = buf[m++];
    quat[2] = buf[m++];
    quat[3] = buf[m++];
    block[0] = buf[m++];
    block[1] = buf[m++];
    block[2] = block[0] / block[1];
    flag_super = ((std::fabs(block[0] - 2) > EPSBLOCK2) && (std::fabs(block[1] - 2) > EPSBLOCK2));
    bonus[nlocal_bonus].radcirc = compute_radcirc(shape, block, flag_super);
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

  // Blockiness exponents can be given optionally for superellipsoids

  double *block = bonus[nlocal_bonus].block;
  bool &flag_super = bonus[nlocal_bonus].flag_super;
  if (ivalue == values.size()) {
    block[0] = block[1] = 2.0;
    flag_super = false;
  }
  else {
    block[0] = utils::numeric(FLERR, values[ivalue++], true, lmp);
    block[1] = utils::numeric(FLERR, values[ivalue++], true, lmp);
    flag_super = ((std::fabs(block[0] - 2) > EPSBLOCK2) && (std::fabs(block[1] - 2) > EPSBLOCK2));
  }
  block[2] = block[0] / block[1]; // ASSUMES EVEN NUMBERS ONLY?

  // reset ellipsoid mass
  // previously stored density in rmass

  rmass[m] *= compute_volume(shape, block, flag_super);

  bonus[nlocal_bonus].radcirc = compute_radcirc(shape, block, flag_super);
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
  double *shape, *block;
  bool flag_super;

  ellipsoid_flag = atom->ellipsoid[ilocal];
  rmass_one = atom->rmass[ilocal];

  if (ellipsoid_flag < 0)
    ellipsoid[ilocal] = 0;
  else
    ellipsoid[ilocal] = 1;

  if (ellipsoid_flag >= 0) {
    shape = bonus[ellipsoid_flag].shape;
    block = bonus[ellipsoid_flag].block;
    flag_super = bonus[ellipsoid_flag].flag_super;
    rmass[ilocal] /= compute_volume(shape, block, flag_super);
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
      buf[m++] = bonus[j].block[0];
      buf[m++] = bonus[j].block[1];
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
    fmt::print(fp, "{} {} {} {} {} {} {} {} {} {}\n", ubuf(buf[i]).i, buf[i + 1], buf[i + 2], buf[i + 3],
               buf[i + 4], buf[i + 5], buf[i + 6], buf[i + 7], buf[i + 8], buf[i + 9]);
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
    double *block = bonus[nlocal_bonus].block;
    bool flag_super = bonus[nlocal_bonus].flag_super;
    shape[0] = shapex;
    shape[1] = shapey;
    shape[2] = shapez;
    quat[0] = 1.0;
    quat[1] = 0.0;
    quat[2] = 0.0;
    quat[3] = 0.0;
    bonus[nlocal_bonus].radcirc = compute_radcirc(shape, block, flag_super);
    bonus[nlocal_bonus].ilocal = i;
    ellipsoid[i] = nlocal_bonus++;
  } else if (shapex == 0.0 && shapey == 0.0 && shapez == 0.0) {
    copy_bonus_all(nlocal_bonus - 1, ellipsoid[i]);
    nlocal_bonus--;
    ellipsoid[i] = -1;
  } else {
    double *shape = bonus[ellipsoid[i]].shape;
    double *block = bonus[ellipsoid[i]].block;
    bool flag_super = bonus[ellipsoid[i]].flag_super;
    shape[0] = shapex;
    shape[1] = shapey;
    shape[2] = shapez;
    bonus[ellipsoid[i]].radcirc = compute_radcirc(shape, block, flag_super);
  }
}

/* ----------------------------------------------------------------------
   set block values in bonus data for particle I
   oriented aligned with xyz axes
   this may create entry in bonus data
------------------------------------------------------------------------- */

void AtomVecEllipsoid::set_block(int i, double blockn1, double blockn2)
{
  if (ellipsoid[i] < 0) {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    double *shape = bonus[nlocal_bonus].shape;
    double *quat = bonus[nlocal_bonus].quat;
    double *block = bonus[nlocal_bonus].block;
    bool &flag_super = bonus[nlocal_bonus].flag_super;
    block[0] = blockn1;
    block[1] = blockn2;
    block[2] = blockn1 / blockn2;
    quat[0] = 1.0;
    quat[1] = 0.0;
    quat[2] = 0.0;
    quat[3] = 0.0;
    bonus[nlocal_bonus].ilocal = i;
    flag_super = ((std::fabs(blockn1 - 2) > EPSBLOCK2) && (std::fabs(blockn2 - 2) > EPSBLOCK2));
    bonus[nlocal_bonus].radcirc = compute_radcirc(shape, block, flag_super);
    ellipsoid[i] = nlocal_bonus++;
  } else {
    double *shape = bonus[ellipsoid[i]].shape;
    double *block = bonus[ellipsoid[i]].block;
    bool &flag_super = bonus[ellipsoid[i]].flag_super;
    block[0] = blockn1;
    block[1] = blockn2;
    block[2] = blockn1 / blockn2;
    flag_super = ((std::fabs(blockn1 - 2) > EPSBLOCK2) && (std::fabs(blockn2 - 2) > EPSBLOCK2));
    bonus[ellipsoid[i]].radcirc = compute_radcirc(shape, block, flag_super);
  }
}

/* ----------------------------------------------------------------------
   compute the circumscribed radius to the ellipsoid
------------------------------------------------------------------------- */

double AtomVecEllipsoid::compute_radcirc(double *shape, double *block, bool flag_super)
{
  if (!flag_super) return std::max(std::max(shape[0], shape[1]), shape[2]);

  // Super ellipsoid
  double a = shape[0], b = shape[1], c = shape[2];
  double n1 = block[0], n2 = block[1], n1divn2 = block[2];
  if (shape[0] < shape[1]) {a = shape[1]; b = shape[0];}

  // Cylinder approximation for n2=2

  if (std::fabs(n2 - 2.0) < EPSBLOCK2) return sqrt(a * a + c * c);

  // Ellipsoid approximation for n1=2

  if (std::fabs(n1 - 2.0) < EPSBLOCK2) return std::max(c, sqrt(a * a + b * b));

  // Bounding box approximation when n1>2 and n2>2

  return sqrt(a * a + b * b + c * c);

  // General super-ellipsoid, Eq. (12) of Podlozhnyuk et al. 2017
  // Not sure if exact solution worth it compared to boundig box diagonal
  // If both blockiness exponents are greater than 2, the exact radius does not
  // seem significantly smaller than the bounding box diagonal. At most sqrt(3)~ 70% too large
  /*
  double x, y, z, alpha, beta, gamma, xtilde;
  double small = 0.1; // TO AVOID OVERFLOW IN POW

  alpha = std::fabs(n2 - 2.0) > small ? std::pow(b / a, 2.0 / (n2 - 2.0)) : 0.0;
  gamma = std::fabs(n1divn2 - 1.0) > small ? std::pow((1.0 + std::pow(alpha, n2)), n1divn2 - 1.0) : 1.0;
  beta = std::pow(gamma * c * c / (a * a), 1.0 / std::max(n1 - 2.0, small));
  xtilde = 1.0 / std::pow(std::pow(1.0 + std::pow(alpha, n2), n1divn2) + std::pow(beta, n1), 1.0 / n1);
  x = a * xtilde;
  y = alpha * b * xtilde;
  z = beta * c * xtilde;
  return sqrt(x * x + y * y + z * z);
  */
}

/* ----------------------------------------------------------------------
   compute the volume of the ellipsoid
------------------------------------------------------------------------- */

double AtomVecEllipsoid::compute_volume(double *shape, double *block, bool flag_super)
{
  double unitvol = 4.0 * MY_PI / 3.0;

  // super-ellipsoid, Eq. (12) of Jaklic and Solina, 2003, for p = q = r = 0

  if (flag_super) {
    double e1 = 2.0 / block[0], e2 = 2.0 / block[1];
    unitvol = e1 * e2 * std::tr1::beta(0.5 * e1, 1.0 + e1) *
                        std::tr1::beta(0.5 * e2, 0.5 * e2); // CAN'T GET std::beta to be recognized, need help
  }
  return unitvol * shape[0] * shape[1] * shape[2];
}
