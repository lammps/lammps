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

#include "atom_vec_line.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using MathConst::MY_PI;

static constexpr double EPSILON = 0.001;

/* ---------------------------------------------------------------------- */

AtomVecLine::AtomVecLine(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::ATOMIC;
  bonus_flag = 1;

  size_forward_bonus = 1;
  size_border_bonus = 3;
  size_restart_bonus_one = 3;
  size_data_bonus = 5;

  atom->line_flag = 1;
  atom->molecule_flag = atom->rmass_flag = 1;
  atom->radius_flag = atom->omega_flag = atom->torque_flag = 1;

  nlocal_bonus = nghost_bonus = nmax_bonus = 0;
  bonus = nullptr;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = {"molecule", "radius", "rmass", "omega", "torque", "line"};
  fields_copy = {"molecule", "radius", "rmass", "omega"};
  fields_comm_vel = {"omega"};
  fields_reverse = {"torque"};
  fields_border = {"molecule", "radius", "rmass"};
  fields_border_vel = {"molecule", "radius", "rmass", "omega"};
  fields_exchange = {"molecule", "radius", "rmass", "omega"};
  fields_restart = {"molecule", "radius", "rmass", "omega"};
  fields_create = {"molecule", "radius", "rmass", "omega", "line"};
  fields_data_atom = {"id", "molecule", "type", "line", "rmass", "x"};
  fields_data_vel = {"id", "v", "omega"};

  setup_fields();
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
    error->all(FLERR, "Atom_style line can only be used in 2d simulations");
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecLine::grow_pointers()
{
  line = atom->line;
  radius = atom->radius;
  rmass = atom->rmass;
  omega = atom->omega;
}

/* ----------------------------------------------------------------------
   grow bonus data structure
------------------------------------------------------------------------- */

void AtomVecLine::grow_bonus()
{
  nmax_bonus = grow_nmax_bonus(nmax_bonus);
  if (nmax_bonus < 0) error->one(FLERR, "Per-processor system is too big");

  bonus = (Bonus *) memory->srealloc(bonus, nmax_bonus * sizeof(Bonus), "atom:bonus");
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecLine::copy_bonus(int i, int j, int delflag)
{
  // if deleting atom J via delflag and J has bonus data, then delete it

  if (delflag && line[j] >= 0) {
    copy_bonus_all(nlocal_bonus - 1, line[j]);
    nlocal_bonus--;
  }

  // if atom I has bonus data, reset I's bonus.ilocal to loc J
  // do NOT do this if self-copy (I=J) since I's bonus data is already deleted

  if (line[i] >= 0 && i != j) bonus[line[i]].ilocal = j;
  line[j] = line[i];
}

/* ----------------------------------------------------------------------
   copy bonus data from I to J, effectively deleting the J entry
   also reset line that points to I to now point to J
------------------------------------------------------------------------- */

void AtomVecLine::copy_bonus_all(int i, int j)
{
  line[bonus[i].ilocal] = j;
  memcpy(&bonus[j], &bonus[i], sizeof(Bonus));
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

/* ---------------------------------------------------------------------- */

int AtomVecLine::pack_comm_bonus(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (line[j] >= 0) buf[m++] = bonus[line[j]].theta;
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecLine::unpack_comm_bonus(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (line[i] >= 0) bonus[line[i]].theta = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecLine::pack_border_bonus(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (line[j] < 0)
      buf[m++] = ubuf(0).d;
    else {
      buf[m++] = ubuf(1).d;
      buf[m++] = bonus[line[j]].length;
      buf[m++] = bonus[line[j]].theta;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecLine::unpack_border_bonus(int n, int first, double *buf)
{
  int i, j, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    line[i] = (int) ubuf(buf[m++]).i;
    if (line[i] == 0)
      line[i] = -1;
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

int AtomVecLine::pack_exchange_bonus(int i, double *buf)
{
  int m = 0;

  if (line[i] < 0)
    buf[m++] = ubuf(0).d;
  else {
    buf[m++] = ubuf(1).d;
    int j = line[i];
    buf[m++] = bonus[j].length;
    buf[m++] = bonus[j].theta;
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecLine::unpack_exchange_bonus(int ilocal, double *buf)
{
  int m = 0;

  line[ilocal] = (int) ubuf(buf[m++]).i;
  if (line[ilocal] == 0)
    line[ilocal] = -1;
  else {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    bonus[nlocal_bonus].length = buf[m++];
    bonus[nlocal_bonus].theta = buf[m++];
    bonus[nlocal_bonus].ilocal = ilocal;
    line[ilocal] = nlocal_bonus++;
  }

  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecLine::size_restart_bonus()
{
  int i;

  int n = 0;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) {
    if (line[i] >= 0)
      n += size_restart_bonus_one;
    else
      n++;
  }

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVecLine::pack_restart_bonus(int i, double *buf)
{
  int m = 0;

  if (line[i] < 0)
    buf[m++] = ubuf(0).d;
  else {
    buf[m++] = ubuf(1).d;
    int j = line[i];
    buf[m++] = bonus[j].length;
    buf[m++] = bonus[j].theta;
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including bonus data
------------------------------------------------------------------------- */

int AtomVecLine::unpack_restart_bonus(int ilocal, double *buf)
{
  int m = 0;

  line[ilocal] = (int) ubuf(buf[m++]).i;
  if (line[ilocal] == 0)
    line[ilocal] = -1;
  else {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    bonus[nlocal_bonus].length = buf[m++];
    bonus[nlocal_bonus].theta = buf[m++];
    bonus[nlocal_bonus].ilocal = ilocal;
    line[ilocal] = nlocal_bonus++;
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack one line from Lines section of data file
------------------------------------------------------------------------- */

void AtomVecLine::data_atom_bonus(int m, const std::vector<std::string> &values)
{
  if (line[m]) error->one(FLERR, "Assigning line parameters to non-line atom");

  if (nlocal_bonus == nmax_bonus) grow_bonus();

  int ivalue = 1;
  double x1 = utils::numeric(FLERR, values[ivalue++], true, lmp);
  double y1 = utils::numeric(FLERR, values[ivalue++], true, lmp);
  double x2 = utils::numeric(FLERR, values[ivalue++], true, lmp);
  double y2 = utils::numeric(FLERR, values[ivalue++], true, lmp);

  // convert x1/y1 and x2/y2 from general to restricted triclniic
  // x is already restricted triclinic

  double coords[3];

  if (domain->triclinic_general) {
    coords[0] = x1; coords[1] = y1; coords[2] = 0.0;
    domain->general_to_restricted_coords(coords);
    x1 = coords[0]; y1 = coords[1];
    coords[0] = x2; coords[1] = y2; coords[2] = 0.0;
    domain->general_to_restricted_coords(coords);
    x2 = coords[0]; y2 = coords[1];
  }

  // remap end points to be near x
  // necessary if atom x was remapped into periodic box

  coords[0] = x1; coords[1] = y1; coords[2] = 0.0;
  domain->remap_near(coords,x[m]);
  x1 = coords[0]; y1 = coords[1];
  coords[0] = x2; coords[1] = y2; coords[2] = 0.0;
  domain->remap_near(coords,x[m]);
  x2 = coords[0]; y2 = coords[1];

  // calculate length and theta
  // error if segment center is not within EPSILON of atom x
  // reset atom x to center point

  double dx = x2 - x1;
  double dy = y2 - y1;
  double length = sqrt(dx * dx + dy * dy);

  bonus[nlocal_bonus].length = length;
  if (dy >= 0.0)
    bonus[nlocal_bonus].theta = acos(dx / length);
  else
    bonus[nlocal_bonus].theta = -acos(dx / length);

  double xc = 0.5 * (x1 + x2);
  double yc = 0.5 * (y1 + y2);
  dx = xc - x[m][0];
  dy = yc - x[m][1];
  double delta = sqrt(dx * dx + dy * dy);

  if (delta / length > EPSILON) error->one(FLERR, "Inconsistent line segment in data file");

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
   return # of bytes of allocated bonus memory
------------------------------------------------------------------------- */

double AtomVecLine::memory_usage_bonus()
{
  double bytes = 0;
  bytes += (double) nmax_bonus * sizeof(Bonus);
  return bytes;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecLine::create_atom_post(int ilocal)
{
  double radius_one = 0.5;
  radius[ilocal] = radius_one;
  rmass[ilocal] = 4.0 * MY_PI / 3.0 * radius_one * radius_one * radius_one;
  line[ilocal] = -1;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecLine::data_atom_post(int ilocal)
{
  line_flag = line[ilocal];
  if (line_flag == 0)
    line_flag = -1;
  else if (line_flag == 1)
    line_flag = 0;
  else
    error->one(FLERR, "Invalid line flag in Atoms section of data file");
  line[ilocal] = line_flag;

  if (rmass[ilocal] <= 0.0) error->one(FLERR, "Invalid density in Atoms section of data file");

  if (line_flag < 0) {
    double radius_one = 0.5;
    radius[ilocal] = radius_one;
    rmass[ilocal] *= 4.0 * MY_PI / 3.0 * radius_one * radius_one * radius_one;
  } else
    radius[ilocal] = 0.0;

  omega[ilocal][0] = 0.0;
  omega[ilocal][1] = 0.0;
  omega[ilocal][2] = 0.0;
}

/* ----------------------------------------------------------------------
   modify values for AtomVec::pack_data() to pack
------------------------------------------------------------------------- */

void AtomVecLine::pack_data_pre(int ilocal)
{
  line_flag = line[ilocal];
  rmass_one = rmass[ilocal];

  if (line_flag < 0)
    line[ilocal] = 0;
  else
    line[ilocal] = 1;

  if (line_flag < 0) {
    double radius_one = radius[ilocal];
    rmass[ilocal] /= 4.0 * MY_PI / 3.0 * radius_one * radius_one * radius_one;
  } else
    rmass[ilocal] /= bonus[line_flag].length;
}

/* ----------------------------------------------------------------------
   unmodify values packed by AtomVec::pack_data()
------------------------------------------------------------------------- */

void AtomVecLine::pack_data_post(int ilocal)
{
  line[ilocal] = line_flag;
  rmass[ilocal] = rmass_one;
}

/* ----------------------------------------------------------------------
   pack bonus line info for writing to data file
   if buf is nullptr, just return buffer size
------------------------------------------------------------------------- */

int AtomVecLine::pack_data_bonus(double *buf, int /*flag*/)
{
  int i, j;
  double length, theta;
  double xc, yc, x1, x2, y1, y2;
  double coords[3];

  int triclinic_general = domain->triclinic_general;

  double **x_bonus;
  if (triclinic_general) x_bonus = x_hold;
  else x_bonus = x;

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int m = 0;
  for (i = 0; i < nlocal; i++) {
    if (line[i] < 0) continue;
    if (buf) {
      buf[m++] = ubuf(tag[i]).d;
      j = line[i];
      length = bonus[j].length;
      theta = bonus[j].theta;

      xc = x_bonus[i][0];
      yc = x_bonus[i][1];
      x1 = xc - 0.5 * cos(theta) * length;
      y1 = yc - 0.5 * sin(theta) * length;
      x2 = xc + 0.5 * cos(theta) * length;
      y2 = yc + 0.5 * sin(theta) * length;
      buf[m++] = x1;
      buf[m++] = y1;
      buf[m++] = x2;
      buf[m++] = y2;

      // if triclinic_general:
      // rotate 4 buf values from restricted to general triclinic
      // output by write_data_bonus() as x1/y1 and x2/y2

      if (triclinic_general) {
        coords[0] = buf[m-4]; coords[1] = buf[m-3]; coords[2] = 0.0;
        domain->restricted_to_general_coords(coords);
        buf[m-4] = coords[0]; buf[m-3] = coords[1];
        coords[0] = buf[m-2]; coords[1] = buf[m-1]; coords[2] = 0.0;
        domain->restricted_to_general_coords(coords);
        buf[m-2] = coords[0]; buf[m-1] = coords[1];
      }

    } else
      m += size_data_bonus;
  }
  return m;
}

/* ----------------------------------------------------------------------
   write bonus line info to data file
------------------------------------------------------------------------- */

void AtomVecLine::write_data_bonus(FILE *fp, int n, double *buf, int /*flag*/)
{
  int i = 0;
  while (i < n) {
    fmt::print(fp, "{} {} {} {} {}\n", ubuf(buf[i]).i, buf[i + 1], buf[i + 2], buf[i + 3],
               buf[i + 4]);
    i += size_data_bonus;
  }
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
    copy_bonus_all(nlocal_bonus - 1, line[i]);
    nlocal_bonus--;
    line[i] = -1;
  } else
    bonus[line[i]].length = value;

  // also set radius = half of length
  // unless value = 0.0, then set diameter = 1.0

  radius[i] = 0.5 * value;
  if (value == 0.0) radius[i] = 0.5;
}
